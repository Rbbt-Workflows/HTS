require 'rbbt/sources/organism'

module HTS

  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
  input :library_name, :string, "LIBRARY_NAME BAM field", "DefaultLibraryName"
  input :platform_unit, :string, "PLATFORM_UNIT BAM field", "DefaultPlatformUnit"
  input :platform, :string, "PLATFORM BAM field", "DefaultPlatform"
  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field", "DefaultSequencingCenter"
  extension :ubam
  task :uBAM => :binary do |fastq1, fastq2, read_group_name, sample_name, library_name, platform_unit, platform, sequencing_center, run_date|
    sample_name ||= self.clean_name
    read_group_name ||= sample_name

    args = {}
    args["FASTQ"] = fastq1
    args["FASTQ2"] = fastq2
    args["OUTPUT"] = self.tmp_path
    args["READ_GROUP_NAME"] = read_group_name
    args["SAMPLE_NAME"] = sample_name
    args["LIBRARY_NAME"] = library_name
    args["PLATFORM_UNIT"] = platform_unit
    args["PLATFORM"] = platform
    args["SEQUENCING_CENTER"] = sequencing_center
    args["RUN_DATE"] = run_date

    GATK.run_log("FastqToSam", args)
  end

  dep :uBAM, :compute => :produce
  extension :ubam
  task :mark_adapters => :binary do
    args = {}
    args["INPUT"] = step(:uBAM).path
    args["METRICS"] = file('metrics.txt')
    args["OUTPUT"] = self.tmp_path

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkIlluminaAdapters", args)
  end

  dep :mark_adapters
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :bwa_mem_args, :string, "Arg string", "-M -p"
  extension :bam
  task :BAM => :binary do |reference, bwa_mem_args|

    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    Open.rm file('SamToFastq')
    Misc.with_fifo(file('SamToFastq')) do |s2f_path|

      args = {}
      args["INPUT"] = step(:mark_adapters).path
      args["FASTQ"] = s2f_path
      args["CLIPPING_ATTRIBUTE"] = "XT"
      args["CLIPPING_ACTION"] = "2"
      args["INTERLEAVE"] = "true"
      args["NON_PF"] = "true"

      io_s2f = GATK.run("SamToFastq", args)
      t_s2f = Thread.new do
        while line = io_s2f.gets
          Log.debug line
        end
      end

      bwa_mem_args += " -t " << config('cpus', 'bwa', :default => 8) 
      io_bwa = BWA.mem([s2f_path], reference, bwa_mem_args)

      args = {}
      args["ALIGNED_BAM"] = "/dev/stdin"
      args["UNMAPPED_BAM"] = step('uBAM').path
      args["OUTPUT"] = self.tmp_path
      args["R"] = reference
      args["CREATE_INDEX"] = "true"
      args["ADD_MATE_CIGAR"] = "true"
      args["CLIP_ADAPTERS"] = "false"
      args["CLIP_OVERLAPPING_READS"] = "true"
      args["INCLUDE_SECONDARY_ALIGNMENTS"] = "true"
      args["MAX_INSERTIONS_OR_DELETIONS"] = "-1"
      args["PRIMARY_ALIGNMENT_STRATEGY"] = "MostDistant"
      args["ATTRIBUTES_TO_RETAIN"] = "XS"
      GATK.run_log("MergeBamAlignment", args, io_bwa)

      FileUtils.rm_rf s2f_path
    end
    nil
  end

  dep :BAM
  extension :bam
  task :BAM_duplicates => :binary do
    args = {}
    args["I"] = step(:BAM).path
    args["O"] = self.tmp_path
    args["M"] = file('metrics.txt')

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkDuplicates", args)
  end


  dep :BAM_duplicates
  extension :bam
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_rescore => :binary do |interval_list|

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference
    reference_code = self.recursive_inputs[:reference]

    DbSNP["All.vcf.gz.tbi"].produce.find
    db_SNP = DbSNP["All.vcf.gz"].produce.find

    bam_file = interval_list ? Samtools.prepare_BAM(step(:BAM_duplicates)) : step(:BAM_duplicates).path

    args = {}
    args["input"] = bam_file
    args["reference"] = reference
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    args["output"] = file('recal_data.table')

    known_sites = [] 
    ["dbsnp_138.vcf.gz","Miller_1000G_indels.vcf.gz", "1000G_phase1.indels.vcf.gz"].each do |file|
      known_sites << GATK.known_sites[reference_code][file].produce.find
    end
    args["known-sites"] = known_sites

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("BaseRecalibrator", args)

    args = {}
    args["input"] = bam_file
    args["output"] = self.tmp_path
    args["bqsr-recal-file"] = file('recal_data.table')
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    GATK.run_log("ApplyBQSR", args)
  end
end

require 'HTS/tasks/BAM/plumbing'
require 'HTS/tasks/BAM/other_aligners'
require 'HTS/tasks/BAM/util'
