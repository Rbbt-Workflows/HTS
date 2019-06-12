require 'rbbt/sources/organism'

module HTS

  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil
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

    gatk("FastqToSam", args)
  end

  dep :uBAM, :compute => :produce
  extension :ubam
  task :mark_adapters => :binary do
    args = {}
    args["INPUT"] = step(:uBAM).path
    args["METRICS"] = file('metrics.txt')
    args["OUTPUT"] = self.tmp_path

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    gatk("MarkIlluminaAdapters", args)
  end

  dep :mark_adapters
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :bwa_mem_args, :string, "Arg string", "-M -p"
  extension :bam
  task :BAM_bwa => :binary do |reference, bwa_mem_args|

    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    Open.rm file('SamToFastq')
    Misc.with_fifo(file('SamToFastq.fastq')) do |s2f_path|
      Misc.with_fifo(file('bwa.bam')) do |bwa_bam|

        args = {}
        args["INPUT"] = step(:mark_adapters).path
        args["FASTQ"] = s2f_path
        args["CLIPPING_ATTRIBUTE"] = "XT"
        args["CLIPPING_ACTION"] = "2"
        args["INTERLEAVE"] = "true"
        #args["NON_PF"] = "true"

        io_s2f = gatk_io("SamToFastq", args)
        t_s2f = Thread.new do
          while line = io_s2f.gets
            Log.debug line
          end
        end

        bwa_mem_args += " -t " << config('cpus', 'bwa', :default => 8) 
        io_bwa = BWA.mem([s2f_path], reference, bwa_mem_args)
        Misc.consume_stream io_bwa, true, bwa_bam

        uBAM = step('uBAM').path

        args = {}
        args["ALIGNED_BAM"] = bwa_bam
        args["UNMAPPED_BAM"] = uBAM
        args["OUTPUT"] = self.tmp_path
        args["R"] = reference
        args["CREATE_INDEX"] = "false"
        args["ADD_MATE_CIGAR"] = "true"
        args["CLIP_ADAPTERS"] = "false"
        args["CLIP_OVERLAPPING_READS"] = "true"
        args["INCLUDE_SECONDARY_ALIGNMENTS"] = "true"
        args["MAX_INSERTIONS_OR_DELETIONS"] = "-1"
        args["PRIMARY_ALIGNMENT_STRATEGY"] = "MostDistant"
        args["ATTRIBUTES_TO_RETAIN"] = "XS"
        args["SORT_ORDER"] = "queryname"

        gatk("MergeBamAlignment", args)

        FileUtils.rm_rf bwa_bam
      end
      FileUtils.rm_rf s2f_path
    end
    nil
  end

  dep :BAM_bwa
  extension :bam
  task :BAM_duplicates => :binary do
    Open.mkdir files_dir 
    output = file('out.bam')

    args = {}
    args["INPUT"] = step(:BAM_bwa).path
    args["OUTPUT"] = output
    args["METRICS_FILE"] = file('metrics.txt') 
    args["ASSUME_SORT_ORDER"] = 'queryname'
    args["CREATE_INDEX"] = 'false'

    gatk("MarkDuplicates", args)

    Open.mv output, self.path
    nil
  end

  dep :BAM_duplicates
  extension :bam
  task :BAM_sorted => :binary do
    Open.mkdir files_dir 
    sorted = file('sorted.bam')
    
    args = {}
    args["INPUT"] = step(:BAM_duplicates).path
    args["OUTPUT"] = sorted
    args["SORT_ORDER"] = 'coordinate'
    gatk("SortSam", args)
    Open.mv sorted, self.path
    nil
  end

  dep :BAM_sorted
  extension :bam
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_rescore => :binary do |interval_list|

    interval_list = nil if interval_list == "none"

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference

    args = {}
    args["reference"] = reference
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list
    args["output"] = file('recal_data.table')

    known_sites = [] 
    ["miller_indels", "dbsnp", "1000G_indels"].each do |file|
      vcf = vcf_file reference, file
      vcf = GATK.prepare_VCF_AF_only vcf
      known_sites << vcf
    end
    args["known-sites"] = known_sites

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    #{{{ Recalibration
    shard = config('shard', :gatk, :rescore, :baserecalibrator, :BaseRecalibrator)
    if shard == 'true'
      contigs = Samtools.reference_contigs reference
      bam_file = Samtools.prepare_BAM(step(:BAM_sorted))
      args["input"] = bam_file

      cpus = config('cpus', :shard, :rescore, :baserecalibrator, :BaseRecalibrator)
      args["intervals"] ||= nil
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing BaseRecalibrator sharded")

      outfiles = Path.setup(file('outfiles_recall'))
      GATKShard.cmd("BaseRecalibrator", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
        bar.tick
        Open.mv ioutfile, outfiles[File.basename(ioutfile + '.report')]
        nil
      end
      bar.remove 

      args = {}
      args["I"] = outfiles.glob("*.report")
      args["O"] = file('recal_data.table')
      gatk("GatherBQSRReports", args)
      Open.rm_rf outfiles
      nil
    else
      bam_file = interval_list ? Samtools.prepare_BAM(step(:BAM_sorted)) : step(:BAM_sorted).path

      args["input"] = bam_file
      gatk("BaseRecalibrator", args)
    end

    #{{{ Apply
    output = file('out.bam')
    args = {}
    args["input"] = bam_file
    args["output"] = output
    args["bqsr-recal-file"] = file('recal_data.table')
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list

    shard = config('shard', :gatk, :rescore, :apply_rescore, :apply_bqsr, :ApplyBQSR)

    if shard == 'true'
      contigs = Samtools.reference_contigs reference
      args["input"] = bam_file

      cpus = config('cpus', :shard, :rescore, :apply_rescore, :apply_bqsr, :ApplyBQSR)
      args["intervals"] ||= nil
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      intervals = (interval_list || intervals_for_reference(reference))

      outfiles = Path.setup(file('outfiles'))
      Open.mkdir outfiles

      bar = self.progress_bar("Processing ApplyBQSR sharded")
      GATKShard.cmd("ApplyBQSR", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
        bar.tick
        chr, pos = Samtools.BAM_start ioutfile
        target = outfiles[File.basename(ioutfile) + '.bam']
        Open.mv ioutfile, target if ! (chr.nil? or chr.empty?)
        nil
      end
      bar.remove 

      contigs = Samtools.reference_contigs reference
      sorted_parts = outfiles.glob("*.bam").sort{|a,b| Misc.genomic_location_cmp_contigs(File.basename(a).split(",").first, File.basename(b).split(",").first, contigs, '__')}

      args = {}
      args["I"] = sorted_parts
      args["O"] = output
      args["CREATE_INDEX"] = 'false'
      args["CREATE_MD5_FILE"] = 'false'
      gatk("GatherBamFiles", args)

      Open.rm_rf outfiles
    else
      bam_file = interval_list ? Samtools.prepare_BAM(step(:BAM_sorted)) : step(:BAM_sorted).path

      args["input"] = bam_file
      gatk("ApplyBQSR", args)
    end

    FileUtils.mv output, self.path
    nil
  end

  extension :bam
  dep_task :BAM, HTS, :BAM_rescore
end

require 'HTS/tasks/BAM/plumbing'
require 'HTS/tasks/BAM/other_aligners'
require 'HTS/tasks/BAM/util'
require 'HTS/tasks/BAM/filter'
