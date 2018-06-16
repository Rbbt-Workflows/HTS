require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/NGS'

require 'tools/BWA'
require 'tools/GATK'
require 'tools/samtools'

module HTS
  extend Workflow
  
  helper :reference_file do |reference|
    case reference
    when 'hg19'
      BWA.references.hg19["reference.fa"].produce.find
    when 'hg38'
      BWA.references.hg38["reference.fa"].produce.find
    when 'b37'
      BWA.references.b37["reference.fa"].produce.find
    else
      reference
    end
  end

  input :fastq1, :file, "FASTQ 1 file"
  input :fastq2, :file, "FASTQ 2 file"
  input :read_group_name, :string, "READ_GROUP_NAME BAM field"
  input :sample_name, :string, "SAMPLE_NAME BAM field"
  input :library_name, :string, "LIBRARY_NAME BAM field"
  input :platform_unit, :string, "PLATFORM_UNIT BAM field"
  input :platform, :string, "PLATFORM BAM field"
  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field"
  extension :bam
  task :uBAM => :binary do |fastq1, fastq2, read_group_name, sample_name, library_name, platform_unit, platform, sequencing_center, run_date|
    args = {}
    args["FASTQ"] = fastq1
    args["FASTQ2"] = fastq2
    args["OUTPUT"] = self.path
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
  extension :bam
  task :mark_adapters => :binary do
    args = {}
    args["INPUT"] = step(:uBAM).path
    args["METRICS"] = file('metrics.txt')
    args["OUTPUT"] = self.path

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkIlluminaAdapters", args)
  end

  dep :mark_adapters
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38)
  input :bwa_mem_args, :string, "Arg string", "-M -p"
  extension :bam
  task :BAM => :binary do |reference, bwa_mem_args|

    reference = reference_file(reference)

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
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

      io_bwa = BWA.mem([s2f_path], reference, bwa_mem_args)

      args = {}
      args["ALIGNED_BAM"] = "/dev/stdin"
      args["UNMAPPED_BAM"] = step('uBAM').path
      args["OUTPUT"] = self.path
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
    end
    nil
  end

  dep :BAM
  extension :bam
  task :BAM_duplicates => :binary do
    args = {}
    args["INPUT"] = step(:BAM).path
    args["OUTPUT"] = self.path
    args["METRICS_FILE"] = file('metrics.txt')

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkDuplicates", args)
  end

  dep :BAM_duplicates
  extension :bam
  task :BAM_rescore => :binary do

    reference = reference_file self.recursive_inputs[:reference]
    reference_code = self.recursive_inputs[:reference]

    DbSNP["All.vcf.gz.tbi"].produce.find
    db_SNP = DbSNP["All.vcf.gz"].produce.find

    args = {}
    args["input"] = step(:BAM_duplicates).path
    args["reference"] = reference
    args["output"] = file('recal_data.table')

    known_sites = [DbSNP["All.vcf.gz"].produce.find] 
    ["Miller_1000G_indels.vcf.gz", "1000G_phase1.indels.vcf.gz"].each do |file|
      known_sites << GATK.known_sites[reference_code][file].produce.find
    end
    args["known-sites"] = known_sites

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("BaseRecalibrator", args)

    args = {}
    args["input"] = step(:BAM_duplicates).path
    args["output"] = self.path
    args["bqsr-recal-file"] = file('recal_data.table')
    GATK.run_log("ApplyBQSR", args)
  end

  input :tumor, :file, "Tumor BAM"
  input :normal, :file, "Tumor BAM"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38)
  extension :vcf
  task :mutect2 => :text do |tumor,normal,reference|

    reference ||= self.recursive_inputs[:reference]

    reference = reference_file reference
    reference_code = self.recursive_inputs[:reference]

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    #tumor_sample = CMD.cmd("#{Samtools::Samtools_CMD} view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #tumor_sample = CMD.cmd("samtools view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #normal_sample = CMD.cmd("samtools view -H '#{normal}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal)

    args["input"] = [tumor, normal]
    args["output"] = self.path
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample
    args["bam-output"] = file('haplotype.bam')

    GATK.run_log("Mutect2", args)
  end

  dep :mutect2
  extension :vcf
  task :mutect2_filtered => :tsv do
    args = {}
    FileUtils.mkdir_p files_dir

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2).path, tmp

    args["variant"] = tmp
    args["output"] = self.path
    GATK.run_log("FilterMutectCalls", args)
  end

  input :fastq1, :file, "FASTQ 1 file"
  input :fastq2, :file, "FASTQ 2 file"
  task :HLA => :text do |fastq1,fastq2|
    fastqs = []
    [fastq1, fastq2].compact.each do |fastq|
      tmp_bam = TmpFile.tmp_file + '.bam'
      CMD.cmd("#{Rbbt.software.opt.seqan.bin.razers3.find} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{Rbbt.software.opt.OptiType.data["hla_reference_dna.fasta"].find}' '#{fastq}'")
      tmp_fastq = TmpFile.tmp_file + '.fastq'
      CMD.cmd("samtools bam2fq '#{tmp_bam}' > '#{tmp_fastq}'")
      CMD.cmd("rm '#{tmp_bam}'")
      fastqs << tmp_fastq
    end

    output = files_dir
    CMD.cmd("python #{Rbbt.software.opt.OptiType.produce["OptiTypePipeline.py"].find} -i #{fastqs.collect{|f| "'#{f}'"} * " "} --dna -v -o '#{output}' ")

    res = Dir.glob(files_dir + '/**/*.tsv').first

    FileUtils.cp res, self.path
    nil
  end

  dep :HLA
  task :HLA_alleles => :array do
    step(:HLA).path.read.split("\n").last.split("\t").select{|e| e.include? "*"}.collect{|e| "HLA-" + e}
  end

end

#require 'NGS/tasks/basic.rb'

#require 'rbbt/knowledge_base/NGS'
#require 'rbbt/entity/NGS'

