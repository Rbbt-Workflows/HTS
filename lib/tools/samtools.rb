require 'rbbt-util'
require 'rbbt/resource'

module Samtools
  extend Resource
  self.subdir = 'share/databases/Samtools'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.HTSLib, :install, Rbbt.share.install.software.HTSLib.find
  Rbbt.claim Rbbt.software.opt.Samtools, :install, Rbbt.share.install.software.Samtools.find
  Rbbt.claim Rbbt.software.opt.bcftools, :install, Rbbt.share.install.software.bcftools.find

  CMD.tool :samtools, Rbbt.software.opt.Samtools do
    Rbbt.software.opt.HTSLib.produce
    Rbbt.software.opt.Samtools.produce
  end

  CMD.tool :bcftools, Rbbt.software.opt.Samtools

  #Samtools_CMD='samtools'
  
  #Samtools_CMD=Rbbt.software.opt.Samtools.produce.bin.samtools.find

  def self.samtools_cmd
    Rbbt.software.opt.Samtools.produce.bin.samtools.find
    'samtools'
  end

  def self.run(command)
    CMD.cmd_log(:samtools, command)
  end

  def self.prepare_BAM(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.bam_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".bai") || Persist.newer?(linked + ".bai", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, linked unless File.exists?(linked)
        cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
        if cpus
          Samtools.run("index -@ #{cpus} '#{ linked }'")
        else
          Samtools.run("index '#{ linked }'")
        end
      end
    end

    linked
  end
  
  def self.prepare_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".fai") || Persist.newer?(linked + ".fai", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        Samtools.run("faidx '#{ linked }'")
      end
    end

    linked
  end
  
  def self.BAM_sort(bam_file,to_sam=false, queryname=true)
  	cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    format_str=" -O BAM "
    if to_sam
      format_str=" -O SAM "
    end
    sort_order_str=""
    if queryname
      sort_order_str="-n"
    end
    if cpus
      Samtools.run("sort #{format_str} #{sort_order_str} -@ #{cpus} '#{bam_file}' -o #{bam_file}.sorted")
    else
      Samtools.run("sort #{format_str} #{sort_order_str} '#{bam_file}' -o #{bam_file}.sorted")
    end
  end

  def self.BAM_get_chr_reads(bam_file, chr)
    cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    CMD.cmd("samtools view -@ #{cpus} -Sb -h '#{bam_file}' #{chr} > #{chr}.bam")
  end

  def self.header(bam_file)
    CMD.cmd("samtools view -H #{bam_file}| grep -v \"^@PG\"").read
  end

  def self.reads_number(bam_file)
    CMD.cmd("samtools view -c #{bam_file}").read
  end

  def self.BAM_start(bam_file)
    CMD.cmd("samtools view '#{bam_file}'| head -n 1 | cut -f 3,4").read.strip.split("\t")
  end

  def self.merge(outbam, inbams)
    cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    CMD.cmd("samtools merge -n -@ #{cpus} #{outbam} #{inbams} ")
  end

  def self.viewSam(inBAM, outSAM)
    CMD.cmd("samtools view #{inBAM} > #{outSAM}")
  end

  def self.reference_contigs(reference)
    Open.read(reference + '.fai').split("\n").collect{|line| line.split("\t").first}
  end

  def self.bam_contigs(bam)
    bam = bam.path if Step === bam
    bam = bam.find if Path === bam
    CMD.cmd("samtools view -H '#{bam}' | grep ^@SQ | cut -f 2 |sed 's/SN://'").read.split("\n")
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HTSLib.produce
  iif Rbbt.software.opt.Samtools.produce 
end
