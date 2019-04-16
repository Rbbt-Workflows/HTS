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
  Rbbt.claim Rbbt.software.opt.bcftools, :install, Rbbt.share.install.software.bcftools.find
  Rbbt.claim Rbbt.software.opt.Samtools, :install, Rbbt.share.install.software.Samtools.find

  #Samtools_CMD='samtools'
  
  #Samtools_CMD=Rbbt.software.opt.Samtools.produce.bin.samtools.find

  def self.samtools_cmd
    Rbbt.software.opt.Samtools.produce.bin.samtools.find
    'samtools'
  end

  def self.run(command)
    CMD.cmd_log("'#{samtools_cmd}' " << command)
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
        Samtools.run("index '#{ linked }'")
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

  def self.BAM_start(bam_file)
    CMD.cmd("samtools view '#{bam_file}'| head -n 1 | cut -f 3,4").read.strip.split("\t")
  end

  def self.reference_contigs(reference)
    Open.read(reference + '.fai').split("\n").collect{|line| line.split("\t").first}
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HTSLib.produce
  iif Rbbt.software.opt.Samtools.produce 
end
