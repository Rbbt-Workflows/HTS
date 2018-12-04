require 'rbbt-util'
require 'rbbt/resource'
require 'tools/samtools'

module BWA
  extend Resource
  self.subdir = 'share/databases/BWA'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib
  
  def self.mem(files, reference, args = "")
    CMD.cmd("#{BWA_CMD} mem #{args} '#{reference}' #{files.collect{|f| "'#{f}'"} * " "} ", :pipe => true)
  end

  def self.prepare_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.digest(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".bwt") || Persist.newer?(linked + ".bwt", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd("'#{BWA_CMD}' index '#{ linked }'")
      end
    end

    linked
  end

  Rbbt.claim Rbbt.software.opt.BWA, :install, Rbbt.share.install.software.BWA.find

  BWA_CMD = Rbbt.software.opt.BWA.produce.bwa.find

  BWA.claim BWA.references.hg19["hg19.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
    TmpFile.with_file do |directory|
      Misc.in_dir directory do
        CMD.cmd("wget '#{url}' -O - | tar xvfz -")
        CMD.cmd("cat *.fz > '#{target}' ")
      end
    end
    Misc.in_dir File.dirname(target) do
      CMD.cmd("#{BWA_CMD} index -p reference -a bwtsw #{target}")
      io = GATK.run("CreateSequenceDictionary", {"R" => target})
      while line = io.gets
        Log.debug line
      end
      CMD.cmd("#{Samtools::Samtools_CMD} faidx #{target}")
    end
    nil
  end

  BWA.claim BWA.references.b37["b37.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz"
    CMD.cmd("wget '#{url}' -O  - | gunzip -c > #{target}")
    Misc.in_dir File.dirname(target) do
      CMD.cmd("#{BWA_CMD} index -p reference -a bwtsw #{target}")
      io = GATK.run("CreateSequenceDictionary", {"R" => target})
      while line = io.gets
        Log.debug line
      end
      CMD.cmd("#{Samtools::Samtools_CMD} faidx #{target}")
    end
    nil
  end
 
  BWA.claim BWA.references.hg38["hg38.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
    CMD.cmd("wget '#{url}' -O #{target}")
    nil
  end

  BWA.claim BWA.references.hs37d5["hs37d5.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz"
    CMD.cmd("wget '#{url}' -O  - | gunzip -c > #{target}")
    Misc.in_dir File.dirname(target) do
      CMD.cmd("#{BWA_CMD} index -p reference -a bwtsw #{target}")
      io = GATK.run("CreateSequenceDictionary", {"R" => target})
      while line = io.gets
        Log.debug line
      end
      CMD.cmd("#{Samtools::Samtools_CMD} faidx #{target}")
    end
    nil
  end
end

if __FILE__ == $0
  Log.severity = 0
  Rbbt.software.opt.BWA.produce
  BWA.references.hs37d5["reference.fa"].produce
end

