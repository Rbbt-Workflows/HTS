require 'rbbt-util'
require 'rbbt/resource'
require 'tools/samtools'
require 'rbbt/sources/organism'

module BWA
  extend Resource
  self.subdir = 'organism/databases/BWA'

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

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".bwt") || Persist.newer?(linked + ".bwt", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, linked unless File.exists?(linked)
        FileUtils.ln_s file + '.alt', linked + '.alt' if File.exists?(file + '.alt') || ! File.exists?(linked + '.alt') 
        CMD.cmd("'#{BWA_CMD}' index '#{ linked }'")
      end
    end

    linked
  end

  Rbbt.claim Rbbt.software.opt.BWA, :install, Rbbt.share.install.software.BWA.find

  BWA_CMD = Rbbt.software.opt.BWA.produce.bwa.find

  Organism.claim Organism["Hsa"].hg19["hg19.fa"], :proc do |target|
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
      CMD.cmd("#{Samtools::samtools_cmd} faidx #{target}")
    end
    nil
  end

  Organism.claim Organism["Hsa"].b37["b37.fa"], :proc do |target|
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
      CMD.cmd("#{Samtools::samtools_cmd} faidx #{target}")
    end
    nil
  end
 
  Organism.claim Organism["Hsa"].hg38["hg38.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
    CMD.cmd("wget '#{url}' -O #{target}")
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
    CMD.cmd("wget '#{url}' -O #{target}.alt")
    nil
  end

  Organism.claim Organism["Hsa"].hs37d5["hs37d5.fa"], :proc do |target|
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
      CMD.cmd("#{Samtools::samtools_cmd} faidx #{target}")
    end
    nil
  end

  Organism.claim Organism["Mmu"].GRCm38["GRCm38.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
    CMD.cmd("wget '#{url}' -O #{target}.gz")
    nil
  end
end

if __FILE__ == $0
  Log.severity = 0
  Rbbt.software.opt.BWA.produce
  BWA.references.hs37d5["reference.fa"].produce
end

