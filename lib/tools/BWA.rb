require 'rbbt-util'
require 'rbbt/resource'
require 'tools/samtools'
require 'rbbt/sources/organism'

module BWA
  extend Resource
  self.subdir = 'organism/databases/BWA'

  Rbbt.claim Rbbt.software.opt.BWA, :install, Rbbt.share.install.software.BWA.find

  CMD.tool :bwa, Rbbt.software.opt.BWA, "bash -c 'type bwa'"

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib
  
  def self.mem_pipe(files, reference, args = "", io = nil, file = nil)
    CMD.cmd(:bwa,"mem #{args} '#{reference}' - > #{file} ", :in => io)
  end

  def self.mem(files, reference, args = "", file = nil)
    if file
      CMD.cmd(:bwa,"mem #{args} '#{reference}' #{files.collect{|f| "'#{f}'"} * " "} > #{file} ")
    else
      CMD.cmd(:bwa,"mem #{args} '#{reference}' #{files.collect{|f| "'#{f}'"} * " "} ", :pipe => true)
    end
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
    if ! File.exist?(linked + ".bwt") || Persist.newer?(linked + ".bwt", file)

      Misc.in_dir dir do
        if file != linked
          Open.rm linked
          Open.rm linked + '.alt'
          FileUtils.ln_s file, linked 
          FileUtils.ln_s file + '.alt', linked + '.alt' if File.exist?(file + '.alt') && ! File.exist?(linked + '.alt') && ! file == linked
        end
        CMD.cmd(:bwa, "index '#{ linked }'")
      end
    end

    linked
  end

  #Organism.claim Organism["Hsa"].hg19["hg19.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
  #  target.sub!(/\.gz$/,'')
  #  TmpFile.with_file do |directory|
  #    Misc.in_dir directory do
  #      CMD.cmd("wget '#{url}' -O - | tar xvfz -")
  #      CMD.cmd("cat *.gz | bgzip > '#{target}.gz' ")
  #    end
  #  end
  #  nil
  #end

  #Organism.claim Organism["Hsa"].b37["b37.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz"
  #  target.sub!(/\.gz$/,'')
  #  CMD.cmd("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
  #  nil
  #end
 
  #Organism.claim Organism["Hsa"].hg38_alt["hg38_alt.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
  #  target.sub!(/\.gz$/,'')
  #  CMD.cmd("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
  #  url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
  #  CMD.cmd("wget '#{url}' -O #{target}.alt")
  #  nil
  #end

  #Organism.claim Organism["Hsa"].hg38["hg38.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
  #  target.sub!(/\.gz$/,'')
  #  CMD.cmd("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
  #  nil
  #end

  #Organism.claim Organism["Hsa"].hs37d5["hs37d5.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz"
  #  target.sub!(/\.gz$/,'')
  #  CMD.cmd("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
  #  nil
  #end

  #Organism.claim Organism["Mmu"].GRCm38["GRCm38.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exist? File.dirname(target)
  #  url = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
  #  target.sub!(/\.gz$/,'')
  #  CMD.cmd("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
  #  nil
  #end
end

if __FILE__ == $0
  Log.severity = 0
  #Rbbt.software.opt.BWA.produce
  BWA.references.hs37d5["reference.fa"].produce
end

