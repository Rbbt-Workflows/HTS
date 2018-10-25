require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/workflow'

Workflow.require_workflow "DbSNP"
module GATK
  extend Resource
  self.subdir = 'share/databases/GATK'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  def self.BAM_sample_name(bam_file)
    TmpFile.with_file do |tmp_file|
      GATK.run_log('GetSampleName', :input => bam_file, :output => tmp_file)
      Open.read(tmp_file)
    end
  end

  Rbbt.claim Rbbt.software.opt.GATK, :install, Rbbt.share.install.software.GATK.find

  DbSNP.claim DbSNP["All.vcf.gz.tbi"], :proc do
    args = {}
    args["feature-file"] = DbSNP["All.vcf.gz"].produce.find
    GATK.run_log("IndexFeatureFile", args)
  end

  %w(11682 11683 11136).each do |tutorial_number|
    GATK.claim GATK.tutorials["tutorial_#{tutorial_number}"], :proc do |target|
      url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/tutorials/datasets/tutorial_#{tutorial_number}.tar.gz"
      Misc.in_dir target do
        Misc.untar Open.open(url), target
        CMD.cmd("mv #{File.basename(target)}/* .; rmdir #{File.basename(target)}")
      end
      nil
    end
  end

  # -- Claims for HG38
  GATK.claim GATK.bundle.hg38, :proc do |target|
    Misc.in_dir target do
      CMD.cmd_log("wget -m -nH --cut-dirs=2 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/")
    end
  end

  
  GATK.claim GATK.known_sites.hg38["Miller_1000G_indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.hg38["1000G_phase1.indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.hg38["dbsnp_138.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  # -- Claims for b37
  
  GATK.claim GATK.known_sites.b37["Miller_1000G_indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.b37["1000G_phase1.snps.high_confidence.b37.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.b37["1000G_phase1.indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.b37["dbsnp_138.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  # -- Claims for hg19

  GATK.claim GATK.known_sites.hg19["Miller_1000G_indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.hg19["1000G_phase1.indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  GATK.claim GATK.known_sites.hg19["dbsnp_138.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz"
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    args = {}
    args["feature-file"] = target
    GATK.run_log("IndexFeatureFile", args)
    nil
  end

  def self.hash2args(hash)
    hash.collect do |k,v| 
      k = '--' + k.to_s unless k[0] == "-"
      next if v.nil? || FalseClass === v
      v = nil if TrueClass === v
      vs = Array === v ? v : [v]
      vs.collect do |v|
        v = "'" + v.to_s + "'" unless v[0] == "'" || v[0] == '"' unless v.nil?
        [k,v].compact * " " 
      end
    end.compact.flatten * " "
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
    if ! File.exists?(linked.replace_extension("dict")) || Persist.newer?(linked.replace_extension('dict'), file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd("'#{GATK_CMD}' CreateSequenceDictionary -R '#{ linked }'")
      end
    end

    linked
  end

  def self.run(command, arg_string = "", sin = nil)
    arg_string = hash2args(arg_string) if Hash === arg_string
    tmpdir = Rbbt::Config.get('tmpdir', :gatk)
    if tmpdir
      CMD.cmd("#{GATK_CMD} --java-options '-Djava.io.tmpdir=#{tmpdir}' #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    else
      CMD.cmd("#{GATK_CMD} #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    end
  end

  def self.run_log(*args)
    io = run(*args)
    while line = io.gets
      Log.debug line
    end
    io.join
    nil
  end

  def self.run_log(command, arg_string = "", sin = nil)
    arg_string = hash2args(arg_string) if Hash === arg_string
    tmpdir = Rbbt::Config.get('tmpdir', :gatk)
    if tmpdir
      CMD.cmd_log("#{GATK_CMD} --java-options '-Djava.io.tmpdir=#{tmpdir}' #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    else
      CMD.cmd_log("#{GATK_CMD} #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    end
  end

  GATK_CMD=Rbbt.software.opt.GATK.produce.gatk.find
end

if __FILE__ == $0
  Log.severity = 0
  iif GATK.known_sites.b37["1000G_phase1.snps.high_confidence.b37.vcf"].produce(true)
end
