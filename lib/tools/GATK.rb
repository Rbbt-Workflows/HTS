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

  GATK.claim GATK.known_sites.b37["Miller_1000G_indels.vcf.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
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
  
  def self.hash2args(hash)
    hash.collect do |k,v| 
      k = '--' + k unless k[0] == "-"
      next if v.nil? || FalseClass === v
      v = nil if TrueClass === v
      vs = Array === v ? v : [v]
      vs.collect do |v|
        v = "'" + v + "'" unless v[0] == "'" || v[0] == '"' unless v.nil?
        [k,v].compact * " " 
      end
    end.compact.flatten * " "
  end

  def self.run(command, arg_string = "", sin = nil)
    arg_string = hash2args(arg_string) if Hash === arg_string
    CMD.cmd("#{GATK_CMD} #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
  end

  def self.run_log(*args)
    io = run(*args)
    while line = io.gets
      Log.debug line
    end
    io.join
    nil
  end

  GATK_CMD=Rbbt.software.opt.GATK.produce.gatk.find
end

if __FILE__ == $0
  Log.severity = 0
  iif GATK.tutorials.tutorial_11682.produce
  iif GATK.tutorials.tutorial_11136.produce
end
