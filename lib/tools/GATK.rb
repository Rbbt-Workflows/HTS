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

  def self.tmpdir
    begin
      tmpdir = Rbbt::Config.get('tmpdir', :gatk)
      if tmpdir && ! File.exists?(tmpdir)
        Open.mkdir tmpdir
        File.chmod(0777, tmpdir)
      end
    end
    tmpdir
  end

  def self.BAM_sample_name(bam_file)
    TmpFile.with_file do |tmp_file|
      GATK.run_log('GetSampleName', :input => bam_file, :output => tmp_file)
      Open.read(tmp_file)
    end
  end

  def self.get_VCF(url, target)
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    CMD.cmd("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    nil
  end

  def self.prepare_VCF(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.vcf_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".tbi") || Persist.newer?(linked + '.tbi', file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        args = {}
        args["feature-file"] = linked
        GATK.run_log("IndexFeatureFile", args)
      end
    end

    linked
  end

  def self.prepare_VCF_AF_only(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.vcf_indices_af_only[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".tbi") || Persist.newer?(linked + '.tbi', file)

      Misc.in_dir dir do

        if !File.exists?(linked)
          Open.write(linked + '.tmp') do |fout|
            TSV.traverse file, :type => :array do |line|
              out = if line[0] == "#"
                      line
                    else
                      parts = line.split("\t", -1)
                      info = parts[7].split(";").select{|f| %w(AC AF CAF).include? f.split("=").first } * ";"
                      info = "NOINFO=true" if info.empty?
                      parts[7] = info
                      parts * "\t"
                    end
              fout.puts out
            end
          end

          if Open.gzip?(linked)
            CMD.cmd("cat #{linked + '.tmp'} | bgzip > #{linked}")
            FileUtils.rm linked + '.tmp'
          else
            FileUtils.mv linked + '.tmp', linked
          end
        end

        args = {}
        args["feature-file"] = linked
        GATK.run_log("IndexFeatureFile", args)
      end
    end

    linked
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
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hg38["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hg38["1000G_phase1.snps.high_confidence.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

	GATK.claim GATK.known_sites.hg38["af-only-gnomad.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hg38["dbsnp_146.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for b37
  
  GATK.claim GATK.known_sites.b37["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.b37["1000G_phase1.snps.high_confidence.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.b37["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

	GATK.claim GATK.known_sites.b37["af-only-gnomad.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.b37["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for hg19

  GATK.claim GATK.known_sites.hg19["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hg19["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz"
    GATK.get_VCF(url, target)
  end


#-- Claims for  hs37d5

  GATK.claim GATK.known_sites.hs37d5["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hs37d5["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  GATK.claim GATK.known_sites.hs37d5["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
    GATK.get_VCF(url, target)
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
    if ! File.exists?(linked.replace_extension("dict", true)) || Persist.newer?(linked.replace_extension('dict', true), file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd("'#{GATK_CMD}' CreateSequenceDictionary -R '#{ linked }'")
      end
    end

    linked
  end

  SPARK_COMMANDS = %w(
CollectAllelicCounts
CountBases
CountReads
Pileup
CalcMetadata
CollectBaseDistributionByCycle
CollectInsertSizeMetrics
CollectMultipleMetrics
CollectQualityYieldMetrics
CompareDuplicates
FlagStat
MeanQualityByCycle
QualityScoreDistribution
PathSeqBwa
PathSeqFilter
PathSeqPipeline
PathSeqScore
ParallelCopyGCSDirectoryIntoHDFS
ApplyBQSR
BQSRPipeline
BaseRecalibrator
BwaAndMarkDuplicatesPipeline
Bwa
ExtractOriginalAlignmentRecordsByName
MarkDuplicates
PrintReads
__RevertSam
SortSam
FindBadGenomicKmers
__HaplotypeCaller
ReadsPipeline
CpxVariantReInterpreter
DiscoverVariantsFromContigAlignmentsSAM
ExtractSVEvidence
FindBreakpointEvidence
StructuralVariationDiscoveryPipeline
SvDiscoverFromLocalAssemblyContigAlignments
CountVariants
PrintVariants
)

  def self.hash2args(hash)
    hash.collect do |k,v| 
      if k[0] != "-"
        if k.length == 1
          k = '-' + k.to_s 
        else
          k = '--' + k.to_s 
        end
      end
      next if v.nil? 
      #v = nil if TrueClass === v
      vs = Array === v ? v : [v]
      vs.collect do |v|
        v = "'" + v.to_s + "'" unless v.to_s[0] == "'" || v.to_s[0] == '"' unless v.nil?
        [k,v].compact * " " 
      end
    end.compact.flatten * " "
  end

  def self.run(command, args = {}, sin = nil)

    arg_string = self.hash2args(args) if Hash === args

    tmpdir = self.tmpdir
    if tmpdir
      CMD.cmd("#{GATK_CMD} --java-options '-Djava.io.tmpdir=#{tmpdir}' #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    else
      CMD.cmd("#{GATK_CMD} #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    end
  end

  def self.run_log(command, args = {}, sin = nil)

    arg_string = self.hash2args(args) if Hash === args

    tmpdir = self.tmpdir
    if tmpdir
      CMD.cmd_log("#{GATK_CMD} --java-options '-Djava.io.tmpdir=#{tmpdir}' #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    else
      CMD.cmd_log("#{GATK_CMD} #{command} #{arg_string}", :log => true, :pipe => true, :in => sin)
    end
  end

  GATK_CMD='gatk'
end

if __FILE__ == $0
  Log.severity = 0
  iif GATK.known_sites.b37["1000G_phase1.snps.high_confidence.b37.vcf"].produce(true)
end
