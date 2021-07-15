require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/workflow'
require 'rbbt/sources/organism'

Workflow.require_workflow "DbSNP"
module GATK
  extend Resource
  self.subdir = 'share/databases/GATK'

  Rbbt.claim Rbbt.software.opt.GATK, :install, Rbbt.share.install.software.GATK.find

  CMD.tool :gatk, Rbbt.software.opt.GATK, "gatk"


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

  def self.BAM_sample_name(bam_file, reference = nil)
    TmpFile.with_file do |tmp_file|
      GATK.run_log('GetSampleName', :input => bam_file, :output => tmp_file, :reference => reference)
      Open.read(tmp_file)
    end
  end

  def self.get_VCF(url, target)
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    begin
      CMD.cmd_log("wget '#{url}'  -O - | gunzip - -c | bgzip -c > '#{target}'")
    rescue
      FileUtils.rm target if File.exists?(target)
    end
    nil
  end

  def self.sort_VCF(file, dir=nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.vcf_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    Misc.in_dir dir do
      args ={}
      args["INPUT"] = file
      args["OUTPUT"] = linked
      GATK.run_log("SortVcf",args)
    end
    linked
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
        if linked != file
          Open.rm linked
          FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        end
        args = {}
        args["input"] = linked
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
                      if line =~ /^#CHROM/
                        "##INFO=<ID=NOINFO,Number=0,Type=Flag,Description=\"\">" + "\n" + line
                      else
                        line
                      end
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
        args["input"] = linked
        GATK.run_log("IndexFeatureFile", args)
      end
    end

    linked
  end

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

  
  Organism.claim Organism["Hsa"].hg38.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["small_exac_common_3.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["1000G_phase1.snps.high_confidence.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["af-only-gnomad.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["dbsnp_146.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for b37
  
  Organism.claim Organism["Hsa"].b37.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].b37.known_sites["small_exac_common_3.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].b37.known_sites["1000G_phase1.snps.high_confidence.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].b37.known_sites["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].b37.known_sites["af-only-gnomad.vcf.gz"], :proc do |target|
		url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].b37.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for hg19

  Organism.claim Organism["Hsa"].hg19.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg19.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz"
    GATK.get_VCF(url, target)
  end


#-- Claims for  hs37d5

  Organism.claim Organism["Hsa"].hs37d5.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hs37d5.known_sites["1000G_phase1.indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hs37d5.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for mm10

  Organism.claim Organism["Mmu"].GRCm38.known_sites["Ensembl.vcf.gz"], :proc do |target|
    url = "ftp://ftp.ensembl.org/pub/release-101/variation/vcf/mus_musculus/mus_musculus.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Mmu"].GRCm38.known_sites["Ensembl.structural.vcf.gz"], :proc do |target|
    url = "ftp://ftp.ensembl.org/pub/release-101/variation/vcf/mus_musculus/mus_musculus_structural_variations.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Rno"]["Rnor_6.0"].known_sites["Ensembl.vcf.gz"], :proc do |target|
    url = "ftp://ftp.ensembl.org/pub/release-101/variation/vcf/rattus_norvegicus/rattus_norvegicus.vcf.gz"
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
    dict = linked.replace_extension('dict', 'gz')
    if ! File.exists?(dict) || Persist.newer?(dict, file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        Open.rm dict if File.exists?(dict)
        CMD.cmd(:gatk, "CreateSequenceDictionary -R '#{ linked }'")
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
RevertSam
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

  def self.run(command, args = {}, sin = nil, tmp_dir = nil)

    progress_bar = args.delete(:progress_bar)

    arg_string = self.hash2args(args) if Hash === args

    log = false
    pipe = true

    tmp_dir ||= self.tmpdir
    if tmp_dir
      CMD.cmd(:gatk, "--java-options '-Dsamjdk.compression_level=1 -Djava.io.tmpdir=#{tmp_dir}' #{command} #{arg_string}", :log => log, :pipe => pipe, :in => sin, :progress_bar => progress_bar)
    else
      CMD.cmd(:gatk, "#{command} #{arg_string}", :log => log, :pipe => pipe, :in => sin, :progress_bar => progress_bar)
    end
  end

  def self.run_log(command, args = {}, sin = nil, tmp_dir = nil)

    progress_bar = args.delete(:progress_bar)

    arg_string = self.hash2args(args) if Hash === args

    tmp_dir ||= self.tmpdir
    if tmp_dir
      CMD.cmd_log(:gatk, "--java-options '-Dsamjdk.compression_level=1 -Djava.io.tmpdir=#{tmp_dir}' #{command} #{arg_string}", :in => sin, :progress_bar => progress_bar)
    else
      CMD.cmd_log(:gatk, "#{command} #{arg_string}", :in => sin, :progress_bar => progress_bar)
    end
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif GATK.known_sites.b37["1000G_phase1.snps.high_confidence.b37.vcf"].produce(true)
end
