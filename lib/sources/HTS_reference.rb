require 'rbbt/sources/organism'
require_relative '../tools/GATK'

module HTS
  def self.uncompress_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.digest(Open.realpath(file))
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    
    if linked =~ /.gz$/
      unzipped_linked = linked.sub(/.gz$/,'')
      CMD.cmd("zcat #{linked} > #{unzipped_linked}") unless File.exists?(unzipped_linked)
      dir.glob("*.gz.*").each do |file|
        unzipped_file = file.sub('.gz.', '.')
        Open.ln_s file, unzipped_file unless File.exists?(unzipped_file)
      end
      unzipped_linked
    else
      linked
    end
  end

end

module Organism

  # -- Claims for hg38
  
  Organism.claim Organism["Hsa"].hg38["hg38.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
    CMD.cmd("wget '#{url}' -O #{target}")
    nil
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
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
  
  Organism.claim Organism["Hsa"].b37["b37.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz"
    CMD.cmd("wget '#{url}' -O  - | gunzip -c > #{target}")
    nil
  end
 
  Organism.claim Organism["Hsa"].b37.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
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

  Organism.claim Organism["Hsa"].hg19["hg19.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
    TmpFile.with_file do |directory|
      Misc.in_dir directory do
        CMD.cmd("wget '#{url}' -O - | tar xvfz -")
        CMD.cmd("cat *.fz > '#{target}' ")
      end
    end
    nil
  end

  Organism.claim Organism["Hsa"].hg19.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    GATK.get_VCF(url, target)
  end

  Organism.claim Organism["Hsa"].hg19.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz"
    GATK.get_VCF(url, target)
  end

  # -- Claims for  hs37d5

  Organism.claim Organism["Hsa"].hs37d5["hs37d5.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz"
    CMD.cmd("wget '#{url}' -O  - | gunzip -c > #{target}")
    nil
  end

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

end