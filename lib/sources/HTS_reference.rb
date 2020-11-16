require 'rbbt/sources/organism'
require_relative '../tools/GATK'

module HTS
  def self.uncompress_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
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

  def self.unfold_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find

    dir = linked + ".by_contig"
    if ! File.exists?(dir) || Persist.newer?(dir, file)
      Open.mkdir dir
      Misc.in_dir dir do
        CMD.cmd_log(CMD.bash("awk '/^>/ {OUT=substr(\\$1,2) \".fa\"}; {print >> OUT; close(OUT)}' #{CMD.gzip_pipe(linked)}"))
      end
    end

    dir
  end

  def self.gtf_file(organism)
    reference, organism = [organism, nil] if %w(hg19 hg38 b37 mm10).include?(organism)

    reference = case Organism.hg_build(organism)
                when 'hg19'
                  'b37'
                else
                  Organism.hg_build(organism)
                end if reference.nil?

    organism = Organism.organism_for_build(reference) if organism.nil?
    file = Organism.gene_set(organism).produce.find

    return file unless reference == "hg38"

    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.gtf_files[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir
    Open.mkdir dir

    linked = dir[basename].find
    Open.ln_s file, linked

    if ! File.exists?(linked + ".fixed.gtf") || Persist.newer?(linked + '.fixed.gtf', file)
      Open.open(linked) do |sout|
        Open.open(linked + '.fixed.gtf', :mode => 'w') do |sin|
          TSV.traverse sout, :type => :array do |line|
            if line =~ /^[0-9A-Z]/
              sin.puts 'chr' << line
            else 
              sin.puts line
            end
          end
        end
      end
    end

    linked + '.fixed.gtf'
  end

  def self.gsutil(url, target)
    CMD.cmd_log(:gsutil, "cp '#{url}' '#{target}'")
  end

end

module Organism

  # -- Claims for hg38
  
  Organism.claim Organism["Hsa"].hg38_noalt["hg38_noalt.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz"
    TmpFile.with_file do |tmpdir|
      Open.mkdir tmpdir
      Misc.in_dir(tmpdir) do
        CMD.cmd_log("wget '#{url}' -O file.tar.gz; tar xvfz file.tar.gz")
        CMD.cmd("cat */*.fa | bgzip > #{target}.gz")
      end
    end
    nil
  end

  Organism.claim Organism["Hsa"].hg38["hg38.fa.gz"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
    target.sub!(/\.gz$/,'')
    CMD.cmd_log("wget '#{url}' -O - | bgzip > #{target}.gz")
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt"
    CMD.cmd_log("wget '#{url}' -O #{target}.alt")
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

  Organism.claim Organism["Hsa"].hg38.known_sites["panel_of_normals.vcf"], :proc do |target|
    url = "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
    HTS.gsutil(url, target)
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["exome_capture.interval_list"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Broad.human.exome.hg38.interval_list.gz"
    target = target + '.gz' unless target =~ /\.gz$/
    CMD.cmd_log("wget '#{url}' -O  #{target}")
    nil
  end

  Organism.claim Organism["Hsa"].hg38.known_sites["wgs_calling_regions.interval_list"], :proc do |target|
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list"
    CMD.cmd_log("wget '#{url}' -O  #{target}")
    nil
  end




  # -- Claims for b37
  
  Organism.claim Organism["Hsa"].b37["b37.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    target.sub!(/\.gz$/,'')
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz"
    CMD.cmd_log("wget '#{url}' -O  - | gunzip -c > #{target}")
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

  Organism.claim Organism["Hsa"].b37.known_sites["panel_of_normals.vcf"], :proc do |target|
    url = "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf"
    HTS.gsutil(url, target)
    nil
  end

  Organism.claim Organism["Hsa"].b37.known_sites["panel_of_normals.vcf"], :proc do |target|
    url = "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf"
    HTS.gsutil(url, target)
    nil
  end

  Organism.claim Organism["Hsa"].b37.known_sites["exome_capture.interval_list"], :proc do |target|
    url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Broad.human.exome.b37.interval_list.gz"
    target = target + '.gz' unless target =~ /\.gz$/
    CMD.cmd_log("wget '#{url}' -O  #{target}")
    nil
  end

  Organism.claim Organism["Hsa"].b37.known_sites["wgs_calling_regions.interval_list"], :proc do |target|
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/b37/v0/wgs_calling_regions.hg38.interval_list"
    CMD.cmd_log("wget '#{url}' -O  #{target}")
    nil
  end





  # -- Claims for hg19

  #Organism.claim Organism["Hsa"].hg19["hg19.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
  #  url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
  #  TmpFile.with_file do |directory|
  #    Misc.in_dir directory do
  #      CMD.cmd_log("wget '#{url}' -O - | tar xvfz -")
  #      CMD.cmd("cat *.fz > '#{target}' ")
  #    end
  #  end
  #  nil
  #end

  #Organism.claim Organism["Hsa"].hg19.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
  #  GATK.get_VCF(url, target)
  #end

  #Organism.claim Organism["Hsa"].hg19.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz"
  #  GATK.get_VCF(url, target)
  #end

  ## -- Claims for  hs37d5

  #Organism.claim Organism["Hsa"].hs37d5["hs37d5.fa"], :proc do |target|
  #  FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
  #  target.sub!(/\.gz$/,'')
  #  url = "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz"
  #  CMD.cmd_log("wget '#{url}' -O  - | gunzip -c > #{target}")
  #  nil
  #end

  #Organism.claim Organism["Hsa"].hs37d5.known_sites["Miller_1000G_indels.vcf.gz"], :proc do |target|
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
  #  GATK.get_VCF(url, target)
  #end

  #Organism.claim Organism["Hsa"].hs37d5.known_sites["1000G_phase1.indels.vcf.gz"], :proc do |target|
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
  #  GATK.get_VCF(url, target)
  #end

  #Organism.claim Organism["Hsa"].hs37d5.known_sites["dbsnp_138.vcf.gz"], :proc do |target|
  #  url = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz"
  #  GATK.get_VCF(url, target)
  #end

  Organism.claim Organism["Mmu"].GRCm38["GRCm38.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
    target.sub!(/\.gz$/,'')
    CMD.cmd_log("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
    nil
  end

  Organism.claim Organism["Rno"]["Rnor_6.0"]["Rnor_6.0.fa"], :proc do |target|
    FileUtils.mkdir_p File.dirname(target) unless File.exists? File.dirname(target)
    url = "ftp://ftp.ensembl.org/pub/release-101/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fa.gz"
    target.sub!(/\.gz$/,'')
    CMD.cmd_log("wget '#{url}' -O - | gunzip -c | bgzip > #{target}.gz")
    nil
  end
end
