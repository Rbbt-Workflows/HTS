module HTS

  dep :BAM_pileup_sumaries
  task :contamination => :text do
    args = {}
    args["input"] = step(:BAM_pileup_sumaries).path
    args["output"] = self.tmp_path
    GATK.run_log("CalculateContamination", args)
  end


  input :bam, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension 'pileup.gz'
  task :pileup => :text do |bam,reference|
    orig_reference = reference_file(reference)
    reference = Samtools.prepare_FASTA orig_reference

    monitor_cmd_genome "samtools mpileup -f '#{reference}' -Q 20 '#{bam}'", false, true
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 hs37d5), :nofile => true
  extension :vcf
  task :BAM_pileup_sumaries_known_biallelic => :tsv do |reference|
    variants_file = vcf_file reference, "1000g_snps"

    variants_file = GATK.prepare_VCF_AF_only variants_file

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference

    args = {}
    args["reference"] = reference
    args["variant"] = variants_file
    args["restrict-alleles-to"] = 'BIALLELIC'
    args["output"] = tmp_path
    GATK.run_log("SelectVariants", args)
    nil
  end

  dep :BAM_pileup_sumaries_known_biallelic, :jobname => "Default"
  input :BAM, :file, "BAM file", nil, :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_pileup_sumaries => :text do |bam,interval_list|

    variants_file = step(:BAM_pileup_sumaries_known_biallelic).path

    args = {}
    args["feature-file"] = variants_file
    GATK.run_log("IndexFeatureFile", args)

    args = {}
    args["input"] = Samtools.prepare_BAM bam 
    args["variant"] = variants_file
    args["output"] = self.tmp_path
    args["intervals"] = interval_list ? interval_list : variants_file
    args["interval-padding"] = 100 if interval_list
    GATK.run_log("GetPileupSummaries", args)
    nil
  end
end