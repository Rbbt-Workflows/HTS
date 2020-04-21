require 'tools/SomaticSeq'
module HTS

  input :tumor_bam_file, :file, "Tumor BAM", nil, :nofile => true
  input :normal_bam_file, :file, "Normal BAM (optional)", nil, :nofile => true
  input :genome_reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :mutect2_vcf, :file, "Mutect2 vcf file",nil, :nofile => true
  input :varscan_snv, :file, "Varscan vcf file containing SNVs",nil, :nofile => true
  input :varscan_indel, :file, "Varscan vcf file containing SNVs",nil, :nofile => true
  input :somaticsniper_vcf, :file, "Somatic Sniper vcf file", nil, :nofile => true
  input :muse_vcf, :file, "Muse vcf file", nil, :nofile => true
  input :strelka_snv, :file, "Strelka vcf file", nil, :nofile => true
  input :strelka_indel, :file, "Strelka vcf file", nil, :nofile => true
  extension :vcf
  task :somatic_seq => :text do |tumor_bam_file,normal_bam_file,genome_reference,mutect2_vcf,varscan_snv,varscan_indel,somaticsniper_vcf,muse_vcf,strelka_snv,strelka_indel|
    reference = reference_file genome_reference

    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    args = {}

    args["tumor-bam-file"] = Samtools.prepare_BAM(tumor_bam_file)
    args["normal-bam-file"] = Samtools.prepare_BAM(normal_bam_file)
    args["genome-reference"] = reference
    args["mutect2-vcf"] = mutect2_vcf
    args["varscan-snv"] = varscan_snv
    args["varscan-indel"] = varscan_indel
    args["somaticsniper-vcf"] = somaticsniper_vcf
    args["muse-vcf"] = muse_vcf
    args["strelka-snv"] = strelka_snv
    args["strelka-indel"] = strelka_indel

    outdir = self.path.to_s[0..self.path.to_s.rindex(".")-1]
    Open.mkdir outdir
    args["outdir"] = outdir
    SomaticSeq.run(args)
    FileUtils.ln_s outdir + "/Consensus.sSNV.vcf", self.path
  end
end
