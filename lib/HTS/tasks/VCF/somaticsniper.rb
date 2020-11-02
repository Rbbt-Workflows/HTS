require 'tools/fpfilter'
module HTS
  
  CMD.tool "bam-somaticsniper".to_sym, nil, "bam-somaticsniper -v" do
    CMD.cmd('conda install somatic-sniper -c bioconda')
  end

  desc <<-EOF
Somatic sniper variant caller

somatic_score is set to 40 according to documentation in the web page
  EOF
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :quality, :integer, "Mapping quality filter threshold", 1
  input :somatic_score, :integer, "Filtering somatic snv output with somatic quality less than", 40 
  extension :vcf
  task :somatic_sniper => :text do |tumor, normal, reference, quality, somatic_score|
    reference = reference_file reference

    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal) if normal


    CMD.cmd_log("bam-somaticsniper".to_sym,"-F vcf -L -G \
                -s 0.01 -T 0.85 -N 2 -r 0.001 \
                -q #{quality} -Q #{somatic_score}  \
                -n '#{normal_sample}' -t '#{tumor_sample}' \
                -f '#{reference}' '#{tumor}' '#{normal}' '#{self.tmp_path}'")

    nil
  end

  dep :somatic_sniper
  extension :vcf
  task :somatic_sniper_filtered => :text do
    orig_reference = reference_file step(:somatic_sniper).inputs[:reference]
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    args = {}
    args["bam-file"] = step(:somatic_sniper).inputs[:tumor]
    args["vcf-file"] = step(:somatic_sniper).path
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["sample"] = GATK.BAM_sample_name(step(:somatic_sniper).inputs[:tumor])
    FPFilter.filter(args.to_hash)

    args = {}
    args["reference"] = reference
    args["variant"] = self.tmp_path
    args["output"] = self.path
    args["exclude-filtered"] = true
    args["exclude-non-variants"] = true
    GATK.run_log("SelectVariants", args)
  end
end
