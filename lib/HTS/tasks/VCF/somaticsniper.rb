module HTS
  
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :quality, :integer, "Mapping quality filter threshold", 1
  extension :vcf
  task :somatic_sniper => :text do |tumor, normal, reference, quality|
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

    CMD.cmd_log("bam-somaticsniper -L -G -Q 15 -s 0.01 -T 0.85 -N 2 -r 0.001 -F vcf  \
                -q #{quality} \
                -n '#{normal_sample}' -t '#{tumor_sample}' -f '#{reference}' \
                '#{tumor}' '#{normal}' '#{self.tmp_path}'")

    nil
  end
end
