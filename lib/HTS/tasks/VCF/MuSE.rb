module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :germline_resource, :file, "Germline resource", :dbsnp, :nofile => true
  extension :vcf
  task :muse => :text do |tumor,normal, reference, germline_resource|
    germline_resource = vcf_file reference, germline_resource
    germline_resource = GATK.prepare_VCF_AF_only germline_resource

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    reference = HTS.uncompress_FASTA orig_reference

    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    Open.mkdir files_dir
    int_fie = file('intermediate')
    muse = Rbbt.software.opt.MuSE.MuSE.find
    CMD.cmd_log("#{muse} call -O #{int_fie} -f #{reference} #{tumor} #{normal}")

    CMD.cmd_log("#{muse} sump -I #{int_fie}.MuSE.txt -E -O #{self.tmp_path}")
    nil
  end
end
