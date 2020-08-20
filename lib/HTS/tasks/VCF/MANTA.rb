module HTS
  CMD.tool :Manta, "configManta.py", "configManta.py -h" do
    CMD.cmd('conda install manta -c bioconda')
  end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :vcf
  task :somatic_config do |tumor,normal,reference|

    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal
    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal
    
    reference = reference_file reference
    orig_reference = reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference
    args = "--tumorBam #{tumor} --referenceFasta #{reference}" 
    args += "--normalBam #{normal}" if normal

    CMD.cmd_log(:Manta, args)
  end

  input :bams, :array, "BAM files", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :vcf
  task :germline_config do |bams,reference|

    bams_str = ""
    bams.each |bam| do
      bam = bam.path if Step === bam
      bam = Samtools.prepare_BAM(bam)
      bams_str += "--bam #{bam} "
    end
    
    reference = reference_file reference
    orig_reference = reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference
    args += "--referenceFasta #{reference} " 
    args += bams_str

    CMD.cmd_log(:Manta, args)
  end
end
