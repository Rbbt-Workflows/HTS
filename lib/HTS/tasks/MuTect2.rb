module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38), :nofile => true
  extension :vcf
  task :mutect2 => :text do |tumor,normal,reference|

    reference = reference_file reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal)

    #tumor_sample = CMD.cmd("#{Samtools::Samtools_CMD} view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #tumor_sample = CMD.cmd("samtools view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #normal_sample = CMD.cmd("samtools view -H '#{normal}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal)

    args["input"] = [tumor, normal]
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample
    FileUtils.mkdir_p files_dir unless File.exists? files_dir
    args["bam-output"] = file('haplotype.bam')

    GATK.run_log("Mutect2", args)
  end

  dep :mutect2
  extension :vcf
  task :mutect2_filtered => :tsv do
    args = {}
    FileUtils.mkdir_p files_dir

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2).path, tmp

    args["variant"] = tmp
    args["output"] = self.tmp_path
    GATK.run_log("FilterMutectCalls", args)
  end

  dep :mutect2_filtered
  extension :vcf
  task :mutect2_clean => :tsv do
    TSV.traverse step(:mutect2_filtered), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end
end
