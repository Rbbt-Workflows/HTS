module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  extension :vcf
  task :mutect2 => :text do |tumor,normal,reference,interval_list|

    reference = reference_file reference

    reference = GATK.prepare_FASTA reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    #tumor_sample = CMD.cmd("#{Samtools::Samtools_CMD} view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #tumor_sample = CMD.cmd("samtools view -H '#{tumor}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    #normal_sample = CMD.cmd("samtools view -H '#{normal}' | grep '@RG'").read.match(/SM:([^\t]*)/)[1]
    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal) if normal

    args["input"] = [tumor, normal].compact
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample if normal_sample
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    FileUtils.mkdir_p files_dir unless File.exists? files_dir
    args["bam-output"] = file('haplotype.bam')

    GATK.run_log("Mutect2", args)
  end

  dep :mutect2
  dep :contamination, :BAM => :tumor
  extension :vcf
  task :mutect2_filtered => :tsv do
    args = {}
    FileUtils.mkdir_p files_dir

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2).path, tmp

    args["variant"] = tmp
    args["output"] = self.tmp_path
    #args["normal-artifact-lod"] = '10.0'
    #args["tumor-lod"] = '0.0'
    args["contamination-table"] = step(:contamination).path
    GATK.run_log("FilterMutectCalls", args)
  end

  dep :mutect2_filtered
  extension :vcf
  task :mutect2_clean => :tsv do
    TSV.traverse step(:mutect2_filtered), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end
end
