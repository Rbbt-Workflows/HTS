module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :pon, :file, "Panel of normals", nil, :nofile => true
  input :germline_resource, :file, "Germline resource", nil, :nofile => true
  input :af_not_in_resource, :float, "Allele frequency of alleles not in resource", nil
  extension :vcf
  task :mutect2 => :text do |tumor,normal,reference,interval_list,pon,germline_resource,af_not_in_resource|

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal) if normal

    args["input"] = [tumor, normal].compact
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample if normal_sample
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    args["panel-of-normals"] = pon if pon
    FileUtils.mkdir_p files_dir unless File.exists? files_dir
    args["bam-output"] = file('haplotype.bam')
    args["germline-resource"] = germline_resource
    args["af-of-alleles-not-in-resource"] = af_not_in_resource.to_s if af_not_in_resource
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
