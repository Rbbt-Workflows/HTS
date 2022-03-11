module HTS
  CMD.tool "configManta.py", nil, "configManta.py -h" do
    CMD.conda(:manta, :python2)
  end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :somatic_config => :text do |tumor,normal,reference|
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal
    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal
    
    reference = reference_file reference
    orig_reference = reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference
    args = "--tumorBam #{tumor} --referenceFasta #{reference} --runDir #{files_dir} " 
    args += "--normalBam #{normal}" if normal
    Open.mkdir files_dir
    if !File.exists? (files_dir + "/runWorkflow.py") 
      CMD.cmd_log("configManta.py", args)
    end
    files_dir + "/runWorkflow.py"
  end

  input :bams, :array, "BAM files", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :germline_config => :text do |bams,reference|

    bams_str = ""
    bams.each {|bam| 
      bam = bam.path if Step === bam
      bam = Samtools.prepare_BAM(bam)
      bams_str += "--bam #{bam} "
    }
    
    reference = reference_file reference
    orig_reference = reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference
    args += "--referenceFasta #{reference} --runDir #{files_dir} " 
    args += bams_str

    Open.mkdir files_dir
    if !File.exists? (files_dir + "/runWorkflow.py")
      CMD.cmd_log("configManta.py", args)
    end
    files_dir + "/runWorkflow.py"
  end

  dep :somatic_config do |jobname, options|
    if options[:type].to_s == "somatic"
      {:inputs => options, :jobname => jobname}
    else
      nil
    end
  end
  dep :germline_config do |jobname, options|
    if options[:type].to_s == "germline"
      {:inputs => options, :jobname => jobname}
    else
      nil
    end
  end
  input :type, :select, "Type of analysis", :somatic, :select_options => %w(germline somatic)
  extension "vcf.gz"
  task :manta_pre => :text do 
    CMD.cmd_log(dependencies.flatten.first.files_dir + "/runWorkflow.py")
    FileUtils.ln_s dependencies.flatten.first.files_dir + "/results/variants/somaticSV.vcf.gz", self.path
  end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension "vcf.gz"
  task :manta_somatic => :text do |tumor,normal,reference|
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal
    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    reference = reference_file reference
    orig_reference = reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    args = "--tumorBam #{tumor} --referenceFasta #{reference} --runDir #{files_dir} " 
    args += "--normalBam #{normal}" if normal

    Open.mkdir files_dir
    if !File.exists? file("runWorkflow.py") 
      CMD.cmd_log("configManta.py", args)
    end
    cmd = file("runWorkflow.py")

    cpus = config :cpus, :manta, :Manta, :default => 1

    CMD.cmd_log(cmd, "-j #{cpus}")
    FileUtils.ln_s file("results/variants/somaticSV.vcf.gz"), self.tmp_path
    nil
  end
end
