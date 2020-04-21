require 'tools/Strelka'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :interval_list, :file, "interval list bed file", nil, :nofile => true
  extension :vcf
  task :strelka_pre => :text do |tumor,normal,reference,interval_list|
    output = file('output')
    reference = reference_file reference
    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    if interval_list and interval_list.include? '.bed'
      interval_list = HTS.prepare_BED interval_list if interval_list
    else
      message "Interval list not in bed format" if interval_list
      interval_list = nil
    end

    reference = Samtools.prepare_FASTA(reference)

    cpus = config :cpus, :strelka, :default => 3
    Strelka.runSomatic(tumor, normal, reference, output, cpus, interval_list)
    
    Open.read(output.results.variants["somatic.snvs.vcf.gz"])
  end

  dep :strelka_pre
  extension :vcf
  task :strelka_pre_indels => :text do 
    Open.read(step(:strelka_pre).file('output').results.variants["somatic.indels.vcf.gz"])
  end

  dep :strelka_pre
  extension :vcf
  input :strelka_filter_tier, :boolean, "Take only tier 1 variants from strelka", true
  input :strelka_filter_evs, :integer, "Strelka EVS minimum score filter", 15
  input :strelka_filter_qs, :integer, "Strelka QS(S|I) minimum score filter", 15
  task :strelka_filtered => :text do |tier,evs,qs|
    TSV.traverse step(:strelka_pre), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/

      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/
      qs_nt = line.split(";").select{|d| d =~ /^QS(S|I)_NT/}.first.split("=").last.to_i
      next unless qs_nt >= qs
      tqs_nt = line.split(";").select{|d| d =~ /^TQS(S|I)_NT/}.first.split("=").last.to_i
      next if tier and tqs_nt > 1 
      tqs = line.split(";").select{|d| d =~ /^TQS(S|I)/}.first.split("=").last.to_i
      next if tier and tqs > 1 
      somatic_evs = line.split(";").select{|d| d =~ /^SomaticEVS/}.first.split.first.split("=").last.to_f
      next unless somatic_evs >= evs

      line
    end
  end

  dep :strelka_pre_indels
  input :strelka_filter_evs, :integer, "Strelka EVS minimum score filter", 0
  input :strelka_filter_qs, :integer, "Strelka QS(S|I) minimum score filter", 0
  dep_task :strelka_filtered_indels, HTS, :strelka_filtered do |jobname,options,dependencies|
    strelka_indels = dependencies.select{|dep| dep.task_name == :strelka_pre_indels}.first
    options["HTS#strelka_pre"] = strelka_indels
    {:inputs => options, :jobname => jobname}
  end

  dep :strelka_filtered, :compute => :produce
  dep :strelka_filtered_indels do |jobname, options|
    job = HTS.job(:strelka_filtered_indels, jobname, options)
    job.overriden = false
    job
  end
  extension :vcf
  task :strelka => :text do
    
    args = {}
    args["INPUT"] = dependencies.collect{|dep| dep.path}
    args["OUTPUT"] = self.tmp_path
    gatk("MergeVcfs", args)
    nil
  end
end
