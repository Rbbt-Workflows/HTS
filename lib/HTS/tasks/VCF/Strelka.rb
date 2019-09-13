require 'tools/Strelka'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :interval_list, :file, "interval list bed file", nil, :nofile => true
  extension :vcf
  task :strelka => :text do |tumor,normal,reference,interval_list|
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


  dep :strelka
  extension :vcf
  input :strelka_filter_tier, :boolean, "Take only tier 1 variants from strelka", true
  input :strelka_filter_evs, :integer, "Strelka EVS minimum score filter", 15
  input :strelka_filter_qss, :integer, "Strelka QSS minimum score filter", 15
  task :strelka_filtered => :text do |tier,evs,qss|
    TSV.traverse step(:strelka), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/

      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/
      qss_nt = line.split(";").select{|d| d =~ /^QSS_NT/}.first.split("=").last.to_i
      next unless qss_nt >= qss
      tqss_nt = line.split(";").select{|d| d =~ /^TQSS_NT/}.first.split("=").last.to_i
      next if tier and tqss_nt > 1 
      tqss = line.split(";").select{|d| d =~ /^TQSS/}.first.split("=").last.to_i
      next if tier and tqss > 1 
      somatic_evs = line.split(";").select{|d| d =~ /^SomaticEVS/}.first.split.first.split("=").last.to_f
      next unless somatic_evs >= evs

      line
    end
  end
end
