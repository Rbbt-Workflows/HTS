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
    interval_list = "/dev/shm/prueba.bed.gz"
    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)

    cpus = config :cpus, :strelka, :default => 3
    Strelka.runSomatic(tumor, normal, reference, output, cpus, interval_list)
    
    Open.read(output.results.variants["somatic.snvs.vcf.gz"])
  end


  dep :strelka
  extension :vcf
  task :strelka_filtered => :text do
	TSV.traverse step(:strelka), :into => :stream, :type => :array do |line|
        next line if line[0] =~ /^#/

        chr = line.split("\t").first
        next unless chr =~ /^(chr)?[0-9MTXY]+$/
        next unless line =~ /PASS/
        qss_nt = line.split(";").select{|d| d =~ /^QSS_NT/}.first.split("=").last.to_i
        next unless qss_nt > 15
        tqss_nt = line.split(";").select{|d| d =~ /^TQSS_NT/}.first.split("=").last.to_i
        next unless tqss_nt == 1
        tqss = line.split(";").select{|d| d =~ /^TQSS/}.first.split("=").last.to_i
        next unless tqss == 1
        somatic_evs = line.split(";").select{|d| d =~ /^SomaticEVS/}.first.split.first.split("=").last.to_f
        next unless somatic_evs > 15
    
        line
	end
  end
end
