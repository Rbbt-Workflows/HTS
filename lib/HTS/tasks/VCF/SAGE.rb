module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :vcf
  task :sage_pre => :text do |tumor,normal,reference|
    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    options = {}
    options[:cpus] = config :cpus, :SAGE, :sage, :default => 3
    SAGE.run(tumor, normal, self.tmp_path, reference, options)
    
    nil
  end

  dep :sage_pre
  extension :vcf
  task :sage => :text do 
    TSV.traverse step(:sage_pre), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line.split("\t")[6].split("\W").include? "PASS"

      parts = line.split("\t")
      parts.last.sub!(/^\.\/\./,'1/0')

      parts * "\t"
    end

  end
end
