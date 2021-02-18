require 'rbbt/util/R'
module HTS

  SEQUENZA_UTILS = 'sequenza-utils'

  CMD.tool "sequenza-utils" do
    CMD.cmd('pip install sequenza-utils')
  end
  
  extension 'gc.gz'
  task :GC_windows => :text do 
    orig_reference = reference_file(self.clean_name)
    reference = Samtools.prepare_FASTA orig_reference
    CMD.cmd("'#{SEQUENZA_UTILS}' gc_wiggle -w 50 -f '#{reference}' -o - | gzip > #{self.tmp_path}")
    nil
  end

  #input :normal, :file, "Normal BAM", nil, :nofile => true
  #input :tumor, :file, "Tumor BAM", nil, :nofile => true
  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  #dep :pileup, :bam => :placeholder do |jobname,options|
  #  deps = []

  #  deps << {:inputs => options.merge(:bam => options[:normal]), :jobname => jobname + '.normal'}
  #  deps << {:inputs => options.merge(:bam => options[:tumor]), :jobname => jobname + '.tumor'}

  #  deps
  #end
  #dep :GC_windows, :jobname => :reference
  #extension 'seqz.gz'
  #task :seqz => :text do |normal,tumor,reference|
  #  Step.wait_for_jobs dependencies
  #  normal = dependencies.select{|dep| dep.task_name.to_s == 'pileup'}.select{|dep| dep.clean_name =~ /\.normal/}.first
  #  tumor = dependencies.select{|dep| dep.task_name.to_s == 'pileup'}.select{|dep| dep.clean_name =~ /\.tumor/}.first

  #  orig_reference = reference_file(reference)
  #  reference = Samtools.prepare_FASTA orig_reference

  #  CMD.cmd("sequenza-utils", "bam2seqz -gc '#{step(:GC_windows).path}' -p -n '#{normal.path}' -t '#{tumor.path}' --fasta '#{reference}' | gzip > #{self.tmp_path}")
  #  nil
  #end

  input :normal, :file, "Normal BAM", nil, :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  dep :GC_windows, :jobname => :reference
  extension 'seqz.gz'
  task :seqz => :text do |normal,tumor,reference|
    Step.wait_for_jobs dependencies

    orig_reference = reference_file(reference)
    reference = Samtools.prepare_FASTA orig_reference

    CMD.cmd("sequenza-utils", "bam2seqz -gc '#{step(:GC_windows).path}'  -n '#{normal}' -t '#{tumor}' --fasta '#{reference}' | bgzip > #{self.tmp_path}")
    nil
  end

  dep :seqz
  task :sequenza => :array do

    Open.mkdir self.files_dir
    seq_binned = file('seqz.gz')
    CMD.cmd("'#{SEQUENZA_UTILS}' seqz_binning -w 50 -s '#{step(:seqz).path}' -o - | grep -v ^GL | gzip > '#{seq_binned}'")
            
    output_dir = file('output')
    R.run <<-EOF
rbbt.require('sequenza')
data.file = '#{seq_binned}'

data <- sequenza.extract(data.file)
CP.data <- sequenza.fit(data)

sequenza.results(sequenza.extract = data, cp.table = CP.data, sample.id = "#{self.clean_name}", out.dir="#{output_dir}")
    EOF

    Open.rm seq_binned
    output_dir.glob("*")
  end

  dep :sequenza
  task :sequenza_purity => :float do
    step(:sequenza).file('output').glob('*_confints_CP.txt').first.read.split("\n")[2].split("\t")[0]
  end

  dep :sequenza
  task :sequenza_ploidy => :float do
    step(:sequenza).file('output').glob('*_confints_CP.txt').first.read.split("\n")[2].split("\t")[1]
  end

  dep :sequenza
  task :sequenza_CNV => :tsv do
    parser = TSV::Parser.new step(:sequenza).file('output').glob('*_segments.txt').first, :header_hash => "", :type => :list
    dumper = TSV::Dumper.new :key_field => "Chromosome Range", :fields => parser.fields[2..-1], :type => :list
    dumper.init
    TSV.traverse parser, :into => dumper do |chr, values|
      start, eend, *rest = values
      [[chr, start, eend] * ":", rest]
    end
  end
end
