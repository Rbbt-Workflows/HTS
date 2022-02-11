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
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  dep :GC_windows, :jobname => :reference
  extension 'seqz.gz'
  task :seqz => :text do |normal,tumor,reference|
    Step.wait_for_jobs dependencies

    orig_reference = reference_file(reference)
    reference = Samtools.prepare_FASTA orig_reference

    cpus = config :cpus, :sequenza_bam2seqz, :bam2seqz, :sequenza, :seqz
    if cpus && cpus.to_i > 1
      chromosomes = Open.read(reference.replace_extension('dict', true)).split("\n").select{|l| 
        (m = l.match(/LN:(\d+)/)) && m[1].to_i > 1_000_000
      }.collect{|l| l.match(/SN:([^\s]+)/)[1]}
        
      normal = Samtools.prepare_BAM normal
      tumor = Samtools.prepare_BAM tumor
      Open.mkdir files_dir
      output = file('output')
      CMD.cmd("sequenza-utils", "bam2seqz -gc '#{step(:GC_windows).path}' --chromosome #{chromosomes.collect{|chr| "'#{chr}'" } * " "} --parallel #{cpus.to_i}  -n '#{normal}' -t '#{tumor}' --fasta '#{reference}' -o #{output}")

      joined = file('joined')

      first = true
      chromosomes.each do |chr|
        file = output + "_" + chr
        if first
          CMD.cmd("cat #{file} > #{joined}")
        else
          CMD.cmd("tail -n +2 #{file} >> #{joined}")
        end
        first = false
      end

      CMD.cmd("bgzip #{joined}")

      FileUtils.mv joined + '.gz', self.tmp_path

      FileUtils.rm_rf files_dir
      nil
    else
      CMD.cmd("sequenza-utils", "bam2seqz -gc '#{step(:GC_windows).path}'  -n '#{normal}' -t '#{tumor}' --fasta '#{reference}' | bgzip > #{self.tmp_path}")
    end
    nil
  end

  dep :seqz
  task :sequenza => :array do

    reference = self.recursive_inputs[:reference]
    reference = reference_file(reference)
    reference = Samtools.prepare_FASTA reference
    chromosomes = reference.replace_extension('dict', true).read.split("\n")
      .select{|line| line.split("\t")[1].include?('SN:') }
      .collect{|line| line.split("\t")[1].sub('SN:', '') }

    chromosomes -= ["MT", "M"]

    Open.mkdir self.files_dir
    seq_binned = file('seqz.gz')
    CMD.cmd("'#{SEQUENZA_UTILS}' seqz_binning -w 50 -s '#{step(:seqz).path}' -o - | grep -v ^GL | gzip > '#{seq_binned}'")
            

    output_dir = file('output')
    R.run <<-EOF
rbbt.require('sequenza')
data.file = '#{seq_binned}'

chromosomes = #{R.ruby2R chromosomes}
data <- sequenza.extract(data.file, chromosome.list=chromosomes)
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
