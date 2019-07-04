
module HTS
  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :kmer, :integer, "K-mer size", nil
  task :FuSeq => :text do |fastq1,fastq2,organism,kmer|
    output = file('output')
    cpus = config(:cpus, :FuSeq, :fuseq, :default => 1) 
    FuSeq.run [fastq1, fastq2].compact, organism, kmer, cpus, output
    "Done"
  end

  dep :FuSeq
  task :FuSeq_process => :tsv do
    organism = step(:FuSeq).recursive_inputs[:organism]
    output = file('output')

    FuSeq.process step(:FuSeq).file('output'), organism, output

    TSV.open file('output/fusions.FuSeq'), :header_hash => '', :merge => true, :key_field => "fusionName"
  end
end
