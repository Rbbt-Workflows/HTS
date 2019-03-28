module HTS

  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  extension :bam
  task :RNASeqBAM => :binary do |fastq1, fastq2|
    if fastq2
      args = "-1 #{fastq1} -2 #{fastq2}"
    else
      args = "-U #{fastq1}"
    end
    CMD.cmd_log("hisat -S #{self.tmp_path} #{args}")
    nil
  end

  dep :RNASeqBAM
  task :RNASeq => :tsv do
    CMD.cmd_log("stringtie #{step(:RNASeqBAM).path} -o #{self.tmp_path}")
    nil
  end

  dep :RNASeqBAM
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :htseq_counts => :tsv do |organism|
    CMD.cmd_log("htseq-count -f bam '#{step(:RNASeqBAM).path}' '#{Organism.gene_set(organism).produce.find}' > #{self.tmp_path}")
    nil
  end


end
