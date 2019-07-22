require 'tools/BAMSurgeon'
module HTS

  input :bed, :file, "variants file in bed format", :nofile => true
  input :bam, :file, "BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :snv_frac, :float, "Maximum allowable linked SNP MAF (for avoiding haplotypes)", 1
  input :mut_frac, :float, "Allelic fraction at which to make SNVs", 0.5
  input :num_snvs, :float, "maximum number of mutations to try (default: entire input)", nil
  input :cnv_file, :file, "tabix-indexed list of genome-wide absolute copy number values", nil, :nofile => true
  input :cover_diff, :float, "allow difference in input and output coverage", 0.1
  input :procs, :float, "split into multiple processes", 1
  input :picard_jar, :file, "path to picard.jar", nil, :nofile => true
  input :min_depth, :float, "minimum read depth to make mutation", 10
  input :max_depth, :float, "maximum read depth to make mutation", 2000
  input :min_mut_threads, :float , "minimum number of mutated reads to output per site", nil
  input :avoid_reads, :file, "file of read names to avoid (mutations will be skipped if overlap)", nil, :nofile => true
  input :aligner, :select, "aligner", "mem", :select_options => %w(mem backtrack novoalign gsnap STAR bowtie2 tmap bwakit), :nofile => true

  task :BAMSurgeon_add_indels => :binary do |bed,bam,reference,snv_frac,mut_frac,num_snvs,cnv_file,cover_diff,procs, picard_jar,min_depth,max_depth,min_mut_threads,avoid_reads,aligner|
    output = self.tmp_path
    reference = reference_file reference
    bam = Samtools.prepare_BAM(bam) if bam
    reference = Samtools.prepare_FASTA(reference)
    BAMSurgeon.add_indels(bed,bam,reference,output,snv_frac,mut_frac,num_snvs,cnv_file,cover_diff,procs,picard_jar,min_depth,max_depth,min_mut_threads,avoid_reads,aligner)  
  end

end
