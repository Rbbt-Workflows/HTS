require 'tools/BAMSurgeon'
module HTS

  input :varfile, :file, "variants file in bed format", :nofile => true
  input :bamfile, :file, "BAM", :nofile => true
  input :outbam, :file, "Output BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :snvfrac, :float, "Maximum allowable linked SNP MAF (for avoiding haplotypes)", 1
  input :mutfrac, :float, "Allelic fraction at which to make SNVs", 0.5
  input :numsnvs, :float, "maximum number of mutations to try (default: entire input)", nil
  input :cnvfile, :file, "tabix-indexed list of genome-wide absolute copy number values", nil, :nofile => true
  input :coverdiff, :float, "allow difference in input and output coverage", 0.1
  input :procs, :float, "split into multiple processes", 1
  input :picardjar, :file, "path to picard.jar", nil, :nofile => true
  input :mindepth, :string, "minimum read depth to make mutation", "10"
  input :maxdepth, :float, "maximum read depth to make mutation", 2000
  input :minmutreads, :float , "minimum number of mutated reads to output per site", nil
  input :avoidreads, :file, "file of read names to avoid (mutations will be skipped if overlap)", nil, :nofile => true
  input :aligner, :select, "aligner", "mem", :select_options => %w(mem backtrack novoalign gsnap STAR bowtie2 tmap bwakit), :nofile => true

  task :BAMSurgeon_add_indels => :binary do |varfile,bamfile,outbam,reference,snvfrac,mutfrac,numsnvs,cnvfile,coverdiff,procs,picardjar,mindepth,maxdepth,minmutreads,avoidreads,aligner|
    if outbam.nil?
      inputs[:outbam]= self.tmp_path
    end
    reference = reference_file reference
    bam = Samtools.prepare_BAM(bam) if bam
    reference = Samtools.prepare_FASTA(reference)
    inputs[:reference] = reference
    iii inputs.to_hash
    BAMSurgeon.add_indels(inputs)
  end

  input :varfile, :file, "variants file in bed format", nil, :nofile => true
  input :bamfile, :file, "BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :snvfrac, :float, "Maximum allowable linked SNP MAF (for avoiding haplotypes)", 1
  input :mutfrac, :float, "Allelic fraction at which to make SNVs", 0.5
  input :numsnvs, :float, "maximum number of mutations to try (default: entire input)", nil
  input :cnvfile, :file, "tabix-indexed list of genome-wide absolute copy number values", nil, :nofile => true
  input :coverdiff, :float, "allow difference in input and output coverage", 0.1
  input :haplosize, :float, "haplotype size", 0
  input :picardjar, :file, "path to picard.jar", nil, :nofile => true
  input :mindepth, :integer, "minimum read depth to make mutation", 10
  input :maxdepth, :float, "maximum read depth to make mutation", 2000
  input :minmutreads, :float , "minimum number of mutated reads to output per site", nil
  input :avoidreads, :file, "file of read names to avoid (mutations will be skipped if overlap)", nil, :nofile => true
  input :aligner, :select, "aligner", "mem", :select_options => %w(mem backtrack novoalign gsnap STAR bowtie2 tmap bwakit), :nofile => true
  task :BAMSurgeon_add_snvs => :binary do |varfile,bamfile,reference,snvfrac,mutfrac,numsnvs,cnvfile,coverdiff,haplosize,picardjar,mindepth,maxdepth,minmutreads,avoidreads,aligner|
    reference = reference_file reference
    reference = Samtools.prepare_FASTA(reference)

    i = inputs.to_hash

    i[:outbam]= self.tmp_path
    i[:bamfile] = Samtools.prepare_BAM(bamfile) 
    i[:reference] = reference
    i[:procs] = config :cpus, :BAMSurgeon, :bamsurgeon, :bam_surgeon, :default => 1

    BAMSurgeon.add_snvs(i)
    nil
  end


  input :varfile, :file, "variants file in bed format", :nofile => true
  input :bamfile, :file, "BAM", :nofile => true
  input :outbam, :file, "Output BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :maxlibsize, :float, "imaximum fragment length of seq. library"
  input :kmer, :float, "kmer size for assembly", 31
  input :svfrac, :float, "allele fraction of variant", 1.0
  input :require_exact, :boolean, "drop mutation if breakpoints cannot be made exactly as input", false
  input :minctglen, :float, "minimum length for contig generation, also used to pad assembly", 4000
  input :maxmuts, :float, "maximum number of mutations to make ", nil
  input :cnvfile, :file, "tabix-indexed list of genome-wide absolute copy number values", nil, :nofile => true
  input :aligner, :select, "aligner", "mem", :select_options => %w(mem backtrack novoalign gsnap STAR bowtie2 tmap bwakit), :nofile => true

  task :BAMSurgeon_add_struct_vars => :binary do |varfile,bamfile,outbam,reference,maxlibsize,kmer,svfrac,require_exact,minctglen,maxmuts,cnvfile,aligner|
    if outbam.nil?
      inputs[:outbam]= self.tmp_path
    end
    reference = reference_file reference
    bam = Samtools.prepare_BAM(bam) if bam
    reference = Samtools.prepare_FASTA(reference)
    inputs[:reference] = reference
    BAMSurgeon.add_struct_vars(inputs)
  end
end
