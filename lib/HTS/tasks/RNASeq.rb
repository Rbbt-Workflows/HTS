require 'tools/FuSeq'

module Salmon
  def self.prepare_CDNA_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.digest(Open.realpath(file))
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".idx") || Persist.newer?(linked + '.idx', file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd_log("salmon index -t '#{linked}' -i #{linked}.idx")
      end
    end

    linked
  end

end

module HTS


  input :fastq1, :array, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :array, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :reference, :select, "Reference code", nil, :select_options => %w(b37 hg38), :nofile => true
  extension :bam
  task :hisat => :binary do |fastq1, fastq2,organism,reference|
    cpus = config("cpus", :hisat_build, :hisat)
    samtools_cpus = config("cpus", :samtools_index, :samtools, :index, :default => nil)

    reference = 'b37' if reference.nil? && organism.nil?
    organism = Organism.organism_for_build(reference) if organism.nil? and reference
    reference = case Organism.hg_build(organism)
                when 'hg19'
                  'b37'
                when 'hg38'
                  'hg38'
                end

    reference = reference_file(reference) 

    reference = GATK.prepare_FASTA reference
    reference = Samtools.prepare_FASTA reference
    reference = HTS.uncompress_FASTA reference

    index = HISAT.build_gft_index(organism, reference, cpus)

    set_info :organism, organism

    sam = file('out.sam')
    CMD.cmd_log("hisat2 -p #{cpus || 1} --dta -x #{index} -1 #{fastq1*","} -2 #{fastq2*""} -S #{sam}")
    CMD.cmd_log("samtools -@ #{samtools_cpus || 1} -o #{self.tmp_path} #{sam} ")
    Open.rm sam
    nil
  end

  dep :hisat
  extension :gtf
  task :stringtie => :tsv do
    organism = self.recursive_inputs[:organism] || step(:hisat).info[:organism]
    cpus = config("cpus", :hisat_build, :hisat)
    gft_file = HTS.gtf_file(Organism.gene_set(organism).produce.find)
    CMD.cmd_log("stringtie -e -B -p #{cpus} -G #{gft_file} -o #{self.tmp_path} #{step(:hisat).path}")
    nil
  end

  dep :stringtie
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :htseq_counts => :tsv do |organism|
    CMD.cmd_log("htseq-count -f bam '#{step(:stringtie).path}' '#{Organism.gene_set(organism).produce.find}' > #{self.tmp_path}")
    nil
  end

  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :salmon => :tsv do |fastq1,fastq2,organism|
    cdna = Salmon.prepare_CDNA_FASTA Organism.cdna_fasta(organism).produce.find
    Open.mkdir files_dir
    output = file('output')
    CMD.cmd_log("salmon quant -l A --index #{cdna}.idx -1 #{fastq1} -2 #{fastq2} --output #{output}")
    Open.cp output["quant.sf"], self.tmp_path
    nil
  end


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
