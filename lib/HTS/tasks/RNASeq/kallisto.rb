module Kallisto
  def self.prepare_CDNA_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.digest(Open.realpath(file))
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    Misc.lock linked do

      if ! File.exists?(linked + ".kallisto_idx") || Persist.newer?(linked + '.kallisto_idx', file)
        Misc.in_dir dir do
          FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
          CMD.cmd_log("kallisto index '#{linked}' -i #{linked}.kallisto_idx")
        end
      end

    end

    linked
  end

end

module HTS
  input :fastq1, :array, "FASTQ 1 files", nil, :nofile => true
  input :fastq2, :array, "FASTQ 2 files", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :kallisto => :tsv do |fastq1,fastq2,organism|
    cdna = Kallisto.prepare_CDNA_FASTA Organism.cdna_fasta(organism).produce.find
    Open.mkdir files_dir
    output = file('output')
    cpus = config :cpus, :kallisto, :kallisto_quant
    CMD.cmd_log("kallisto", "quant --plaintext --index #{cdna}.idx #{fastq1.zip(fastq2).flatten * " "} --output #{output}", "--threads" => cpus )
    Open.cp output["abundance.tsv"], self.tmp_path
    nil
  end

end
