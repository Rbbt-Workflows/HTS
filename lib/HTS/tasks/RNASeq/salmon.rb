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
  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :salmon => :tsv do |fastq1,fastq2,organism|
    cdna = Salmon.prepare_CDNA_FASTA Organism.cdna_fasta(organism).produce.find
    Open.mkdir files_dir
    output = file('output')
    cpus = config :cpus, :salmon, :salmon_quant
    CMD.cmd_log("salmon", "quant -l A --index #{cdna}.idx -1 #{fastq1} -2 #{fastq2} --output #{output}", "--threads" => cpus )
    Open.cp output["quant.sf"], self.tmp_path
    nil
  end
end
