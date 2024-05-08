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
    Misc.lock linked do
      if ! File.exist?(linked + ".salmon_idx") || Persist.newer?(linked + '.salmon_idx', file)

        Misc.in_dir dir do
          FileUtils.ln_s file, dir[basename] unless File.exist?(linked)
          CMD.cmd_log("salmon index -t '#{linked}' -i #{linked}.salmon_idx")
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
  task :salmon => :tsv do |fastq1,fastq2,organism|
    cdna = Salmon.prepare_CDNA_FASTA Organism.cdna_fasta(organism).produce.find
    Open.mkdir files_dir
    output = file('output')
    cpus = config :cpus, :salmon, :salmon_quant
    #CMD.cmd_log("salmon", "quant -l A --validateMappings --index #{cdna}.salmon_idx -1 #{fastq1 * " "} -2 #{fastq2 * " "} --output #{output}", "--threads" => cpus )
    #CMD.cmd_log("salmon", "quant -l A --index #{cdna}.salmon_idx -1 #{fastq1 * " "} -2 #{fastq2 * " "} --output #{output} --validateMappings ", "--threads" => cpus )
    CMD.cmd_log("salmon", "quant -l A -1 #{fastq1 * " "} -2 #{fastq2 * " "}",
                "no-version-check" => true,
                "output" => output,
                "threads" => cpus,
                "index" => "#{cdna}.salmon_idx",
                :add_option_dashes => true)
    Open.cp output["quant.sf"], self.tmp_path
    nil
  end
end
