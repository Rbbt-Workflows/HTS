require 'pathname'
module HTS

  input :input, :file, "Input folder with mafs and dataset.bginfo files", nil, :nofile => true
  input :intogen_reference, :select, "Reference code", "hg38", :select_options => %w(hg38 hg19), :nofile => true
  #resumable
  task :intogen_pre => :text do |input, intogen_reference|
    nextflow_script = config('NXFScript', :intogen)
    conda_env = ENV["CONDA_PREFIX"]
    if conda_env.nil?
      raise "An activated conda environment with intogen dependencies is necessary to continue the execution"
    end
    input = File.expand_path(input)

    Misc.in_dir file('output') do
      FileUtils.mkdir_p("work/conda/")
      FileUtils.ln_s  "#{conda_env}", file('output') + "/work/conda" unless Dir.exist?file('output') + "/work/conda/#{File.basename(conda_env)}"
      cmd_str ="export INTOGEN_VEP=vep92; export INTOGEN_RELEASE=v20191009; export INTOGEN_HOME=#{File.dirname(nextflow_script)}; \
      export INTOGEN_GENOME=#{intogen_reference}; export LC_ALL=C.UTF-8; export LANG=C.UTF-8; "
      cmd_str += "nextflow run  #{nextflow_script} -resume "
      args = {}
      args["--input"] = input
      args["--output"] = file("results")

      CMD.cmd_log(cmd_str, args)
    end
    FileUtils.touch file("results") + "/done"
    FileUtils.ln_s file("results") + "/done", self.path
    nil
  end

  dep :intogen_pre
  task :intogen_combination => :text do
    intogen_pre_path = step(:intogen_pre).path
    nextflow_script = config('NXFScript', :intogen)
    results_path = intogen_pre_path + ".files/results"
    output_path = results_path + "/combination"

    FileUtils.mkdir_p(output_path)
    datasets_path = Dir.glob(File.dirname(nextflow_script) + "/datasets/#{step(:intogen_pre).inputs[:intogen_reference]}*").first
    combination_container_path = Dir.glob(File.dirname(nextflow_script) + "/containers/#{step(:intogen_pre).inputs[:intogen_reference]}*/combination.simg").first
    run_identifier = File.basename(Dir.glob(results_path + "/signature/*.in.gz")[0],".in.gz")

    Misc.in_dir  output_path do
      cmd_str = "export INTOGEN_DATASETS=#{datasets_path}; export INTOGEN_CONTAINERS=#{Pathname.new(combination_container_path).dirname.to_s}; "
      cmd_str += "#{combination_container_path} #{output_path} #{run_identifier}"
      CMD.cmd_log(cmd_str)
    end
    FileUtils.ln_s output_path + "/#{run_identifier}.stouffer.out.gz", self.path
    nil
  end

end
