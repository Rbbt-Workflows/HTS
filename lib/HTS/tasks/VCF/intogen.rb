module HTS 

  input :input, :file, "Input folder", nil, :nofile => true
  input :resume, :bool, "Resume previous execution or start a new one", true
  input :nextflow_script, :file,"Location of intogen nextflow script", nil, :nofile => true
  input :intogen_reference, :select, "Reference code", "hg38", :select_options => %w(hg38 hg19), :nofile => true
  input :conda_env, :file, "Condat environment path", nil, :nofile => true
  task :intogen_pre => :text do |input, resume, nextflow_script, intogen_reference, conda_env|
    FileUtils.mkdir_p(self.path + ".files/work/conda")
    files_dir= self.path + ".files"
    nextflow_script = File.expand_path(nextflow_script)
    conda_env = File.expand_path(conda_env)
    input= File.expand_path(input)
    Dir.chdir(files_dir) do
      FileUtils.mkdir_p("work/conda/")
      FileUtils.mkdir_p(self.path + ".results/")
      begin
        FileUtils.ln_s  "#{conda_env}", files_dir + "/work/conda" unless Dir.exist?files_dir + "/work/conda/#{File.basename(conda_env)}"
      rescue
      end
      cmd_str ="export INTOGEN_VEP=vep92; export INTOGEN_RELEASE=v20191009; export INTOGEN_HOME=#{File.dirname(nextflow_script)}; export INTOGEN_GENOME=#{intogen_reference}; export LC_ALL=C.UTF-8; export LANG=C.UTF-8; "
      cmd_str.concat("~/nextflow run $(realpath #{nextflow_script})")
      args = {}
      args["--input"] = input
      args["--output"] = self.path + ".results"
      if resume
        cmd_str.concat(" --resume")
      end
      CMD.cmd_log(cmd_str, args)
    end 
    FileUtils.ln_s Dir.glob(self.path + ".results/oncodrivefml/*.out.gz")[0], self.path
    nil
  end
  
  dep :intogen_pre
  task :intogen_combination => :text do 
    intogen_pre_path=step(:intogen_pre).path
    inputs = step(:intogen_pre).inputs
    container = File.expand_path(Dir.glob(File.dirname(inputs[:nextflow_script]) + "/containers/#{inputs[:intogen_reference]}*/combination.simg")[0])
    output_path = intogen_pre_path + ".results/combination"
    FileUtils.mkdir_p(output_path)
    run_identifier = File.basename(Dir.glob(intogen_pre_path + ".results/signature/*.in.gz")[0],".in.gz")
    datasets_dir = Dir.glob(File.dirname(File.expand_path(inputs[:nextflow_script]))+"/datasets/#{inputs[:intogen_reference]}*").first
    containers_dir =  Dir.glob(File.dirname(File.expand_path(inputs[:nextflow_script]))+"/containers/#{inputs[:intogen_reference]}*").first
    Dir.chdir(output_path) do
      cmd_str = "export INTOGEN_DATASETS=#{datasets_dir}; export INTOGEN_CONTAINERS=#{containers_dir}; "
      cmd_str += "#{container} #{output_path} #{run_identifier}"
      CMD.cmd_log(cmd_str)
    FileUtils.mkdir_p(self.path+".files")
    FileUtils.cp Dir.glob(output_path+"/*.gz"), self.path+".files/"
    end
    FileUtils.ln_s output_path + "/#{run_identifier}.stouffer.out.gz", self.path
    nil
  end
end
