module HTS 

  input :input, :file, "Input folder", nil, :nofile => true
  input :output, :file, "Output folder", nil, :nofile => true
  input :resume, :bool, "Resume previous execution or start a new one", true
  input :nextflow_script, :file,"Location of intogen nextflow script", nil, :nofile => true
  task :intogen_pre => :binary do |input, output, resume, nextflow_script|
    cmd_str="~/nextflow #{nextflow_script}"
    if output.nil?
      output = self.path
    end
    args = {}
    args["--input"] = input
    args["--output"] = output   
    if !resume.nil?
      cmd_str.concat(" --resume")
    end
    CMD.cmd_log(cmd_str,args)
    nil
  end 
end

