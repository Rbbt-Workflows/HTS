module HTS 

  input :input, :file, "Input folder", nil, :nofile => true
  input :output, :file, "Output folder", nil, :nofile => true
  input :resume, :bool, "Resume previous execution or start a new one", true
  input :nextflow_script, :file,"Location of intogen nextflow script", nil, :nofile => true
  input :intogen_reference, :select, "Reference code", "hg38", :select_options => %w(hg38 hg19), :nofile => true
  task :intogen_pre => :binary do |input, output, resume, nextflow_script, intogen_reference|
    cmd_str ="export INTOGEN_VEP=vep92; export INTOGEN_RELEASE=v20191009; export INTOGEN_HOME=#{File.dirname(nextflow_script)}; export INTOGEN_GENOME=#{intogen_reference};  export LC_ALL=C.UTF-8;export LANG=C.UTF-8 "

    cmd_str.concat("~/nextflow #{nextflow_script}")
    
    if output.nil?
      output = self.path
    end
    args = {}
    args["--input"] = input
    args["--output"] = output   
    if resume
      cmd_str.concat(" --resume")
    end
    CMD.cmd_log(cmd_str,args)
    nil
  end 
end

