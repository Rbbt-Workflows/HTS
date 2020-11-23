require 'rbbt-util'
require 'rbbt/workflow'
require_relative '../tools/samtools'

module Strelka
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.Strelka_bin, :proc do |filename|
    url = "https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2"


    TmpFile.with_file do |tmpfile|
      Open.mkdir tmpfile
      CMD.cmd_log("wget '#{url}' -O - |tar xvjf - -C #{tmpfile}")
      codedir = Dir.glob(File.join(tmpfile, "*")).first
      Open.mv codedir, filename
    end
  end

  Rbbt.claim Rbbt.software.opt.Strelka, :install, Rbbt.share.install.software.Strelka.find

  CMD.tool "configureStrelkaSomaticWorkflow.py", Rbbt.software.opt.Strelka

  def self.runSomatic(tumor, normal, reference, output, cpus, interval_list)
    #cmd_config = Rbbt.software.opt.Strelka.produce.bin["configureStrelkaSomaticWorkflow.py"].find 
    cmd_config = "configureStrelkaSomaticWorkflow.py"

    cmd_string = "--tumorBam='#{tumor}' --ref='#{reference}' --runDir='#{output}' "
    cmd_string += " --normalBam='#{normal}' " unless normal.nil?
    cmd_string += " --callRegions='#{interval_list}' " unless interval_list.nil?

    CMD.cmd_log(cmd_config, cmd_string)    
    cmd_workflow = File.join(output, "runWorkflow.py")

    CMD.cmd_log("'#{ cmd_workflow }' --mode local -j #{cpus}")
  end

end

if __FILE__ == $0
  Log.severity = 0 

  normal = "~/.rbbt/software/opt/Strelka_bin/data/normal.bam"
  tumor = "~/.rbbt/software/opt/Strelka_bin/data/tumor.bam"
  reference = "/data/rbbt/share/databases/BWA/references/b37/reference.fa"
  TmpFile.with_file do |tmpfile|
    Strelka.runSomatic(normal, tumor, reference, tmpfile)
  end

end


