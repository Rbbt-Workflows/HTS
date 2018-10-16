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


  def self.runSomatic(normal, tumor, reference, output, cpus = 3)
    cmd_config = Rbbt.software.opt.Strelka_bin.bin["configureStrelkaSomaticWorkflow.py"].find


    CMD.cmd_log("'#{ cmd_config }' --normalBam='#{normal}' --tumorBam='#{tumor}' --ref='#{reference}' --runDir='#{output}'")
    
    cmd_workflow = File.join(output, "runWorkflow.py")

    CMD.cmd_log("'#{ cmd_workflow }' --mode local -j #{cpus}")
  end

  #Rbbt.claim Rbbt.software.opt.Strelka, :install, Rbbt.share.install.software.Strelka.find
end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.Strelka_bin.produce(true).path

  normal = "~/.rbbt/software/opt/Strelka_bin/data/normal.bam"
  tumor = "~/.rbbt/software/opt/Strelka_bin/data/tumor.bam"
  reference = "/data/rbbt/share/databases/BWA/references/b37/reference.fa"
  TmpFile.with_file do |tmpfile|
    Strelka.runSomatic(normal, tumor, reference, tmpfile)
  end

end


