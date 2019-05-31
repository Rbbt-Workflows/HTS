require 'rbbt-util'
require 'rbbt/workflow'
require_relative '../tools/samtools'
module ControlFREEC
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.ControlFREEC, :install, Rbbt.share.install.software.ControlFREEC.find
	
  FREEC_CMD = Rbbt.software.opt.ControlFREEC.src.freec.find
  FREEC_PLOT_SCRIPT = Rbbt.software.opt.ControlFREEC.scripts["makeGraph.R"].find
  def self.run(config_file)
    CMD.cmd_log("'#{FREEC_CMD}' -conf '#{config_file}'")
  end

  def self.makegraphs(results_path, output)
    Open.mkdir output
    CMD.cmd_log("ln -s '#{results_path}'/*_ratio.txt '#{results_path}'/*_BAF.txt #{output}")
    CMD.cmd_log("cat '#{FREEC_PLOT_SCRIPT}' | R --slave --args 2 '#{output}'/*_ratio.txt '#{output}'/*_BAF.txt")
  end
end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.ControlFREEC.produce(true).path
end


