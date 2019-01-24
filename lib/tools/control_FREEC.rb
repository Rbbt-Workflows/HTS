require 'rbbt-util'
require 'rbbt/workflow'
require_relative '../tools/samtools'
module ControlFREEC
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.ControlFREEC, :install, Rbbt.share.install.software.ControlFREEC.find
	
	FREEC_CMD = Rbbt.software.opt.ControlFREEC.src.freec.find

  def self.run(config_file)
		CMD.cmd_log("'#{FREEC_CMD}' -conf '#{config_file}'")
  end
end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.ControlFREEC.produce(true).path
end


