require 'rbbt-util'
require 'rbbt/workflow'
module SvABA
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.SvABA, :install, Rbbt.share.install.software.SvABA.find

  def self.call(tumor, normal, reference, svABA_options, output)

    Misc.in_dir output do
      CMD.cmd_log("svaba run -t #{ tumor } -n #{ normal } -G '#{reference}' -a svABA #{svABA_options}")
    end
  end
end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.SvABA.produce(true).path
end


