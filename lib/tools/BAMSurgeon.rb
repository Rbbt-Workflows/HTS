require 'rbbt-util'
require 'rbbt/workflow'
module BAMSurgeon
  extend Resource
  # ToDo pull files out of directory
  #
  
  Rbbt.claim Rbbt.software.opt.BAMSurgeon, :install, Rbbt.share.install.software.BAMSurgeon.find

  def self.add_indels(args)
    addindels_cmd = Rbbt.software.opt.BAMSurgeon.bin["addindel.py"].find
    CMD.cmd_log("#{addindels_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end

  def self.add_snvs(args)
    addindels_cmd = Rbbt.software.opt.BAMSurgeon.bin["addsnv.py"].find
    CMD.cmd_log("#{addindels_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end

  def self.add_struct_vars(args)
    addindels_cmd = Rbbt.software.opt.BAMSurgeon.bin["addsv.py"].find
    CMD.cmd_log("#{addindels_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end
end

if __FILE__ == $0
  Log.severity = 0 
  #iii Rbbt.software.opt.BAMSurgeon.produce(true)
end


