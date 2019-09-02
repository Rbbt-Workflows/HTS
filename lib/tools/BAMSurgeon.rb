require 'rbbt-util'
require 'rbbt/workflow'
module BAMSurgeon
  extend Resource
  # ToDo pull files out of directory
  #
  
  Rbbt.claim Rbbt.software.opt.Velvet, :install, "https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz"
  Rbbt.claim Rbbt.software.opt.BAMSurgeon, :install, Rbbt.share.install.software.BAMSurgeon.find
  Rbbt.claim Rbbt.software.opt.PicardTools, :install, Rbbt.share.install.software.PicardTools.find
  Rbbt.claim Rbbt.software.opt.Exonerate, :install do
    url = "https://github.com/adamewing/exonerate.git"
    commands =<<-EOF
get_git "$name" "$url"
cd $(opt_dir "$name")
automake --add-missing
build "$name" $extra
    EOF
    {:url => url, :commands => commands}
  end

  def self.add_indels(args)
    addindels_cmd = Rbbt.software.opt.BAMSurgeon.produce.bin["addindel.py"].find
    CMD.cmd_log("python2 #{addindels_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end

  def self.add_snvs(args)
    addsnvs_cmd = Rbbt.software.opt.BAMSurgeon.produce.bin["addsnv.py"].find
    CMD.cmd_log("python2 #{addsnvs_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end

  def self.add_struct_vars(args)
    addsvs_cmd = Rbbt.software.opt.BAMSurgeon.produce.bin["addsv.py"].find
    CMD.cmd_log("python2 #{addsvs_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end
end

if __FILE__ == $0
  Log.severity = 0 
  #iii Rbbt.software.opt.BAMSurgeon.produce(true)
  iii Rbbt.software.opt.Exonerate.produce(true)
end


