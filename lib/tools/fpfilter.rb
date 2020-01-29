require 'rbbt-util'
require 'rbbt/workflow'
module FPFilter
  extend Resource
  
  Rbbt.claim Rbbt.software.opt.fpfilter, :install, Rbbt.share.install.software.fpfilter.find

  CMD.tool "bam-readcount", nil, "bash -c 'type bam-readcount'" do
    CMD.cmd('conda install bam-readcount -c bioconda')
  end

  def self.filter(args)
    fpfilter_cmd = Rbbt.software.opt.fpfilter.produce["fpfilter.pl"].find
    CMD.get_tool 'bam-readcount'
    CMD.cmd_log("perl #{fpfilter_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end
end

if __FILE__ == $0
  Log.severity = 0
end


