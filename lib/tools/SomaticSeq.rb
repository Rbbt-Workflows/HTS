require 'rbbt-util'
require 'rbbt/workflow'

module SomaticSeq

  Rbbt.claim Rbbt.software.opt.SomaticSeq, :install, Rbbt.share.install.software.SomaticSeq.find
  SOMATICSEQ_CMD = Rbbt.software.opt.SomaticSeq.produce["somaticseq_parallel.py"].find

  def self.run(args)
    args = args.to_hash.merge('add_option_dashes' => true)
    cmd_string = "python3 '#{SOMATICSEQ_CMD}' -ref #{args["genome-reference"]} -outdir #{args["outdir"]} paired "
    args.delete("genome-reference")
    args.delete("outdir")
    CMD.cmd_log(cmd_string, args.to_hash.merge('add_option_dashes' => true))
  end
end

if __FILE__ == $0
  Log.severity = 0
  #iii Rbbt.software.opt.SomaticSeq.produce(true).path
end

