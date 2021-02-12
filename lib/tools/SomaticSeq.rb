require 'rbbt-util'
require 'rbbt/workflow'

module SomaticSeq

  Rbbt.claim Rbbt.software.opt.SomaticSeq, :install, Rbbt.share.install.software.SomaticSeq.find
  SOMATICSEQ_CMD = Rbbt.software.opt.SomaticSeq["somaticseq_parallel.py"]

  CMD.tool :somaticseq, Rbbt.software.opt.SomaticSeq, "python3 #{SOMATICSEQ_CMD}" #, "python3 #{SOMATICSEQ_CMD}"
    
  def self.run(args)
    args = args.to_hash.merge('add_option_dashes' => true)
    cmd_string = "python3 '#{SOMATICSEQ_CMD}' -ref #{args["genome-reference"]} -outdir #{args["outdir"]} paired "
    args.delete("genome-reference")
    args.delete("outdir")
    CMD.cmd_log(:somaticseq, args.to_hash.merge('add_option_dashes' => true))
  end
end
