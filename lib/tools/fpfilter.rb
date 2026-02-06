require 'rbbt-util'
require 'rbbt/workflow'
module FPFilter
  extend Resource
  
  Rbbt.claim Rbbt.software.opt.fpfilter, :install, Rbbt.share.install.software.fpfilter.find

  CMD.tool "bam-readcount", nil, "bash -c 'type bam-readcount'" do
    CMD.cmd('conda install bam-readcount -c bioconda')
  end

  def self.filter(args)
    begin
      fpfilter_cmd = Rbbt.software.opt.fpfilter.produce["fpfilter.pl"].find
      CMD.get_tool 'bam-readcount'
      CMD.cmd_log("perl #{fpfilter_cmd}", args.to_hash.merge('add_option_dashes' => true))
    rescue
      TmpFile.with_file nil, extension: :sam do |header|
        CMD.cmd_log("samtools view -H --no-PG #{args['bam-file']} |  sed 's/\\t\\(DT\\|....\\):[^\\t]*//g' > #{header}")
        TmpFile.with_file nil, extension: :bam do |tmp_file|
          CMD.cmd_log("samtools reheader #{header} #{args['bam-file']} > #{tmp_file}")
          args['bam-file'] = tmp_file
          CMD.cmd_log("perl #{fpfilter_cmd}", args.to_hash.merge('add_option_dashes' => true))
        end
      end
    end
  end
end

if __FILE__ == $0
  Log.severity = 0
end


