require 'rbbt-util'
require 'rbbt/workflow'
module SvABA
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.SvABA, :install, Rbbt.share.install.software.SvABA.find

  CMD.tool :svaba, Rbbt.software.opt.SvABA

  def self.call(tumor, normal, reference, svABA_options, output, cpus = nil)

    cpus ||= Rbbt::Config.get :cpus, :svaba, :svABA, :default => 1
    Misc.in_dir output do
      CMD.cmd_log(:svaba, "run -n #{ normal } -t #{ tumor } -G '#{reference}' -p #{cpus} -a svABA #{svABA_options}")
    end
  end

  def self.fix_vcf_sample_names(file, tumor_bam, normal_bam)
    fields = TSV.parse_header(file).fields
    tumor = fields.pop
    normal = fields.pop
    tumor_sample = GATK.BAM_sample_name tumor_bam
    normal_sample = GATK.BAM_sample_name normal_bam
    TSV.traverse file, :type => :array, :into => :stream do |line|
      next line unless line =~ /^(?:#|CHR)/
      line.sub! tumor, tumor_sample
      line.sub! normal, normal_sample
      if line =~ /^#?CHR/
        "##tumor_sample=#{tumor_sample}" + "\n" +
          "##normal_sample=#{normal_sample}" + "\n" +
          line
      else
        line
      end
    end
  end
end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.SvABA.produce(true).path
end


