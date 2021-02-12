require 'rbbt-util'
require 'rbbt/workflow'
module Delly
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.Delly, :install, Rbbt.share.install.software.Delly.find

  CMD.tool :delly, Rbbt.software.opt.Delly

  def self.exclude(reference)
    case 
    when reference =~ /hg19|b37/
      Rbbt.software.opt.Delly.excludeTemplates["human.hg19.excl.tsv"].find
    when reference =~ /hg38/
      Rbbt.software.opt.Delly.excludeTemplates["human.hg38.excl.tsv"].find
    end
  end

  def self.call(tumor, normal, reference, exclude, output)
    #delly_cmd = Rbbt.software.opt.Delly.bin.delly.find
    delly_cmd = :delly

    if normal.nil?
      CMD.cmd_log(delly_cmd, "call -x #{exclude} -o #{output}.bcf -g #{reference} #{tumor}")
      CMD.cmd_log(:bcftools, "view #{output}.bcf > #{output}.vcf")
    else
      CMD.cmd_log(delly_cmd, "call -x #{exclude} -o #{output}.call -g #{reference} #{tumor} #{normal}")
      tumor_name = GATK.BAM_sample_name(tumor)
      normal_name = GATK.BAM_sample_name(normal)
      TmpFile.with_file do |sample_file|
        txt = [[tumor_name, 'tumor'], [normal_name, 'control']].collect{|p| p * "\t"} * "\n"
        Open.write(sample_file, txt)
        CMD.cmd_log(delly_cmd, "filter -f somatic -o #{output}.pre -s #{sample_file} #{output}.call ")
        CMD.cmd_log(delly_cmd, "call -g #{reference} -v #{output}.pre -o #{output}.bcf -x #{exclude} #{tumor} #{normal}")
        Misc.sensiblewrite(output + '.vcf', HTS.add_vcf_sample_header(CMD.cmd(:bcftools, "view #{output}.bcf", :pipe => true), tumor_name, normal_name))
      end
    end
  end

end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.Delly.produce(true).path
end


