require 'rbbt-util'
require 'rbbt/workflow'
module Delly
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.Delly, :install, Rbbt.share.install.software.Delly.find

  def self.exclude(reference)
    case 
    when reference =~ /hg19|b37/
      Rbbt.software.opt.Delly.excludeTemplates["human.hg19.excl.tsv"].find
    when reference =~ /hg38/
      Rbbt.software.opt.Delly.excludeTemplates["human.hg38.excl.tsv"].find
    end
  end

  def self.call(tumor, normal, reference, exclude, output)
    delly_cmd = Rbbt.software.opt.Delly.bin.delly.find

    if normal.nil?
      CMD.cmd_log("#{delly_cmd} call -x #{exclude} -o #{output}.bcf -g #{reference} #{tumor}")
      CMD.cmd_log("bcftools view #{output}.bcf > #{output}.vcf")
    else
      CMD.cmd_log("#{delly_cmd} call -x #{exclude} -o #{output}.call -g #{reference} #{tumor} #{normal}")
      tumor_name = GATK.BAM_sample_name(tumor)
      normal_name = GATK.BAM_sample_name(normal)
      TmpFile.with_file do |sample_file|
        txt = [[tumor_name, 'tumor'], [normal_name, 'normal']].collect{|p| p * "\t"} * "\n"
        Open.write(txt, sample_file)
        CMD.cmd_log("#{delly_cmd} filter -f somatic -o #{output}.pre -s #{sample_file} #{output}.call ")
      end
    end
  end

end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.Delly.produce(true).path
end


