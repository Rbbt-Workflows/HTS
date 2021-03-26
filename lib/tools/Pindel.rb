require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/util/cmd'
module Pindel

  #Rbbt.claim Rbbt.software.opt.Pindel, :install do
  #  {:git => "https://github.com/genome/pindel.git"}
  #end

  CMD.tool :pindel, nil, "bash -c 'type pindel'" do
    CMD.cmd('conda install pindel -c bioconda')
  end

  CMD.tool :pindel2vcf, nil, "bash -c 'type pindel2vcf'" do
    CMD.cmd('conda install pindel -c bioconda')
  end
    
  def self.call(tumor, normal, reference, insert_size, output, cpus = nil)

    reference = reference.remove_extension('.gz')
    cpus ||= Rbbt::Config.get :cpus, :pindel, :Pindel, :default => 1

    if insert_size.nil?
      insert_size = CMD.cmd(:samtools, "view #{tumor} | head -n 10000 | awk '{if ($9 > 0 && $9 < 2000) {S+=$9; T+=1}}END{print S/T}'").read.to_i
    end

    Misc.in_dir output do
      TmpFile.with_file do |bam_files|
        Open.write(bam_files) do |f|
          f.puts [tumor, insert_size, "tumor"] * " "
          f.puts [normal, insert_size, "normal"] * " "
        end
        CMD.cmd_log(:pindel, "-f #{reference} -i #{bam_files} -o #{output}/pindel -T #{cpus}")
      end
    end
  end
end


if __FILE__ == $0
  Log.severity = 0
  CMD.cmd_log(:pindel)
end
