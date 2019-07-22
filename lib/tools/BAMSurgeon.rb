require 'rbbt-util'
require 'rbbt/workflow'
module BAMSurgeon
  extend Resource
  # ToDo pull files out of directory
  #
  
  Rbbt.claim Rbbt.software.opt.BAMSurgeon, :install, Rbbt.share.install.software.BAMSurgeon.find

  def self.add_indels(bed,bam,reference,output,snv_frac,mut_frac,num_snvs,cnv_file,cover_diff,procs,picard_jar,min_depth,max_depth,min_mut_threads,avoid_reads,aligner)
    addindels_cmd = Rbbt.software.opt.BAMSurgeon.bin["addindel.py"].find

    CMD.cmd_log("#{addindels_cmd} -v #{bed} -f #{bam} -r #{reference} -o #{output} -s #{snv_frac} -m #{mut_frac} -n #{num_snvs} -c #{cnv_file} -d #{cover_diff} -p #{procs} --picardjar #{picard_jar} --mindepth #{min_depth} --maxdepth #{max_depth}  --minmutreads #{min_mut_threads} --avoidreads #{avoid_reads} --aligner #{aligner}")
  end


end

if __FILE__ == $0
  Log.severity = 0 
  #iii Rbbt.software.opt.BAMSurgeon.produce(true)
end


