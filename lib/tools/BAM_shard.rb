require 'tools/samtools'
require 'digest'

require 'pry'

class BAMShard
  DIGEST_SIZE = 128
  CHUNK_LIMIT = 100_000_000 
  def self.create_chunk_files(samfile)
    chunks ||= Rbbt::Config.get('cpus', 'bamshard', default: 3)
    chunk_size = 2 ** DIGEST_SIZE / chunks.to_i
    header = Samtools.header(samfile)
    i=0
    reads = {}
    File.foreach(samfile) do |fd|
      if fd.start_with?("@")
        next
      end
      read_name = fd.split()[0]
      chunk_selected = Digest::MD5.hexdigest(read_name).to_i(16)/chunk_size
        
      reads[chunk_selected] = [] if reads[chunk_selected].nil?
      reads[chunk_selected] << fd

      if reads[chunk_selected].length > CHUNK_LIMIT
        if File.zero?("chunk#{chunk_selected}.reads")
          File.open("chunk#{chunk_selected}.reads","w"){|f| f.write header}
        end
        File.open("chunk#{chunk_selected}.reads","w"){|f| reads[chunk_selected].each {|read| f.write read}}          
        reads[chunk_selected] = []
      end
    
    end
    reads.keys.each do |key|    
      if File.zero?("chunk#{key}.reads")
        File.open("chunk#{key}.reads","w"){|f| f.write header}
      end
      File.open("chunk#{key}.reads","w"){|f| reads[key].each {|read| f.write read}}
    end
  end

  def self.revert_BAM(bamfile, output, &callback)

    cpus ||= Rbbt::Config.get('cpus', 'bamshard', :default => 3)

    Open.mkdir output + ".files"
    TmpFile.with_file do |workdir|
      Open.mkdir workdir
      #Dir.chdir workdir

      Misc.in_dir workdir do
        Dir.mkdir "temp"
        outsam = workdir + "/temp/" + File.basename(bamfile) + ".sam"
        binding.pry
        Samtools.run("view -h #{bamfile} > #{outsam}")
        self.create_chunk_files(outsam)

      end
      binding.pry
      q = RbbtProcessQueue.new cpus
      q.callback &callback
      q.init do |bam|
          args = {}
          args["INPUT"] = workdir + "/" + bam
          args["OUTPUT"] = output + ".files" + "/" + "#{File.basename(bam)}.ubam"
          args["OUTPUT_BY_READGROUP"] = "false"
          args["SANITIZE"] = "true"
          args["SORT_ORDER"] = "queryname"
          args["ATTRIBUTE_TO_CLEAR"] = ["XA", "BD", "XS", "BI"]
          args["RESTORE_ORIGINAL_QUALITIES"] = "true"
          args["REMOVE_DUPLICATE_INFORMATION"] = "true"
          args["REMOVE_ALIGNMENT_INFORMATION"] = "true"
          args["VALIDATION_STRINGENCY"] = "SILENT"

          GATK.run_log("RevertSam", args)
      end
      bams =  Dir.glob("#{workdir}/*reads").each {|f| !File.directory? f}
      for bam in bams do
        iii bam
        binding.pry
        q.process bam
      end

      q.join
      binding.pry
      inbams = Dir.glob(output + ".files/*reads").map(&File.method(:realpath))
      args = {}
      args["OUTPUT"] = output
      args["INPUT"] = inbams
      args["USE_THREADING"] = "true"
      args["MERGE_SEQUENCE_DICTIONARIES"] = "true"
      args["SORT_ORDER"] = "queryname"

      GATK.run_log("MergeSamFiles", args)

      nil
    end
  end
end
