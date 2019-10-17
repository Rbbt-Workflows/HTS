require 'tools/samtools'

class BAMShard

  GAP_SIZE = 1_000
  CHUNK_SIZE = 10_000_000

  def self.create_chunk_files(bamfile)
    chunks ||= Rbbt::Config.get('cpus', 'bamshard', :default => 3)
    num_reads = Samtools.reads_number(bamfile)
    chunk_size = num_reads.to_i / (chunks.to_i - 1)
    header = Samtools.header(bamfile)
    Samtools.BAM_sort(bamfile, to_sam=true)
    bamfile = bamfile + ".sorted"
    chunk_idx = 0
    idx = 0
    chfd = open("chunk#{chunk_idx}.reads", "w")
    chfd << header
    previous_line = ""
    File.foreach(bamfile) do |fd|
      if fd.start_with?("@")
        next
      end
      if idx < chunk_size
        idx += 1
        previous_line = fd
        chfd << fd
        next
      else
        if  previous_line.split()[0] == fd.split()[0]
          chfd << fd
          previous_line = fd
        else
          idx = 0
          previous_line = fd
          chfd.close()
          chunk_idx += 1
          chfd = open("chunk#{chunk_idx}.reads","w")
          chfd << header
          chfd << fd
        end
      end
    end
  end

  def self.revert_BAM(bamfile, output, &callback)

    cpus ||= Rbbt::Config.get('cpus', 'bamshard', :default => 3)

    q = RbbtProcessQueue.new cpus
    q.callback &callback
    Open.mkdir output + ".files"
    TmpFile.with_file do |workdir|
      Open.mkdir workdir
      #Dir.chdir workdir

      Misc.in_dir workdir do
        self.create_chunk_files(bamfile)
      end

      q.init do |bam|
          args = {}
          args["INPUT"] = workdir + "/" + bam
          args["OUTPUT"] = output + ".files" + "/" + "#{bam}.ubam"
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

      bams =  Dir.entries(workdir).select {|f| !File.directory? f}
      for bam in bams do
        q.process bam
      end

      q.join
      inbams = Dir.glob(output + ".files/*").map(&File.method(:realpath)).join(" ")
      Samtools.merge(output,inbams)
      nil
    end
  end
end
