require 'tools/samtools'

class BAMShard

  GAP_SIZE = 1_000
  CHUNK_SIZE = 10_000_000

  def self.create_chr_bams(bamfile)
    chr_list = (1..22).to_a.map(&:to_s) + ["X","Y"]
    for chr in chr_list do
      Samtools.BAM_get_chr_reads(bamfile, chr)
    end
  end

  def self.cmd(bamfile, output, &callback)

    cpus ||= Rbbt::Config.get('cpus', 'bamshard', :default => 3)

    q = RbbtProcessQueue.new cpus
    q.callback &callback
    Open.mkdir output + ".files"
    TmpFile.with_file do |workdir|
      Open.mkdir workdir
      #Dir.chdir workdir

      Misc.in_dir workdir do
        self.create_chr_bams(bamfile)
      end

      q.init do |bam|
          args = {}
          args["INPUT"] = workdir + "/" + bam
          args["OUTPUT"] = output + ".files" + "/" + "#{bam}.ubam"
          args["OUTPUT_BY_READGROUP"] = "false"
          args["SANITIZE"] = "false"
          args["SORT_ORDER"] = "queryname"
          args["ATTRIBUTE_TO_CLEAR"] = ["XA", "BD", "XS", "BI"]
          args["RESTORE_ORIGINAL_QUALITIES"] = "true"
          args["REMOVE_DUPLICATE_INFORMATION"] = "true"
          args["REMOVE_ALIGNMENT_INFORMATION"] = "true"

          GATK.run_log("RevertSam", args)
          "DONE"
      end

      bams =  Dir.entries(workdir).select {|f| !File.directory? f}
      for bam in bams do
        q.process bam
      end

      q.join

      bams = Dir.glob(output + ".files/*").map(&File.method(:realpath))
      args = {}
      args["OUTPUT"] = output
      args["INPUT"] = bams
      args["USE_THREADING"] = "true"
      args["MERGE_SEQUENCE_DICTIONARIES"] = "true"

      GATK.run_log("MergeSamFiles", args)
    end
  end
end
