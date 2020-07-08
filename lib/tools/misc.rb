module HTS
  Rbbt.claim Rbbt.software.opt.BAM_readcount, :install, "https://github.com/genome/bam-readcount.git"
  CMD.tool "bam-readcount", Rbbt.software.opt.BAM_readcount

  def self.prepare_BED(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    CMD.cmd_log("bgzip -c #{file} > '#{ linked + '.gz' }'") unless File.exists?(linked + '.gz')
    if ! File.exists?(linked + ".gz.tbi") || Persist.newer?(linked + ".gz.tbi", file)
      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd_log("tabix '#{ linked }.gz'")
      end
    end

    linked + '.gz'
  end
end
