module NovoAlign
  extend Resource
  self.subdir = 'share/databases/NovoAlign'

  Rbbt.claim Rbbt.software.opt.NovoAlign, :install, Rbbt.share.install.software.NovoAlign.find
  CMD.tool :novoindex, Rbbt.software.opt.NovoAlign
  CMD.tool :novoalign, Rbbt.software.opt.NovoAlign

  def self.prepare_FASTA(file, dir = nil)
    file = file.find if Path === file

    digest = Misc.digest(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir
    
    linked = dir[basename].find
    if ! File.exists?(linked + ".nix") || Persist.newer?(linked + ".nix", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd(:novoindex, "'#{ linked }.nix' '#{linked}'")
      end
    end

    linked
  end

end
