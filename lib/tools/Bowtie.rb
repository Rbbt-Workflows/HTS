require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/workflow'

Workflow.require_workflow "DbSNP"
module Bowtie
  extend Resource
  self.subdir = 'share/databases/Bowtie'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  def self.prepare_FASTA(file, dir = nil)
    file = file.find if Path === file

    digest = Misc.digest(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".1.bt2") || Persist.newer?(linked + ".1.bt2", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd("'#{Rbbt.software.opt.Bowtie2.url.local.bin["bowtie2-build"].find}' -f '#{ linked }' '#{linked}'")
      end
    end

    linked
  end

  Rbbt.claim Rbbt.software.opt.Bowtie, :install, Rbbt.share.install.software.Bowtie.find
  Rbbt.claim Rbbt.software.opt.Bowtie2, :install, Rbbt.share.install.software.Bowtie2.find
end

if __FILE__ == $0
  Log.severity = 0
  iii Rbbt.software.opt.Bowtie2.produce(true)
end
