require 'rbbt-util'
require 'rbbt/resource'

module HISAT
  extend Resource
  self.subdir = 'share/databases/HISAT'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.HISAT, :install, Rbbt.share.install.software.HISAT.find


  def self.build_gft_index(organism, cpus = nil)
    cpus ||= Rbbt::Config.get("cpus", :hisat_build, :hisat)

    file = HTS.gtf_file(organism)

    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    reference = case Organism.hg_build(organism)
                when 'hg19'
                  'b37'
                when 'hg38'
                  'hg38'
                else
                  Organism.hg_build(organism)
                end

    reference = HTS.helpers[:reference_file].call(reference) 

    reference = GATK.prepare_FASTA reference
    reference = Samtools.prepare_FASTA reference
    reference = HTS.uncompress_FASTA reference


    dir = Rbbt.var.HISAT_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir


    linked = dir[basename].find
    if ! File.exists?(linked + ".idx.1.ht2") || Persist.newer?(linked + '.idx.1.ht2', file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd_log("bash -c 'extract_splice_sites.py #{Open.gzip?(linked) ? "<(gunzip -c '#{linked}')" : "'#{ linked }'"} > #{linked}.ss'")
        CMD.cmd_log("bash -c 'extract_exons.py #{Open.gzip?(linked) ? "<(gunzip -c '#{linked}')" : "'#{ linked }'"} > #{linked}.exon'")

        if Organism.hg_build(organism) == 'hg38'
          CMD.cmd_log(%(sed -i 's/^\\([0-9A-Z]\\)/chr\\1/' #{linked}.ss))
          CMD.cmd_log(%(sed -i 's/^\\([0-9A-Z]\\)/chr\\1/' #{linked}.exon))
        end

        CMD.cmd_log("hisat2-build -p #{cpus || 1} --ss #{linked}.ss --exon #{linked}.exon #{reference} #{linked}.idx ")
      end
    end

    linked + '.idx'
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HISAT.produce
end
