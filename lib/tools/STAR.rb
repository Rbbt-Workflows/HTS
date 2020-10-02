require 'rbbt-util'
require 'rbbt/resource'

module STAR
  extend Resource
  self.subdir = 'share/databases/STAR'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.STAR, :install, Rbbt.share.install.software.STAR.find

  CMD.tool "STAR", Rbbt.software.opt.STAR

  def self.build_gft_index(organism, read_length = 100, cpus = nil)
    cpus ||= Rbbt::Config.get("cpus", :STAR_build, :STAR, :default => 1)

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
                  'hg38_noalt'
                else 
                  Organism.hg_build(organism)
                end

    reference = HTS.helpers[:reference_file].call(reference) 

    reference = GATK.prepare_FASTA reference
    reference = Samtools.prepare_FASTA reference
    reference = HTS.uncompress_FASTA reference


    dir = Rbbt.var.STAR_indices[digest][read_length.to_s].find if dir.nil?
    Path.setup(dir) unless Path === dir


    linked = dir[basename].find
    if ! File.exists?(dir.Genome) || Persist.newer?(dir.Genome, file)
      Open.ln_s file, linked

      TmpFile.with_file do |unzipped|
        if Open.gzip?(linked)
          CMD.cmd("gunzip -c '#{ linked }' > #{unzipped}")
        else
          Open.ln_s linked, unzipped
        end

        Misc.in_dir dir do
          FileUtils.ln_s file, linked unless File.exists?(linked)
          CMD.cmd_log("STAR", "--runThreadN #{cpus} --runMode genomeGenerate --genomeDir #{dir} --genomeFastaFiles #{reference} --sjdbGTFfile #{unzipped} --sjdbOverhang #{read_length.to_i - 1}")
        end
      end
    end

    dir
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.STAR.produce
end
