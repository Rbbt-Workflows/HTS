require 'rbbt-util'
require 'rbbt/resource'

module Samtools
  extend Resource
  self.subdir = 'share/databases/Samtools'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.HTSLib, :install do
    url="https://github.com/samtools/htslib.git"
    extra = "--disable-s3"
    {:git => url, :extra => extra}
  end
  Rbbt.claim Rbbt.software.opt.Samtools, :install do
    Rbbt.software.opt.HTSLib.produce
    url = "https://github.com/samtools/samtools.git"
    extra = "--with-htslib='#{Rbbt.software.opt.HTSLib.find}'"
    {:git => url, :extra => extra}
  end
  Rbbt.claim Rbbt.software.opt.bcftools, :install do
    Rbbt.software.opt.HTSLib.produce
    url = "https://github.com/samtools/bcftools.git"
    extra = "--with-htslib='#{Rbbt.software.opt.HTSLib.find}'"
    {:git => url, :extra => extra}
  end

  CMD.tool :samtools, Rbbt.software.opt.Samtools

  CMD.tool :bcftools, Rbbt.software.opt.bcftools

  #Samtools_CMD='samtools'
  
  #Samtools_CMD=Rbbt.software.opt.Samtools.produce.bin.samtools.find

  def self.samtools_cmd
    'samtools'
  end

  def self.run_log(command)
    CMD.cmd_log(:samtools, command)
  end

  def self.run(command)
    CMD.cmd(:samtools, command).read
  end

  def self.prepare_BAM(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.bam_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    CMD.get_tool :samtools
    linked = dir[basename].find

    Misc.lock linked do
      if ! (
        (File.exists?(linked + ".bai") && ! Persist.newer?(linked + ".bai", file)) ||
        (File.exists?(linked + ".crai") && ! Persist.newer?(linked + ".crai", file))
      )

        Misc.in_dir dir do
          Open.ln_s file, linked unless Open.exists?(linked)
          cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
          if cpus
            Samtools.run("index -@ #{cpus} '#{ linked }'")
          else
            Samtools.run("index '#{ linked }'")
          end
        end
      end
    end

    if File.exists?(linked + '.crai') && linked =~ /\.bam$/
      linked_cram = linked.sub(/\.bam$/,'.cram')
      Open.ln_s linked, linked_cram
      Open.ln_s linked + '.crai', linked_cram + '.crai'
      linked_cram
    else
      linked
    end
  end
  
  def self.prepare_FASTA(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    CMD.get_tool :samtools

    linked = dir[basename].find
    if ! File.exists?(linked + ".fai") || Persist.newer?(linked + ".fai", file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        Samtools.run("faidx '#{ linked }'")
      end
    end

    linked
  end
  
  def self.BAM_sort(bam_file, to_sam=false, queryname=true)
  	cpus = Rbbt::Config.get("cpus", :samtools_sort, :samtools, :sort, :default => nil)
    format_str=" -O BAM "
    if to_sam
      format_str=" -O SAM "
    end
    sort_order_str=""
    if queryname
      sort_order_str="-n"
    end
    if cpus
      Samtools.run("sort #{format_str} #{sort_order_str} -@ #{cpus} '#{bam_file}' -o #{bam_file}.sorted")
    else
      Samtools.run("sort #{format_str} #{sort_order_str} '#{bam_file}' -o #{bam_file}.sorted")
    end
  end

  def self.BAM_get_chr_reads(bam_file, chr)
    cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    CMD.cmd("samtools view --no-PG -@ #{cpus} -Sb -h '#{bam_file}' #{chr} > #{chr}.bam")
  end

  def self.header(bam_file)
    CMD.cmd("samtools view --no-PG -H #{bam_file}| grep -v \"^@PG\"").read
  end

  def self.reads_number(bam_file)
    CMD.cmd("samtools view --no-PG -c #{bam_file}").read
  end

  def self.BAM_start(bam_file)
    CMD.cmd("samtools view --no-PG '#{bam_file}'| head -n 1 | cut -f 3,4").read.strip.split("\t")
  end

  def self.merge(outbam, inbams)
    cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    CMD.cmd("samtools merge --no-PG -n -@ #{cpus} #{outbam} #{inbams} ")
  end

  def self.BAM_sample_name(bam_file)
    CMD.cmd(:samtools, "view --no-PG -H #{bam_file} | grep \"^@RG\"|grep \"SM:\" | awk '{for(i=1;i<=NF;i++)   if ( $i ~ /SM:.*/ ){ wln=$i; break } ; printf \"%s\\n\",wln}'| cut -d: -f2 | head -n 1").read.strip
  end

  def self.viewSam(inBAM, outSAM)
    cpus = Rbbt::Config.get("cpus", :samtools_index, :samtools, :index, :default => nil)
    CMD.cmd("samtools view --no-PG -@ #{cpus} #{inBAM} > #{outSAM}")
  end

  def self.reference_contigs(reference)
    Open.read(reference + '.fai').split("\n").collect{|line| line.split("\t").first}
  end

  def self.contig_sizes(reference)
    sizes = {}
    Open.read(reference + '.fai').split("\n").collect do |line| 
      name, size, *rest = line.split("\t")
      sizes[name] = size.to_i
    end
    sizes
  end

  def self.bam_contigs(bam)
    bam = bam.path if Step === bam
    bam = bam.find if Path === bam
    CMD.cmd("samtools view --no-PG -H '#{bam}' | grep ^@SQ | cut -f 2 |sed 's/SN://'").read.split("\n")
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HTSLib.produce
  iif Rbbt.software.opt.Samtools.produce 
end
