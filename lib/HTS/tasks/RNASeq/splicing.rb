module HTS

  Rbbt.claim Rbbt.software.opt.gclib, :install, "https://github.com/gpertea/gclib.git"
  Rbbt.claim Rbbt.software.opt.gffcompare, :install do
    Rbbt.software.opt.gclib.produce
    {:git => "https://github.com/gpertea/gffcompare.git"}
  end

  dep :RNA_BAM
  extension :gtf
  input :rna_strandness, :select, "RNA Strandness", 'FR', :select_options => %w(FR)
  input :organism, :string, "Organism code"
  task :stringtie => :tsv do |rna_strandness,organism|
    organism ||= self.recursive_inputs[:organism] || dependencies.first.info[:organism]
    rna_strandness = self.recursive_inputs[:rna_strandness]

    gtf_file_gz = HTS.gtf_file(organism)

    if Open.gzip?(gtf_file_gz)
      gtf_file = file('gtf_file') 
      Open.write(gtf_file, Open.read(gtf_file_gz))
    else
      gtf_file = gtf_file_gz
    end

    strand_arg = rna_strandness == "FR" ? "--fr" : (rna_strandness.nil? ? "" : "--rf")

    cpus = config("cpus", :stringtie, :default => 1)
    CMD.cmd_log("stringtie -p #{cpus} #{strand_arg} -G #{gtf_file} -o #{self.tmp_path} #{step(:RNA_BAM).path}")
    nil
  end

  dep :stringtie
  input :organism, :string, "Organism code"
  task :splicing_variants => :string do |organism|
    organism ||= self.recursive_inputs[:organism] || dependencies.first.info[:organism]

    gtf_file_gz = HTS.gtf_file(organism)

    gtf_file = file('reference') 
    Open.write(gtf_file, Open.read(gtf_file_gz))

    stringtie = file('stringtie') 
    Open.write(stringtie, Open.read(step(:stringtie).path))

    output = file('output')
    Open.mkdir output

    Rbbt.software.opt.gffcompare.produce
    CMD.cmd_log("gffcompare -r '#{gtf_file}' '#{stringtie}' -o '#{output}/#{self.clean_name}'")
    "DONE"
  end
end
