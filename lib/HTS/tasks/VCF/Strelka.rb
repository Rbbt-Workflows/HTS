require 'tools/Strelka'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension :vcf
  task :strelka => :text do |tumor,normal,reference|
    output = file('output')
    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)

    cpus = config :cpus, :strelka, :default => 3

    Strelka.runSomatic(tumor, normal, reference, output, cpus)
    
    Open.read(output.results.variants["somatic.snvs.vcf.gz"])
  end
end