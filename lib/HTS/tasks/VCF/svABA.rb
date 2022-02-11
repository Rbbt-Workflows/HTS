require 'tools/svABA'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :svABA_options, :string, "Additional options"
  task :svABA => :array do |tumor,normal,reference,svABA_options|
    output = file('output')
    orig_reference = reference

    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)

    Open.mkdir files_dir
    cpus = config :cpus, :svaba, :svABA, :default => 1
    SvABA.call(tumor, normal, reference, svABA_options, output, cpus)
    
    Dir.glob(output + '/**')
  end

  dep :svABA
  extension :vcf
  task :svABA_indels => :text do
    file = step(:svABA).join.file('output/svABA.svaba.somatic.indel.vcf')
    tumor_bam = step(:svABA).inputs[:tumor]
    normal_bam = step(:svABA).inputs[:normal]
    SvABA.fix_vcf_sample_names(file, tumor_bam, normal_bam)
  end
end
