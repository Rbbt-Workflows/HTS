require 'tools/Delly'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :vcf
  task :delly => :text do |tumor,normal,reference|
    output = file('output')
    orig_reference = reference

    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)
    exclude = Delly.exclude(orig_reference)

    Open.mkdir files_dir
    Delly.call(tumor, normal, reference, exclude, output)
    
    Open.read(output + '.vcf')
  end
end
