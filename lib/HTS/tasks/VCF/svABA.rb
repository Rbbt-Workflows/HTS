require 'tools/svABA'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :svABA_options, :string, "Additional options"
  task :svABA => :array do |tumor,normal,reference,svABA_options|
    output = file('output')
    orig_reference = reference

    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)
    exclude = Delly.exclude(orig_reference)

    Open.mkdir files_dir
    SvABA.call(tumor, normal, reference, svABA_options, output)
    
    Dir.glob(output + '**')
  end
end
