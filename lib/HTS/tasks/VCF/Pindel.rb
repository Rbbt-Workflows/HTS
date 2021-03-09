require 'tools/Pindel'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :insert_size, :integer, "Insert size of pair ends"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :pindel_pre => :array do |tumor,normal,insert_size,reference|
    output = file('output')
    orig_reference = reference

    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)

    Open.mkdir files_dir
    cpus = config :cpus, :pindel, :Pindel, :default => 1
    Pindel.call(tumor, normal, reference, insert_size, output, cpus)
    
    Dir.glob(File.join(output,'*'))
  end

  dep :pindel_pre
  extension :vcf
  input :min_supporting_reads, :integer, "Minimun number of reads supporting indel", 10
  task :pindel_indels => :text do |min_supporting_reads|
    pindel_pre = step(:pindel_pre)
    reference = pindel_pre.inputs[:reference]
    orig_reference = reference

    reference = reference_file reference
    reference = Samtools.prepare_FASTA(reference)

    reference = reference.remove_extension('.gz')

    reference_date = '20210101'
    si_output = pindel_pre.file('output/pindel_SI').read
    d_output = pindel_pre.file('output/pindel_D').read
    TmpFile.with_file(si_output + "\n" + d_output) do |output|
      CMD.cmd_log(:pindel2vcf, "-G --min_supporting_reads #{ min_supporting_reads } -p #{output} -r #{ reference } -R #{ orig_reference } -d #{reference_date} -v #{self.tmp_path}")
    end
    nil
  end
end
