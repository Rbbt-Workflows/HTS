require 'tools/HapPy'
module HTS

  input :truth_vcf_file, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", nil, :select_options => %w(b37 hg38)
  task :hap_py => :tsv do |truth,input,reference|
    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference.sub!(/\.gz$/,'')
    Misc.in_dir files_dir do
      vcf_sorted = File.join('.', "input.vcf")
      
      truth_sorted = File.join('.', "truth.vcf")
      truth_io = TSV.get_stream truth
      CMD.cmd('bcftools', "sort > #{truth_sorted}", :in => truth_io)

      input_sorted = File.join('.', "input.vcf")
      input_io = TSV.get_stream input
      CMD.cmd('bcftools', "sort > #{input_sorted}", :in => input_io)

      CMD.cmd_log('hap.py', "-o '#{self.clean_name}' -r '#{reference}' #{truth_sorted} #{input_sorted}")
    end

    tsv = TSV.open(file(self.clean_name + '.summary.csv'), :sep => ',', :header_hash => '')
    tsv
  end

end
