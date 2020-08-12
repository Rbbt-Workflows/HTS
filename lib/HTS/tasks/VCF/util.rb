require 'tools/HapPy'
module HTS
  CMD.tool :bcftools, nil, "bcftools -v" do
    CMD.cmd('conda install bcftools -c bioconda')
  end
  CMD.tool :bedtools, nil, "bedtools" do
    CMD.cmd('conda install bedtools -c bioconda')
  end



  input :truth_vcf_file, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
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

  input :input_vcf, :file, "VCF to get Transitions/Tranversions ratio"
  task :get_TSTV => :text do |input_vcf|
    CMD.cmd_log(:bcftools, "stats #{input_vcf} | grep -v '#' | grep TSTV | awk '{ print $5}'")
  end

  input :coverage_file, :file, "Coverage file in bed format"
  input :mutations_file, :file,"Mutations file in bed format"
  task :low_coverage => :text do |coverage_file, mutations_file|
    CMD.cmd(:bedtools, "intersect -a #{coverage_file} -b #{mutations_file}")
  end


end
