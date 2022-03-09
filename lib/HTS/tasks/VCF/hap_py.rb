module HTS
  input :truth_vcf, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :somatic, :boolean, "Do somatic instead of germline", false
  task :hap_py => :tsv do |truth,input,reference,somatic|
    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    truth = File.expand_path(truth) if Misc.is_filename? truth
    input = File.expand_path(input) if Misc.is_filename? input

    Misc.in_dir files_dir do
      input_sorted = File.join('.', "input.vcf")
      truth_sorted = File.join('.', "truth.vcf")
      truth_sorted_orig = File.join('.', "truth.orig.vcf")
      truth_sorted_tmp = File.join('.', "truth.tmp.vcf")
      truth_sorted_tmp2 = File.join('.', "truth.tmp2.vcf")

      input_io = TSV.get_stream input
      CMD.cmd('bcftools', "sort > #{input_sorted}", :in => input_io)

      truth_io = TSV.get_stream truth
      Open.write(truth_sorted_orig, truth_io)

      CMD.cmd("grep '##' #{truth_sorted_orig} > #{truth_sorted_tmp}")
      CMD.cmd("grep '##contig' #{input_sorted} >> #{truth_sorted_tmp}")
      CMD.cmd("grep '##FORMAT' #{input_sorted} >> #{truth_sorted_tmp}")
      CMD.cmd("echo '##FORMAT=<ID=,Number=R,Type=Integer,Description=>' >> #{truth_sorted_tmp}")

      ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the
      CMD.cmd("grep '#CHR' #{truth_sorted_orig} >> #{truth_sorted_tmp}")
      CMD.cmd("grep -v '#' #{truth_sorted_orig} |grep -v _alt| grep -v _random >> #{truth_sorted_tmp}")

      Open.open(truth_sorted_tmp) do |io|
        Open.open(truth_sorted_tmp2, :mode => 'w') do |file|
          TSV.traverse io, :type => :array do |line|
            value = case line
                    when /^##/
                      line
                    when /^#/
                      line + "\t" + "FORMAT" + "\t"  + clean_name
                    else
                      line = (line.split("\t")[0..4] + ["", ".",""]) * "\t"
                      if line =~ /^chr/ || ! reference.include?("hg38")
                        line + "\t" + "GT" + "\t"  + "0/1"
                      else
                        "chr" + line + "\t" + "GT" + "\t"  + "0/1"
                      end
                    end
            file.puts value
          end
        end
      end

      CMD.cmd('bcftools', "sort '#{truth_sorted_tmp2}' > #{truth_sorted}")

      if somatic
        CMD.cmd_log('som.py', " --no-roc -o '#{self.clean_name}' --somatic -r '#{reference}' #{truth_sorted} #{input_sorted}")
      else
        CMD.cmd_log('hap.py', " --no-roc -o '#{self.clean_name}' -r '#{reference}' #{truth_sorted} #{input_sorted}")
      end
    end

    tsv = TSV.open(file(self.clean_name + '.summary.csv'), :sep => ',', :header_hash => '')
    tsv
  end

  input :vcf_file, :file, "Imput VCF file"
  extension :vcf
  task :pre_py => :text do |vcf_file|
    CMD.cmd_log('hap.py', " --no-roc -o '#{self.clean_name}' -r '#{reference}' #{truth_sorted} #{input_sorted}")
  end
end
