require 'tools/HapPy'
module HTS
  CMD.tool :bcftools, nil, "bcftools -v" do
    CMD.cmd('conda install bcftools -c bioconda')
  end
  CMD.tool :bedtools, nil, "bedtools" do
    CMD.cmd('conda install bedtools -c bioconda')
  end

  Rbbt.claim Rbbt.software.opt.RTG, :install do
    {:src => 'https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.11/rtg-tools-3.11-nojre.zip'}
  end

  CMD.tool 'rtg', Rbbt.software.opt.RTG, 'rtg'

  input :truth_vcf, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
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

  input :input_vcf, :file, "VCF to get Transitions/Tranversions ratio"
  task :get_TSTV => :text do |input_vcf|
    CMD.cmd_log(:bcftools, "stats #{input_vcf} | grep -v '#' | grep TSTV | awk '{ print $5}'")
  end

  input :coverage_file, :file, "Coverage file in bed format"
  input :mutations_file, :file,"Mutations file in bed format"
  task :low_coverage => :text do |coverage_file, mutations_file|
    CMD.cmd(:bedtools, "intersect -a #{coverage_file} -b #{mutations_file}")
  end


  input :truth_vcf, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :somatic, :boolean, "Do somatic instead of germline", false
  task :vcfeval => :tsv do |truth,input,reference,somatic|
    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    truth = File.expand_path(truth) if Misc.is_filename? truth
    input = File.expand_path(input) if Misc.is_filename? input

    Misc.in_dir files_dir do
      input_sorted = File.join('.', "input.vcf")
      truth_sorted = File.join('.', "truth.vcf")
      sdf = File.join('.', "sdf")

      truth_sorted_orig = File.join('.', "truth.orig.vcf")
      truth_sorted_tmp = File.join('.', "truth.tmp.vcf")
      truth_sorted_tmp2 = File.join('.', "truth.tmp2.vcf")

      input_io = TSV.get_stream input
      CMD.cmd('bcftools', "sort |cut -f 1-10 > #{input_sorted}", :in => input_io)

      truth_io = TSV.get_stream truth
      Open.write(truth_sorted_orig, truth_io)

      CMD.cmd("echo '##fileformat=VCFv4.2' > #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##' #{truth_sorted_orig} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##contig' #{input_sorted} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##contig' #{input_sorted} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##FORMAT' #{input_sorted} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("echo '##FORMAT=<ID=,Number=R,Type=Integer,Description=>' >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("echo '##FORMAT=<ID=GT,Number=R,Type=String,Description=>' >> #{truth_sorted_tmp}", :nofail => true)

      ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the
      CMD.cmd("grep '#CHR' #{truth_sorted_orig} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep -v '#' #{truth_sorted_orig} |grep -v _alt| grep -v _random >> #{truth_sorted_tmp}", :nofail => true)

      Open.open(truth_sorted_tmp) do |io|
        Open.open(truth_sorted_tmp2, :mode => 'w') do |file|
          TSV.traverse io, :type => :array do |line|
            value = case line
                    when /^##/
                      line
                    when /^#/
                      line + "\t" + "FORMAT" + "\t"  + clean_name
                    else
                      next if line.split("\t")[4].split(",").include? line.split("\t")[3]
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

      if ! Open.read(input_sorted).include? "FORMAT=<ID=GT"
        TmpFile.with_file do |tmpfile|
          Path.setup(tmpfile)
          TSV.traverse input_sorted, :type => :array, :into => tmpfile do |line|
            if line =~ /^#CHR/
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n" + line
            elsif line =~ /^#/
              line
            else
              parts = line.split("\t")
              parts[8] += ":GT"
              parts[9] += ":0/1"
              parts * "\t"
            end
          end.join
          Open.mv tmpfile, input_sorted
        end
      end

      CMD.cmd('bcftools', "sort '#{truth_sorted_tmp2}' > #{truth_sorted}")

      CMD.cmd('bgzip', "#{truth_sorted}")
      CMD.cmd('bgzip', "#{input_sorted}")
      CMD.cmd('tabix', "#{truth_sorted}.gz")
      CMD.cmd('tabix', "#{input_sorted}.gz")

      CMD.cmd_log('rtg', "format #{reference} -o '#{sdf}'")

      text = CMD.cmd('rtg', "vcfeval --all-records -t #{sdf} -o '#{file('output')}' -b '#{truth_sorted}.gz' -c '#{input_sorted}.gz'").read.split("\n").reject{|l| l.include?("---") || l.include?("Selected")}
      text = text.collect{|line| line.gsub(/^  */,'')}
      TSV.open(StringIO.new(text * "\n"), :header_hash => '', :sep => /\s+/, :type => :list)
    end

  end


end
