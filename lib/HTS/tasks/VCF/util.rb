require 'tools/HapPy'
module HTS
  CMD.tool :bcftools, nil, "bcftools -v" do
    CMD.cmd('conda install bcftools -c bioconda')
  end
  CMD.tool :bedtools, nil, "bedtools" do
    CMD.cmd('conda install bedtools -c bioconda')
  end

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

  input :mutations, :array, "Array of genomic mutations"
  input :bam, :file, "BAM file", nil, :nofile => true
  input :add_chr, :boolean, "Add chr prefix to BED contigs (auto-detects if not specified)"
  task :mutation_read_counts => :tsv do |mutations,bam,add_chr|

    bam_contigs =  Samtools.bam_contigs(bam)

    add_chr = bam_contigs.first.include? 'chr' if add_chr.nil?
    add_txt = add_chr ? 'add' : 'remove'

    TmpFile.with_file(Misc.genomic_mutations_to_BED(mutations, add_txt, bam_contigs)) do |bed|
      TmpFile.with_file(CMD.cmd("samtools view -H '#{bam}' | grep -P \"@SQ\tSN:\" | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/'", :pipe => true)) do |genome|
        CMD.cmd_log(:bedtools, "coverage -counts -sorted -a #{bed} -b '#{bam}' -g #{genome} > #{self.tmp_path}")
      end
    end
    nil
  end

  input :mutations, :array, "Array of genomic mutations"
  input :bam, :file, "BAM file", nil, :nofile => true
  input :add_chr, :boolean, "Add chr prefix to BED contigs (auto-detects if not specified)"
  task :mutation_pileup => :tsv do |mutations,bam,add_chr|

    bam_contigs =  Samtools.bam_contigs(bam)

    add_chr = bam_contigs.first.include? 'chr' if add_chr.nil?
    add_txt = add_chr ? 'add' : 'remove'

    mutations = mutations.read.split("\n") if IO === mutations

    mutations = mutations.select{|m| %w(A C T G).include? m.split(":").last}

    bedfile = file('regions.bed')
    pileup = file('pileup')

    Open.write(bedfile, Misc.genomic_mutations_to_BED(mutations, add_txt, bam_contigs))

    CMD.cmd(:samtools, "mpileup -l '#{bedfile}' '#{bam}' -o '#{pileup}' -a")
    CMD.cmd("paste '#{pileup}' '#{bedfile}' > #{self.tmp_path}")

    nil
  end

  input :mutations, :array, "Array of genomic mutations"
  input :bam, :file, "BAM file", nil, :nofile => true
  input :add_chr, :boolean, "Add chr prefix to BED contigs (auto-detects if not specified)"
  task :mutation_support => :tsv do |mutations,bam,add_chr|

    bam_contigs =  Samtools.bam_contigs(bam)

    add_chr = bam_contigs.first.include? 'chr' if add_chr.nil?
    add_txt = add_chr ? 'add' : 'remove'

    pileup = file('pileup')
    mutations = mutations.select{|m| %w(A C T G).include? m.split(":").last}

    tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Coverage", "Variant Reads"], :type => :list, :cast => :to_i)

    mutations.each do |mutation|
      tsv[mutation] = [0, 0]
    end

    TmpFile.with_file(Misc.genomic_mutations_to_BED(mutations, add_txt, bam_contigs)) do |bed|
      io = CMD.cmd(:samtools, "mpileup -l #{bed} '#{bam}' ", :pipe => true)
      TSV.traverse io, :type => :array, :bar => self.progress_bar("Calculating mutation support") do |line|
        chr, pos, ref, depth, bases = line.split("\t")
        depth = depth.to_i
        counts = Misc.counts(bases.upcase.chars)

        %w(A C T G).each do |alt|
          m = [chr, pos, alt] * ":"
          if tsv[m]
            tsv[m] = [depth, counts[alt]]
          end
        end
      end
    end

    tsv
  end



end
