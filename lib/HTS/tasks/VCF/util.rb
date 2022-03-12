require 'tools/HapPy'
module HTS
  CMD.tool :bcftools, nil, "bcftools -v" do
    CMD.cmd('conda install bcftools -c bioconda')
  end
  CMD.tool :bedtools, nil, "bedtools" do
    CMD.cmd('conda install bedtools -c bioconda')
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

  input :vcf1, :file, "File 1", nil, :nofile => true
  input :vcf2, :file, "File 2", nil, :nofile => true
  extension :vcf
  task :join_vcfs => :text do |vcf1,vcf2|
    header_lines = []

    Open.open(vcf1) do |f|
      header_lines += CMD.cmd("grep '^##'", :in => f, :nofail => true).read.split("\n")
    end

    Open.open(vcf2) do |f|
      header_lines += CMD.cmd("grep '^##' | grep -v '##SAMPLE' ", :in => f, :nofail => true).read.split("\n")
    end


    header_lines.uniq!

    Misc.open_pipe do |f|
      f.puts header_lines * "\n"

      Open.open(vcf1) do |v|
        io = CMD.cmd("grep -v '^##' ", :in => v, :pipe => true, :nofail => true)
        Misc.consume_stream(io, false, f, false)
      end

      Open.open(vcf2) do |v|
        io = CMD.cmd("grep -v '^#' ", :in => v, :pipe => true, :nofail => true)
        Misc.consume_stream(io, false, f)
      end
    end
  end


end
