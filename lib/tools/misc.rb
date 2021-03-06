module HTS
  Rbbt.claim Rbbt.software.opt.BAM_readcount, :install, "https://github.com/genome/bam-readcount.git"
  CMD.tool "bam-readcount", Rbbt.software.opt.BAM_readcount

  def self.prepare_BED(file, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)

    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_indices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    CMD.cmd_log("bgzip -c #{file} > '#{ linked + '.gz' }'") unless File.exists?(linked + '.gz')
    if ! File.exists?(linked + ".gz.tbi") || Persist.newer?(linked + ".gz.tbi", file)
      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        CMD.cmd_log("tabix '#{ linked }.gz'")
      end
    end

    linked + '.gz'
  end


  def self.vcf_clean_IUPAC_alleles(input, output)
    iupac = Misc::IUPAC2BASE.select{|k,v| v.length > 1}.collect{|k,v| k}

    Open.write(output) do |foutput|
      TSV.traverse input, :into => foutput, :type => :array do |line|
        next if iupac.include?(line.split("\t")[3])
        line
      end
      foutput.join
    end
    nil
  end

  def self.combine_caller_vcfs(list)
    preambles = {}
    list.each do |name,file|
      preambles[name] = Open.read(file).split("\n").select{|line| line =~ /^#/ }
    end
    preamble = []

    fields = nil

    preambles.each do |name,lines|
      lines.each do |line|
        if line =~ /^#CHR/
          fields ||= line.split("\t")
          next
        end
        line = line.sub("FILTER=<ID=", "FILTER=<ID=#{name}--")
        line = line.sub("FORMAT=<ID=", "FORMAT=<ID=#{name}--")
        line = line.sub("INFO=<ID=", "INFO=<ID=#{name}--")
        preamble << line unless preamble.include?(line) || line =~ /(tumor|normal)_sample/
      end
    end


    list.keys.each do |name|
      preamble << "##FILTER=<ID=#{name}--PASS,Description=\"Passes #{name} filter\">" unless preamble.select{|l| l.include?("FILTER") && l.include?(name + "--PASS")}.any?
    end

    preamble << '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">'
    preamble << '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    preamble << '##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele fractions of alternate alleles in the tumor">'
    preamble << '##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of A alleles used in tiers 1,2">'
    preamble << '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles or variant allele">'
    preamble << '##FORMAT=<ID=BCOUNT,Number=4,Type=Integer,Description="Occurrence count for each base at this site (A,C,G,T)">'

    variants = {}
    list.each do |name,file|

      tumor_sample = HTS.guess_vcf_tumor_sample(file)
      file_fields = TSV.parse_header(file).fields
      sample1, sample2 = file_fields.values_at -2, -1
      swap_samples = sample1 == tumor_sample
      TSV.traverse file, :type => :array do |line|
        next if line =~ /^#/
        parts = line.split("\t")
        chr, pos, rsid, ref, alt, qual, filter, info, format, normal_sample, tumor_sample, *rest = parts
        normal_sample, tumor_sample = tumor_sample, normal_sample if swap_samples
        mutation = [chr, pos, ref, alt] * ":"
        filter = filter.split(";").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ";"
        info = info.split(";").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ";"
        format = format.split(":").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ":"
        variants[mutation] ||= []
        variants[mutation] << [filter, info, format, normal_sample, tumor_sample] + rest
      end
    end

    fields[-2] = "NORMAL"
    fields[-1] = "TUMOR"

    str = preamble * "\n" + "\n"
    str += fields * "\t" + "\n"

    variants.each do |mutation, lists|
      chr, pos, ref, alt = mutation.split(":")
      mfilter = []
      minfo = []
      mformat = []
      msample1 = []
      msample2 = []
      mrest = []
      lists.each do |filter, info, format, sample1, sample2,*rest|
        mfilter << filter
        minfo << info
        mformat += format.split(":")
        msample1 += sample1.split(":")
        msample2 += sample2.split(":")
        mrest << rest
      end

      new_format = []
      new_sample1 = []
      new_sample2 =[]

      common_format = %w(GT AF DP AU AD BCOUNT)
      common_format.each do |key|
        match = mformat.select{|k| k.split("--").last == key }.first
        next unless match
        kpos = mformat.index match
        new_format << match.partition("--").last
        new_sample1 << msample1[kpos]
        new_sample2 << msample2[kpos]
      end

      mformat = new_format + mformat
      msample1 = new_sample1 + msample1
      msample2 = new_sample2 + msample2

      str += ([chr, pos, '.', ref, alt, '.', mfilter * ";", minfo * ";", mformat * ":", msample1 * ":", msample2 * ":"] + Misc.zip_fields(mrest).collect{|l| l * ":"}) * "\t" 
      str += "\n"
    end

    str
  end

  def self.guess_vcf_tumor_sample(vcf)
    begin
      CMD.cmd("grep 'tumor_sample=' '#{vcf}'").read.strip.split("=").last
    rescue
      tsv = TSV.open(vcf, :type => :list)
      entry = tsv.keys.first
      fields = tsv.fields
      sample1, sample2 = fields.values_at 8, 9
      if fields.length == 9
        Log.warn "Could not find tumor_sample field in #{Misc.fingerprint(vcf)}, but only one sample: #{fields.last}"
        fields.last
      elsif fields.include? "TUMOR"
        Log.warn "Could not find tumor_sample field in #{Misc.fingerprint(vcf)}, using TUMOR"
        "TUMOR"
      elsif tsv[entry] && (genotype = tsv[entry][sample1].split(":").select{|p| p == "0/1" || p == "0|1"}.first)
        Log.warn "Could not find tumor_sample field in #{Misc.fingerprint(vcf)}, but #{sample1} has genotype #{genotype}"
        sample1
      else
        Log.warn "Could not find tumor_sample field in #{Misc.fingerprint(vcf)}, using last field"
        fields.last
      end
    end
  end

  def self.add_vcf_sample_header(vcf, tumor_sample, normal_sample)
    TSV.traverse vcf, :type => :array, :into => :stream do |line|
      next line unless line =~ /^(?:#|CHR)/
        if line =~ /^#?CHR/
          "##tumor_sample=#{tumor_sample}" + "\n" +
            "##normal_sample=#{normal_sample}" + "\n" +
            line
      else
        line
      end
    end
  end

  def self.add_vcf_genotype(vcf)

    found = false
    TSV.traverse vcf, :type => :array, :into => :stream do |line|
      found = true if line =~ /^##FORMAT=<ID=GT,/
      next line if line =~ /^##/
      if  line =~ /^#CHR/
        if ! found 
          format_line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' << "\n"
          next format_line + line 
        else
          next line
        end
      end
      parts = line.split("\t")

      format = parts[8].split(":")

      if ! format.include? "GT"
        parts[8] += ":GT"
        (9..parts.length-1).each do |pos|
          donor = parts[pos].split(":")
          hash = Misc.zip2hash(format, donor)

          parts[pos] += ":1/0"
        end
      end

      parts * "\t"
    end
  end

  def self.genome_info(organism)
    Persist.persist('genome_info', :yaml, :organism => organism) do
      sizes = Organism.chromosome_sizes(organism)
      chromosomes = chromosome_sizes.keys.sort{|a,b| Misc.chr_cmp_strict(a,b)}
      chromosome_offset = {}
      chromosomes.each_with_index do |chr,i|
        offset = 0
        chromosomes[0..i-1].each{|c| offset += chromosome_sizes[c] } if i > 0
        chromosome_offset[chr] = offset
      end
      total = Misc.sum(chromosome_sizes.collect{|k,v| v}).to_i

      {:total => total,
       :sizes => sizes,
       :offsets => chromosome_offset}
    end
  end

  def self.chromosome_progress(chr, pos, organism)

    genome_info = genome_info(organism)

    pos = pos.to_i
    chr = chr.gsub('chr', '')

    total, chromosome_offset = genome_info.values_at :total, :offsets
    offset = chromosome_offset[chr]
    [pos + offset, genome_info[:total]]
  end

end
