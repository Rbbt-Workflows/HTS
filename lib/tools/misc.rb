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
          fields ||= line
          next
        end
        line = line.sub("FILTER=<ID=", "FILTER=<ID=#{name}--")
        line = line.sub("FORMAT=<ID=", "FORMAT=<ID=#{name}--")
        line = line.sub("INFO=<ID=", "INFO=<ID=#{name}--")
        preamble << line unless preamble.include?(line)
      end
    end

    list.keys.each do |name|
      preamble << "##FILTER=<ID=#{name}--PASS,Description=\"Passes #{name} filter\">" unless preamble.select{|l| l.include?("FILTER") && l.include?(name + "--PASS")}.any?
    end

    variants = {}
    list.each do |name,file|
      TSV.traverse file, :type => :array do |line|
        next if line =~ /^#/
        chr, pos, rsid, ref, alt, qual, filter, info, format, sample, *rest = line.split("\t")
        mutation = [chr, pos, ref, alt] * ":"
        filter = filter.split(";").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ";"
        info = info.split(";").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ";"
        format = format.split(":").reject{|f| f == '.'}.collect{|f| name + "--" + f} * ":"
        variants[mutation] ||= []
        variants[mutation] << [filter, info, format, sample] + rest
      end
    end

    str = preamble * "\n" + "\n"
    str += fields + "\n"

    variants.each do |mutation, lists|
      chr, pos, ref, alt = mutation.split(":")
      mfilter = []
      minfo = []
      mformat = []
      msample = []
      mrest = []
      lists.each do |filter, info, format, sample,*rest|
        mfilter << filter
        minfo << info
        mformat << format
        msample << sample
        mrest << rest
      end
      str += ([chr, pos, '.', ref, alt, '.', mfilter * ";", minfo * ";", mformat * ":", msample * ":"] + Misc.zip_fields(mrest).collect{|l| l * ":"}) * "\t" 
      str += "\n"
    end

    str
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

end
