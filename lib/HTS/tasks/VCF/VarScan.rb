require 'tools/VarScan'
require 'tools/fpfilter'
require 'pry'
module HTS

  input :normal, :file, "Normal BAM", nil, :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  dep :GC_windows, :jobname => :reference
  input :normal_purity, :float, "Normal sample purity", 1 
  input :tumor_purity, :float, "Tumor sample purity", 1 
  extension 'vcf'
  task :varscan_somatic_alt => :text do |normal,tumor,reference,normal_purity,tumor_purity|
    Open.mkdir files_dir
    output = file('output')
    pileup = file('pileup')

    reference = reference_file reference
    reference = Samtools.prepare_FASTA reference
    CMD.cmd("samtools mpileup -f '#{reference}' -Q 20 '#{normal}' '#{tumor}' > '#{pileup}'")
    io = Misc.in_dir output do
      monitor_cmd_genome ["java -jar #{Rbbt.software.opt.jars["VarScan.jar"].produce.find} somatic '#{pileup}' '#{clean_name}' --mpileup 1 --normal-purity #{normal_purity} --tumor-purity #{tumor_purity} --output-vcf '1' "], output[clean_name + '.snp.vcf']
    end

    ConcurrentStream.setup(io) do
      FileUtils.rm pileup
    end

    io
  end

  input :normal, :file, "Normal BAM", nil, :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  dep :pileup, :bam => :placeholder do |jobname,options|
    deps = []

    deps << {:inputs => options.merge(:bam => options[:normal]), :jobname => jobname + '.normal'}
    deps << {:inputs => options.merge(:bam => options[:tumor]), :jobname => jobname + '.tumor'}

    deps
  end
  dep :GC_windows, :jobname => :reference
  input :normal_purity, :float, "Normal sample purity", 1 
  input :tumor_purity, :float, "Tumor sample purity", 1 
  extension 'vcf'
  task :varscan_somatic => :text do |normal,tumor,reference,normal_purity,tumor_purity|
    normal = dependencies.select{|dep| dep.clean_name.include? '.normal'}.first
    tumor = dependencies.select{|dep| dep.clean_name.include? '.tumor'}.first

    output = file('output')
    Misc.in_dir output do
      io_normal = CMD.cmd("zcat |sed 's/\t/#/;s/\t/#/' ", :pipe => true, :in => TSV.get_stream(normal, :noz => true))
      io_tumor = CMD.cmd("zcat |sed 's/\t/#/;s/\t/#/' ", :pipe => true, :in => TSV.get_stream(tumor, :noz => true))

      pipe = Misc.paste_streams([io_normal, io_tumor]) do |a,b|
        Misc.genomic_location_cmp_strict(a, b, '#')
      end

      monitor_cmd_genome ["sed 's/#/\t/;s/#/\t/' | grep -v '[[:space:]][[:space:]]' | java -jar #{Rbbt.software.opt.jars["VarScan.jar"].produce.find} somatic --mpileup '#{clean_name}' --normal-purity #{normal_purity} --tumor-purity #{tumor_purity} --output-vcf '1' - ", :in => pipe], output[clean_name + '.snp.vcf']
    end
  end
  
  dep :varscan_somatic, :compute => :produce
  extension :vcf
  task :varscan_classify => :text do 
    output = file('output')
    Misc.in_dir output do
      Open.ln_s step(:varscan_somatic).path, "#{clean_name}.vcf"
      CMD.cmd("java -jar #{Rbbt.software.opt.jars["VarScan.jar"].produce.find} processSomatic #{clean_name}.vcf")
    end
    Open.cp output.glob("*.Somatic.vcf").first, self.path
    nil
  end

  dep :varscan_classify
  extension :vcf
  task :varscan_fpfiltered => :text do |jobname, options|
    orig_reference = reference_file recursive_inputs[:reference]
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    CMD.get_tool "bam-readcount"

    args = {}
    tumor = recursive_inputs[:tumor]
    tumor = tumor.path if Step === tumor

    args["bam-file"] = tumor
    args["vcf-file"] = step(:varscan_classify).path
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["sample"] = "TUMOR"
    FPFilter.filter(args.to_hash)

    TSV.traverse self.tmp_path, :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/

      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/
      #ssc = line.split(";").select{|d| d =~ /^SSC=/}.first.split.first.split("=").last.to_f
      #next unless ssc > 10

      line
    end
  end

  input :normal, :file, "Normal BAM", nil, :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  dep :pileup, :bam => :placeholder do |jobname,options|
    deps = []

    deps << {:inputs => options.merge(:bam => options[:normal]), :jobname => jobname + '.normal'}
    deps << {:inputs => options.merge(:bam => options[:tumor]), :jobname => jobname + '.tumor'}

    deps
  end
  dep :GC_windows, :jobname => :reference
  input :normal_purity, :float, "Normal sample purity", 1 
  input :tumor_purity, :float, "Tumor sample purity", 1 
  task :varscan_copynumber => :text do |normal,tumor,reference,normal_purity,tumor_purity|
    Open.mkdir files_dir
    normal = dependencies.select{|dep| dep.clean_name.include? '.normal'}.first
    tumor = dependencies.select{|dep| dep.clean_name.include? '.tumor'}.first
    output = file('output')
    Misc.in_dir output do
      io_normal = CMD.cmd("gunzip '#{normal.join.path}' -c | sed 's/\t/#/;s/\t/#/' ", :pipe => true)
      io_tumor = CMD.cmd("gunzip '#{tumor.join.path}' -c | sed 's/\t/#/;s/\t/#/' ", :pipe => true)

      pipe = Misc.paste_streams([io_normal, io_tumor]) do |a,b|
        Misc.genomic_location_cmp_strict(a, b, '#')
      end

      self.monitor_cmd_genome ["sed 's/#/\t/;s/#/\t/' | grep -v '[[:space:]][[:space:]]'| java -jar #{Rbbt.software.opt.jars["VarScan.jar"].produce.find} copynumber --mpileup '#{clean_name}' --normal-purity #{normal_purity} --tumor-purity #{tumor_purity} --output-vcf '1' - ", {:in => pipe}], clean_name + '.copynumber'
    end
  end
end
