module HTS

  Rbbt.claim Rbbt.software.opt.RTG, :install do
    {:src => 'https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.11/rtg-tools-3.11-nojre.zip'}
  end

  CMD.tool 'rtg', Rbbt.software.opt.RTG, 'rtg'

  input :truth_vcf, :file, "Truth VCF"
  input :input_vcf, :file, "VCF to compare"
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :truth_sample, :string, "Tumor sample name in truth VCF"
  input :input_sample, :string, "Tumor sample name in input VCF"
  input :bed_regions, :file, "Regions to evaluate"
  task :vcfeval => :tsv do |truth,input,reference,truth_sample,input_sample,bed_regions|
    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    truth = File.expand_path(truth) if Misc.is_filename? truth
    input = File.expand_path(input) if Misc.is_filename? input

    Misc.in_dir files_dir do
      input_sorted = File.join('.', "input.vcf")
      truth_sorted = File.join('.', "truth.vcf")
      sdf = File.join('.', "sdf")

      input_io = TSV.get_stream input
      #Misc.sensiblewrite(input_sorted, HTS.add_vcf_genotype(CMD.cmd('bcftools', "sort", :in => input_io)))
      Misc.sensiblewrite(input_sorted, HTS.add_vcf_genotype(CMD.cmd("sort -k 1,1 -k2,2n", :in => input_io)))

      truth_io = TSV.get_stream truth
      #Misc.sensiblewrite(truth_sorted, HTS.add_vcf_genotype(CMD.cmd('bcftools', "sort", :in => truth_io)))
      Misc.sensiblewrite(truth_sorted, HTS.add_vcf_genotype(CMD.cmd('sort -k 1,1 -k2,2n', :in => truth_io)))

      truth_sample = HTS.guess_vcf_tumor_sample(truth_sorted)
      input_sample = HTS.guess_vcf_tumor_sample(input_sorted)

      begin
        CMD.cmd("grep -v '#' #{truth_sorted}")
      rescue
        raise ParameterException, "No mutations in truth set"
      end

      begin
        CMD.cmd("grep -v '#' #{input_sorted}")
      rescue
        raise ParameterException, "No mutations in input set"
      end


      CMD.cmd(:bgzip, truth_sorted)
      CMD.cmd(:bgzip, input_sorted)
      CMD.cmd(:tabix, "#{truth_sorted}.gz")
      CMD.cmd(:tabix, "#{input_sorted}.gz")

      sdf = Persist.persist("RTG SDF #{orig_reference}", :path, :other => {:reference => reference}) do |filename|
        CMD.cmd_log('rtg', "format #{reference} -o '#{filename}'")
      end

      Open.ln_s sdf, './sdf'

      if bed_regions && File.exist?(bed_regions)
        text = CMD.cmd('rtg', "vcfeval --decompose --squash-ploidy --no-roc --all-records -t '#{sdf}' --bed-regions #{bed_regions} --sample '#{truth_sample},#{input_sample}' -o '#{file('output')}' -b '#{truth_sorted}.gz' -c '#{input_sorted}.gz'").read.split("\n").reject{|l| l.include?("---") || l.include?("Selected")}
      else
        text = CMD.cmd('rtg', "vcfeval --decompose --squash-ploidy --no-roc --all-records -t '#{sdf}' --sample '#{truth_sample},#{input_sample}' -o '#{file('output')}' -b '#{truth_sorted}.gz' -c '#{input_sorted}.gz'").read.split("\n").reject{|l| l.include?("---") || l.include?("Selected")}
      end
      text = text.collect{|line| line.gsub(/^  */,'')}
      TSV.open(StringIO.new(text * "\n"), :header_hash => '', :sep => /\s+/, :type => :list)
    end

  end
end
