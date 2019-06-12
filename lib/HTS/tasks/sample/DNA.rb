module Sample
  input :by_group, :boolean, "Separate files by read-group if RevertSam is required", false
  extension :bam
  dep_task :BAM, HTS, :BAM do |sample,options|
    sample_files = Sample.sample_files sample
    raise "Sample #{ sample } not found" if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:FASTQ]
      if Array === fastq_files.first && fastq_files.first.length > 1
        options = options.merge({:fastq1_files => fastq_files.first, :fastq2_files => (fastq_files.last || [])})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      end
    elsif uBAM_files = sample_files[:uBAM]
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({"HTS#uBAM" => uBAM_files})
        {:inputs => options, :jobname => sample}
      end
    elsif bam_files = sample_files[:BAM]
      options = options.merge({"HTS#BAM" => [bam_files].flatten.first})
      {:inputs => options, :jobname => sample}
    elsif orig_bam_files = sample_files["orig.BAM"]
      if options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign_by_group, :inputs => options, :jobname => sample}
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign, :inputs => options, :jobname => sample}
      end
    end
  end

  extension :bam
  dep_task :BAM_normal, Sample, :BAM do |sample,options|
    nsample = nil
    sample_files = nil
    sample_study = Sample.sample_study(sample)
    [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample
      
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files
    end

    {:inputs => options, :jobname => nsample} if sample_files
  end

  #{{{ EXPORT HTS METHODS

  #{{{ genomic_mutations and variant callers
  
  CALLERS = %w(strelka varscan mutect2 muse somatic_sniper delly svABA)
  {
    :strelka => :strelka,
    :varscan => :varscan_somatic,
    :mutect2 => [:mutect2, true],
    :muse => :muse,
    :somatic_sniper => :somatic_sniper,
    :delly => :delly,
    :svABA => :svABA,
    :sequenza_purity => :sequenza_purity,
    :sequenza_ploidy => :sequenza_ploidy,
  }.each do |task,otask|
    if Array === otask
      otask, allow_tumor_only = otask
    else
      allow_tumor_only = false
    end

    if allow_tumor_only
      dep :BAM, :compute => :bootstrap
      dep :BAM_normal, :compute => :bootstrap do |sample,options|
        nsample = nil
        sample_files = nil
        sample_study = Sample.sample_study(sample)
        [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
          break if sample_files
        end

        {:inputs => options, :jobname => sample} if sample_files
      end
    else
      dep :BAM, :compute => :bootstrap
      dep :BAM_normal, :compute => :bootstrap do |sample,options|
        nsample = nil
        sample_files = nil
        sample_study = Sample.sample_study(sample)
        [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
          break if sample_files
        end

        raise ParameterException, "No normal sample found" if sample_files.nil?

        {:inputs => options, :jobname => sample} 
      end
    end
    extension :vcf if CALLERS.include?(task.to_s)
    dep_task task, HTS, otask, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
      options = add_sample_options jobname, options

      if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.empty?
        options[:normal] = nil
      end
      {:inputs => options}
    end
  end

  dep :BAM, :compute => [:produce, :canfail]
  dep :BAM_normal, :compute => [:produce, :canfail]
  dep_task :haplotype, HTS, :haplotype, :BAM => :BAM_normal do |jobname, options, dependencies|
    options = add_sample_options jobname, options
    options[:BAM] = :BAM if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname}
  end

  dep :strelka, :compute => :bootstrap do |jobname, options|
    %w(strelka varscan mutect2 somatic_sniper muse).collect do |var_caller|
      {:task => var_caller, :inputs => options, :jobname => jobname}
    end
  end
  input :only_pass, :boolean, "Consider only PASS variants", true
  task :caller_comparison => :tsv do |only_pass|
    caller_variants = {}
    dependencies.each do |dep|
      var_caller = dep.task_name

      variants = TSV.traverse dep, :into => [], :type => :list do |chr, values|
        pos, id, ref, alt, qual, filter, *rest = values

        next if only_pass and ! (filter.split(";").include?("PASS")  || filter == ".")

        [chr, pos, alt] * ":"
      end
      caller_variants[var_caller] = variants
    end

    var_callers = caller_variants.keys.sort
    tsv = TSV.setup({}, :key_field => "Caller", :fields => ["Unique"] + var_callers, :type => :list, :cast => :to_i)
    var_callers.each do |c1|
      m1 = caller_variants[c1]
      matches = var_callers.collect do |c2|
        next if c1 == c2
        m2 = caller_variants[c2]
        m1 & m2
      end
      unique = m1 - matches.flatten
      tsv[c1] = [unique.length] + matches.collect{|m| m.nil? ? m1.length : m.length}
    end

    tsv
  end

  input :caller, :select, "Caller to use", :mutect2, :select_options => CALLERS
  dep Sample, :mutect2 do |sample,options|
    sample_files = Sample.sample_files sample
    if sample_files.include? "VCF"
      nil
    else
      {:task => :HTS_genomic_mutations, :inputs => options, :jobname => sample}
      vcaller = options[:caller]
      options[:vcf_file] = vcaller.to_sym
      {:task => vcaller, :jobname => sample, :inputs => options}
    end
  end
  extension "vcf"
  task :vcf_file => :text do 
    sample_files = Sample.sample_files self.clean_name
    if sample_files.include? "VCF"
      TSV.get_stream sample_files["VCF"].first
    else
      TSV.get_stream dependencies.first.path
    end
  end

  dep :vcf_file 
  dep_task :genomic_mutations, Sequence, :genomic_mutations, :vcf_file => :vcf_file 

  #{{{

  dep :BAM, :compute => :bootstrap
  dep :BAM_normal, :compute => :bootstrap do |sample,options|
    nsample = nil
    sample_files = nil
    sample_study = Sample.sample_study(sample)
    [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files
    end

    {:inputs => options, :jobname => sample} if sample_files
  end
  dep :genomic_mutations
  dep_task :mutation_BAM_img, HTS, :mutation_BAM_img, :tumor => :BAM, :normal => :BAM_normal, :positions => :genomic_mutations do |jobname,options|
    if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.empty?
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep_task :BAM_coverage, HTS, :BAM_coverage, :bam => :BAM do |jobname,options|
    options = add_sample_options jobname, options
    {:inputs => options}
  end

  dep :BAM
  dep_task :BAM_qualimap, HTS, :BAM_qualimap, :bam => :BAM, :interval_list => :placeholder do |jobname,options|
    options = add_sample_options jobname, options
    {:inputs => options}
  end

  dep :BAM_normal
  dep_task :BAM_qualimap_normal, HTS, :BAM_qualimap, :bam => :BAM_normal, :interval_list => nil do |jobname,options|
    nsample = nil
    sample_files = nil
    sample_study = Sample.sample_study(sample)
    [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files
    end

    if nsample
      options = add_sample_options nsample, options 
      {:inputs => options}
    end
  end


end