module Sample
  dep_task :kallisto, HTS, :kallisto do |sample,options|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:RNA_FASTQ]
      if Array === fastq_files.first && fastq_files.first.length > 1
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      else
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      end
    elsif uBAM_files = sample_files[:uBAM]
      raise "No support for uBAM"
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({"HTS#uBAM" => uBAM_files})
        {:inputs => options, :jobname => sample}
      end
    elsif bam_files = sample_files[:BAM]
      raise "No support for BAM"
      options = options.merge({"HTS#BAM" => [bam_files].flatten.first})
      {:inputs => options, :jobname => sample}
    elsif orig_bam_files = sample_files["orig.BAM"]
      raise "No support for orig.BAM"
      if options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign_by_group, :inputs => options, :jobname => sample}
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign, :inputs => options, :jobname => sample}
      end
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample}
    end
  end


  dep_task :salmon, HTS, :salmon do |sample,options|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:RNA_FASTQ]
      if Array === fastq_files.first && fastq_files.first.length > 1
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      else
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      end
    elsif uBAM_files = sample_files[:uBAM]
      raise "No support for uBAM"
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({"HTS#uBAM" => uBAM_files})
        {:inputs => options, :jobname => sample}
      end
    elsif bam_files = sample_files[:BAM]
      raise "No support for BAM"
      options = options.merge({"HTS#BAM" => [bam_files].flatten.first})
      {:inputs => options, :jobname => sample}
    elsif orig_bam_files = sample_files["orig.BAM"]
      raise "No support for orig.BAM"
      if options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign_by_group, :inputs => options, :jobname => sample}
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign, :inputs => options, :jobname => sample}
      end
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample}
    end
  end

  dep_task :RNA_BAM, HTS, :RNA_BAM, :fastq1 => :placeholder, :fastq2 => :placeholder do |sample,options|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:RNA_FASTQ]
      options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
      {:inputs => options, :jobname => sample}
    elsif uBAM_files = sample_files[:uBAM]
      raise "No support for uBAM"
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({"HTS#uBAM" => uBAM_files})
        {:inputs => options, :jobname => sample}
      end
    elsif bam_files = sample_files[:BAM]
      raise "No support for BAM"
      options = options.merge({"HTS#BAM" => [bam_files].flatten.first})
      {:inputs => options, :jobname => sample}
    elsif orig_bam_files = sample_files["orig.BAM"]
      raise "No support for orig.BAM"
      if options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign_by_group, :inputs => options, :jobname => sample}
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign, :inputs => options, :jobname => sample}
      end
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample} 
    end
  end

  dep :RNA_BAM
  dep_task :stringtie, HTS, :stringtie do |sample,options,dependencies|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    dep_RNA_BAM = dependencies.flatten.select{|dep| dep.task_name == :RNA_BAM }.first
    options = options.merge("HTS#RNA_BAM" => dep_RNA_BAM, :not_overriden => true)
    {:inputs => options, :jobname => sample}
  end

  dep :stringtie
  dep_task :htseq_count, HTS, :htseq_count do |sample,options,dependencies|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    dep_stringtie = dependencies.flatten.select{|dep| dep.task_name == :stringtie }.first
    options = options.merge("HTS#stringtie" => dep_stringtie, :not_overriden => true)
    {:inputs => options, :jobname => sample}
  end

  dep :stringtie
  dep_task :splicing_variants, HTS, :splicing_variants do |sample,options,dependencies|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    dep_stringtie = dependencies.flatten.select{|dep| dep.task_name == :stringtie }.first
    options = options.merge("HTS#stringtie" => dep_stringtie)
    {:inputs => options, :jobname => sample}
  end


  dep_task :RNA_BAM_normal, Sample, :RNA_BAM do |sample,options|
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
      dep :RNA_BAM, :compute => :bootstrap
      dep :RNA_BAM_normal, :compute => :bootstrap do |sample,options|
        nsample = nil
        sample_files = nil
        sample_study = Sample.sample_study(sample)
        [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
          break if sample_files
        end

        if sample_files && sample_files[:RNA_FASTQ]
          {:inputs => options, :jobname => sample} 
        elsif sample_files
          {:task => :BAM_normal, :inputs => options, :jobname => sample} if sample_files
        else
          nil
        end
      end
    else
      dep :RNA_BAM, :compute => :bootstrap
      dep :RNA_BAM_normal, :compute => :bootstrap do |sample,options|
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
    dep_task "RNA_" + task.to_s, HTS, otask, :normal => :RNA_BAM_normal, :tumor => :RNA_BAM do |jobname,options,dependencies|
      options = add_sample_options jobname, options

      if dependencies.flatten.select{|dep| dep.task_name == :RNA_BAM_normal}.empty?
        options[:normal] = nil
      end
      {:inputs => options}
    end
  end

  dep :RNA_mutect2
  dep_task :expanded_vcf_RNA, Sequence, :expanded_vcf, :vcf_file => :RNA_mutect2

  dep_task :FuSeq, HTS, :FuSeq_process do |sample,options|
    sample_files = Sample.sample_files sample
    raise "Sample #{ sample } not found" if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:RNA_FASTQ]
      if Array === fastq_files.first && fastq_files.first.length > 1
        raise "No support for multiple FASTQ"
        options = options.merge({:fastq1_files => fastq_files.first, :fastq2_files => (fastq_files.last || [])})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({:fastq1 => fastq_files.first.first, :fastq2 => fastq_files.last.first})
        {:inputs => options, :jobname => sample}
      end
    elsif uBAM_files = sample_files[:uBAM]
      raise "No support for uBAM"
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        {:task => :BAM_rescore_mutiplex, :inputs => options, :jobname => sample}
      else
        options = options.merge({"HTS#uBAM" => uBAM_files})
        {:inputs => options, :jobname => sample}
      end
    elsif bam_files = sample_files[:BAM]
      raise "No support for BAM"
      options = options.merge({"HTS#BAM" => [bam_files].flatten.first})
      {:inputs => options, :jobname => sample}
    elsif orig_bam_files = sample_files["orig.BAM"]
      raise "No support for orig.BAM"
      if options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign_by_group, :inputs => options, :jobname => sample}
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first})
        {:task => :BAM_rescore_realign, :inputs => options, :jobname => sample}
      end
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample}
    end
  end

end
