module Sample

  def self.load_study_files
    study_files = {}
    Sample.all_studies.each do |study|
      dir = Sample.study_dir study

      sample_files = {}
      dir.glob('WES/*').each do |path|
        file = File.basename path
        sample = file.split(".").first
        fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz")).sort
        if fastq_files.any?
          fastq2_files = fastq_files.reject{|f| File.basename(f) =~ /_2\.fastq/ }
          fastq1_files = fastq_files - fastq2_files
          sample_files[sample] ||= {}
          sample_files[sample]["FASTQ"] = [fastq1_files, fastq2_files]
        end
        
        bam_files = path.glob("*.bam")
        if bam_files.any?
          sample_files[sample] ||= {}
          sample_files[sample]["BAM"] = bam_files
        end

        sample_files[sample] ||= {}
        sample_files[sample]["BAM"] = path if path =~ /.*\.bam/

        ubam_files = path.glob("*.ubam")
        if ubam_files.any?
          sample_files[sample] ||= {}
          sample_files[sample]["uBAM"] = ubam_files
        end
        sample_files[sample] ||= {}
        sample_files[sample]["uBAM"] = path if path =~ /.*\.ubam/

        orig_files = path.glob("*.orig.bam") + path.orig.glob("*.bam")
        if orig_files.any?
          sample_files[sample] ||= {}
          sample_files[sample]["orig.BAM"] = orig_files
        end

        sample_files[sample] ||= {}
        sample_files[sample]["orig.BAM"] = path if path =~ /.*\.orig\.bam/
      end

      dir.glob('WES.orig/*').each do |path|
        file = File.basename path
        sample = file.split(".").first
        orig_files = path.glob("*.bam")
        if orig_files.any?
          sample_files[sample] ||= {}
          sample_files[sample]["orig.BAM"] = orig_files
        end
        sample_files[sample] ||= {}
        sample_files[sample]["orig.BAM"] = path if path =~ /.*\.bam$/
      end

      sample_files.delete_if do |sample,info| info.empty? end
      study_files[study] = sample_files
    end
    study_files
  end

  def self.sample_files(sample)
    load_study_files.each do |study, sample_files|
      files = sample_files[sample]
      return IndiferentHash.setup(files) if files
    end
    nil
  end

  def self.sample_study(sample)
    load_study_files.each do |study, sample_files|
      return study if sample_files[sample]
    end
    nil
  end
  
  
  def self.sample_options(sample)
    load_study_files.each do |study, sample_files|
      next unless Sample.study_dir(study).sample_info.exists?
      sample_info = TSV.open(Sample.study_dir(study).sample_info).to_double
      if sample_info[sample]
        info = sample_info[sample].to_hash
        IndiferentHash.setup info
        info[:sample_name] ||= sample
        return info
      end
    end
    IndiferentHash.setup({:sample_name => sample})
  end

  STUDY_OPTIONS = {:organism => :string, :reference => :string, :interval_list => :file, 
                   :pon => :file, :germline_resource => :file}
  def self.study_options(sample)
    options = {}
    load_study_files.each do |study, sample_files|
      next unless Sample.sample_study(sample) == study
      options_file = Sample.study_dir(study).options.find
      next unless options_file.exists?
      study_options = if File.directory? options_file
                        input_names = STUDY_OPTIONS.keys
                        input_types = STUDY_OPTIONS
                        Workflow.load_inputs(options_file, input_names, input_types)
                      else
                        YAML.load(options_file)
                      end
      options.merge! study_options
    end
    options
  end

  input :by_group, :boolean, "Separate files by read-group if RevertSam is required", true
  extension :bam
  dep_task :BAM, HTS, :BAM_rescore do |sample,options|
    sample_files = Sample.sample_files sample

    options = options.merge(Sample.sample_options(sample))
    options = options.merge(Sample.study_options(sample))
    
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
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise "Normal sample for #{ sample } not found" if sample_files.nil?
    {:inputs => options, :jobname => nsample} if sample_files
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    {:inputs => options, :jobname => sample} if sample_files
  end
  dep_task :mutect2_snv, HTS, :mutect2_clean, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    {:inputs => options, :jobname => sample} if sample_files
  end
  dep_task :strelka, HTS, :strelka, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise ParameterException, "No normal sample found" if sample_files.nil?
    {:inputs => options, :jobname => sample} 
  end
  dep_task :sequenza_purity, HTS, :sequenza_purity, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise ParameterException, "No normal sample found" if sample_files.nil?
    {:inputs => options, :jobname => sample} 
  end
  dep_task :sequenza_ploidy, HTS, :sequenza_ploidy, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise ParameterException, "No normal sample found" if sample_files.nil?
    {:inputs => options, :jobname => sample} 
  end
  dep :sequenza_purity
  dep_task :varscan, HTS, :varscan_somatic_alt, :normal => :BAM_normal, :tumor => :BAM, :tumor_purity => :sequenza_purity do |jobname,options,dependencies|
    if dependencies.flatten.length == 2
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise ParameterException, "No normal sample found" if sample_files.nil?
    {:inputs => options, :jobname => sample} 
  end
  dep_task :delly, HTS, :delly, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end

  dep :BAM
  dep :BAM_normal do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
      break if sample_files
    end

    raise ParameterException, "No normal sample found" if sample_files.nil?
    {:inputs => options, :jobname => sample} 
  end
  dep_task :svABA, HTS, :svABA, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
    if dependencies.flatten.length == 1
      options[:normal] = nil
    end
    {:inputs => options}
  end


end
