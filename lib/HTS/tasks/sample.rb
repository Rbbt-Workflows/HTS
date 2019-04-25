module Sample

  def self.load_study_files
    @@study_files ||= begin
                        study_files = {}
                        Sample.all_studies.each do |study|
                          dir = Sample.study_dir study

                          sample_files = {}
                          dir.glob('WES/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz")).sort
                            if fastq_files.any?
                              fastq2_files = fastq_files.select{|f| File.basename(f) =~ /_2\.fastq/ }
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

  input :by_group, :boolean, "Separate files by read-group if RevertSam is required", false
  extension :bam
  dep_task :BAM, HTS, :BAM do |sample,options|
    sample_files = Sample.sample_files sample
    raise "Sample #{ sample } not found" if sample_files.nil?

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

  #{{{ EXPORT HTS METHODS

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
        [sample + '_normal', 'normal'].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
          break if sample_files
        end

        {:inputs => options, :jobname => sample} if sample_files
      end
    else
      dep :BAM, :compute => :bootstrap
      dep :BAM_normal, :compute => :bootstrap do |sample,options|
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
    end
    extension :vcf if CALLERS.include?(task.to_s)
    dep_task task, HTS, otask, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
      if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.empty?
        options[:normal] = nil
      end
      {:inputs => options}
    end
  end

  dep :BAM
  dep_task :haplotype, HTS, :haplotype, :BAM => :BAM


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
  dep :mutect2 do |jobname,options|
    vcaller = options[:caller]
    {:task => vcaller, :jobname => jobname, :inputs => options}
  end
  dep_task :genomic_mutations, Sequence, :genomic_mutations, :vcf_file => :mutect2 do |jobname,options|
    vcaller = options[:caller]
    options[:vcf_file] = vcaller.to_sym
    {:task => :genomic_mutations, :workflow => Sequence, :jobname => jobname, :inputs => options}
  end

  dep :BAM, :compute => :bootstrap
  dep :BAM_normal, :compute => :bootstrap do |sample,options|
    nsample = nil
    sample_files = nil
    [sample + '_normal', 'normal'].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if Sample.sample_study(sample) == Sample.sample_study(nsample)
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
    options = options.merge(Sample.study_options(jobname))
    {:inputs => options}
  end

end
