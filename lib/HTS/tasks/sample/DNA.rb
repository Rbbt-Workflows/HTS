
module Sample

  def self.can_produce?(job)
    return true if job.done?
    deps = job.dependencies.dup
    while deps.any?
      dep = deps.pop
      return false if dep.task_name == :missing_data
      next if dep.done?
      deps.concat dep.dependencies if dep.dependencies
    end
    true
  end

  task :missing_data => :string do 
    sample = clean_name
    raise "Sample #{ sample } unknown or data missing"
  end

  input :by_group, :boolean, "Separate files by read-group if RevertSam is required", false
  input :bazam, :boolean, "Use bazam instead of RevertSam", false
  extension "bam"
  dep_task :BAM, HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder do |sample,options|
    sample_files = Sample.sample_files sample

    next {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?

    options = add_sample_options sample, options

    if fastq_files = sample_files[:FASTQ]
      if Array === fastq_files.first && fastq_files.first.length > 1
        options = options.merge({:fastq1_files => fastq_files.first, :fastq2_files => (fastq_files.last || [])})
        HTS.job(:BAM_rescore_mutiplex, sample, options)
      else
        options = options.merge({:fastq1 => fastq_files.first, :fastq2 => fastq_files.last})
        {:inputs => options, :jobname => sample}
      end
    elsif uBAM_files = sample_files[:uBAM]
      if Array === uBAM_files && uBAM_files.length > 1
        options = options.merge({:uBAM_files => uBAM_files})
        HTS.job(:BAM_rescore_mutiplex, sample, options)
      else
        uBAM_files = uBAM_files.first if Array === uBAM_files
        options = options.merge({"HTS#uBAM" => uBAM_files, :not_overriden => true})
        HTS.job(:BAM, sample, options)
      end
    elsif bam_files = sample_files[:BAM]
      path = [bam_files].flatten.first
      job = Step.new path
      job.task_name = :BAM
      job.workflow = HTS
      job
    elsif cram_files = sample_files[:CRAM]
      path = [cram_files].flatten.first
      job = Step.new path
      job.task_name = :BAM
      job.workflow = HTS
      job
    elsif orig_bam_files = sample_files["orig.BAM"]
      if options[:bazam]
        options = options.merge({:bam => [orig_bam_files].flatten.first, :not_overriden => true})
        options.delete(:bazam)
        HTS.job(:BAM_rescore_realign_bazam, sample, options)
      elsif options[:by_group]
        options = options.merge({:bam_file => [orig_bam_files].flatten.first, :not_overriden => true})
        HTS.job(:BAM_rescore_realign_by_group, sample, options)
      else
        options = options.merge({:bam_file => [orig_bam_files].flatten.first, :not_overriden => true})
        HTS.job(:BAM_rescore_realign, sample, options)
      end
    elsif orig_cram_files = sample_files["orig.CRAM"]
      if options[:by_group]
        options = options.merge({:bam_file => [orig_cram_files].flatten.first, :not_overriden => true})
        HTS.job(:BAM_rescore_realign_by_group, sample, options)
      else
        options = options.merge({:bam_file => [orig_cram_files].flatten.first, :not_overriden => true})
        HTS.job(:BAM_rescore_realign, sample, options)
      end
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample}
    end
  end

  extension 'bam'
  dep_task :BAM_normal, Sample, :BAM do |sample,options|
    nsample = nil
    sample_files = nil
    sample_study = Sample.sample_study(sample)
    [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample

      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files
    end

    if sample_files
      {:inputs => options, :jobname => nsample}
    else
      j = Sample.job(:BAM, sample + "_normal", options)
      Sample.can_produce?(j) ? j : nil
    end
  end

  #{{{ EXPORT HTS METHODS
  #{{{ genomic_mutations and variant callers
  {
    :strelka => :strelka,
    :mutect2 => [:mutect2, true],
    :mutect2_filters => [:mutect2_filters, true],
    :varscan => :varscan_fpfiltered,
    :somatic_sniper => :somatic_sniper_filtered,
    :muse => :muse,
    :delly => :delly,
    :svABA => :svABA,
    :svABA_indels => :svABA_indels,
    :sequenza_purity => :sequenza_purity,
    :sequenza_ploidy => :sequenza_ploidy,
    :sequenza_CNV => :sequenza_CNV,
    :manta_pre => :manta_pre,
    :manta_somatic => :manta_somatic,
    :pindel_indels => :pindel_indels,
    :haplotype => :haplotype,
    #:somatic_seq => :somatic_seq,
    #:varscan_somatic => :varscan_somatic,
    #:varscan_somatic_alt => :varscan_somatic_alt,
    #:varscan_classify => :varscan_classify,
    #:varscan_fpfiltered => :varscan_fpfiltered,

  }.each do |task,otask|
    if Array === otask
      otask, allow_tumor_only = otask
    else
      allow_tumor_only = false
    end

    if allow_tumor_only
      dep :BAM_normal, :compute => :canfail do |sample,options|
        nsample = nil
        sample_files = nil
        sample_study = Sample.sample_study(sample)
        [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
          break if sample_files
        end

        if sample_files
          {:inputs => options, :jobname => sample}
        else
          j = Sample.job(:BAM_normal, sample, options)
          Sample.can_produce?(j) ? j : nil
        end
      end
    else
      dep :BAM_normal, :compute => :bootstrap do |sample,options|
        nsample = nil
        sample_files = nil
        sample_study = Sample.sample_study(sample)
        [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
          nsample = normal_sample
          sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
          break if sample_files
        end

        if sample_files
          {:inputs => options, :jobname => sample}
        else
          j = Sample.job(:BAM_normal, sample, options)
          Sample.can_produce?(j) ? j : {:workflow => Sample, :task => :missing_data, :jobname => sample}
        end
      end
    end
    dep :BAM, :compute => :bootstrap
    extension :vcf if CALLERS.include?(task.to_s)
    dep_task task, HTS, otask, :normal => :BAM_normal, :tumor => :BAM do |jobname,options,dependencies|
      sample = jobname
      sample_files = Sample.sample_files sample

      if dependencies.empty? && sample_files.nil?
        next {:workflow => Sample, :task => :missing_data, :jobname => sample}
      end

      options = add_sample_options jobname, options

      if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.empty?
        options[:normal] = nil
      end
      {:inputs => options}
    end
  end

  dep :BAM_normal
  dep_task :normal_mutect2_for_panel, HTS, :mutect2_pre, :tumor => :BAM_normal, :normal => nil, :pon => 'none', :max_mnp_distance => 0 do |sample,options,dependencies|
    options = add_sample_options sample, options
    bam_normal = dependencies.flatten.first
    {:inputs => options}
  end



  dep :mutect2, :canfail => true
  dep :strelka, :canfail => true
  dep :muse, :canfail => true
  dep :manta_somatic, :canfail => true
  #dep :pindel_indels, :canfail => true
  #dep :svABA_indels, :canfail => true
  #dep :varscan
  task :caller_cohort => :array do
    dependencies.collect do |dep|
      next if dep.error?
      name = dep.task_name.to_s
      Open.cp(dep.path, file(name))
    end
    Path.setup(files_dir).glob("*")
  end

  dep :caller_cohort, :compute => :produce
  dep Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    cohort = dependencies.flatten.select{|dep| dep.task_name.to_s === "caller_cohort" }.first
    cohort.dependencies.collect{|dep| {:inputs => options.merge(:vcf_file => dep), :jobname => dep.task_name.to_s} }
  end
  task :combined_caller_cohort => :tsv do
    tsv = TSV.setup({}, "Genomic Mutation~Caller#:type=:flat")
    dependencies[1..-1].each do |dep| 
      name = dep.clean_name.to_s
      dep.load.each do |mutation|
        tsv[mutation] ||= []
        tsv[mutation]  << name
      end
    end
    tsv
  end

  dep :caller_cohort, :compute => :produce
  extension :vcf
  task :combined_caller_vcfs => :text do
    list = {}
    step(:caller_cohort).dependencies.each do |dep|
      next if dep.error?
      name = dep.task_name.to_s
      list[name] = dep.path
    end
    
    HTS.combine_caller_vcfs(list)
  end

  dep :combined_caller_vcfs
  input :min_callers, :integer, "Min number of callers to pass variant", 2
  extension :vcf
  task :consensus_somatic_variants => :text do |min_callers|
    TSV.traverse step(:combined_caller_vcfs), :into => :stream, :type => :array do |line|
      next line if line =~ /^#/
      num = line.split("\t")[6].split(";").collect{|f| f.split("--").first}.uniq.length
      next unless num >= min_callers
      line
    end
  end

  dep :BAM_normal do |jobname,options,dependencies|
    sample = jobname
    if Sample.sample_files(sample + "_normal") || (sample.include?(":") && Sample.sample_files(sample.sub(/:.*/, ":normal")))
      {:inputs => options, :jobname => jobname, :task => :BAM_normal}
    elsif sample_files = Sample.sample_files(sample)
      {:inputs => options, :jobname => jobname, :task => :BAM}
    elsif Sample.can_produce?(job = Sample.job(:BAM_normal, jobname))
      job
    elsif Sample.can_produce?(job = Sample.job(:BAM, jobname)).done?
      job
    else
      {:workflow => Sample, :task => :missing_data, :jobname => sample} if sample_files.nil?
    end
  end
  dep_task :haplotype, HTS, :haplotype, :BAM => :BAM_normal do |jobname, options, dependencies|
    options = add_sample_options jobname, options
    options[:BAM] = :BAM if dependencies.flatten.reject{|dep| (dep.dependencies.empty? ) || (dep.error? && ! dep.recoverable_error?) }.select{|dep|  dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname}
  end

  dep :strelka, :compute => [:bootstrap, :canfail] do |jobname, options|
    %w(strelka varscan mutect2 somatic_sniper muse).collect do |var_caller|
      {:task => var_caller, :inputs => options, :jobname => jobname}
    end
  end
  input :only_pass, :boolean, "Consider only PASS variants", true
  task :caller_comparison => :tsv do |only_pass|
    caller_variants = {}
    dependencies.each do |dep|
      next if dep.error?
      var_caller = dep.task_name

      variants = TSV.traverse dep, :into => [], :type => :list do |chr, values|
        pos, _id, _ref, alt, _qual, filter, *_rest = values

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

    if sample_files && sample_files.include?("VCF")
      nil
    else
      vcaller = options[:caller]
      options[:vcf_file] = vcaller.to_sym
      {:task => vcaller, :jobname => sample, :inputs => options}
    end
  end
  extension "vcf"
  task :vcf_file => :text do
    sample_files = Sample.sample_files self.clean_name
    if sample_files && sample_files.include?("VCF")
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
  dep :organism
  input :positions, :array, "List of positions to image. Use list of genomic mutations by default", :genomic_mutations
  dep_task :mutation_BAM_img, HTS, :mutation_BAM_img, :organism => :organism, :tumor => :BAM, :normal => :BAM_normal do |jobname,options|
    options = add_sample_options jobname, options
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
    sample = jobname
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

  dep :vcf_file
  dep_task :expanded_vcf, Sequence, :expanded_vcf, :vcf_file => :vcf_file


  dep :BAM
  dep :BAM_normal
  dep_task :plot_freec_results, HTS, :plot_freec_results, :sample_mateFile => :BAM, :control_mateFile => :BAM_normal

  dep Sample, :BAM
  dep Sample, :BAM_normal
  dep Sample, :mutect2
  dep Sample, :somatic_sniper
  dep Sample, :varscan
  #dep Sample, :muse
  dep Sample, :strelka
  dep_task :somatic_seq, HTS, :somatic_seq, :mutect2_vcf => :mutect2, :varscan_snv => :varscan, :somaticsniper_vcf => :somatic_sniper, :strelka_vcf => :strelka do |jobname, options,deps|
    tumor = deps.flatten.select{|d| d.task_name == :BAM}.first.path
    normal = deps.flatten.select{|d| d.task_name == :BAM_normal}.first.path

    options = options.merge({:tumor_bam_file => tumor, :normal_bam_file => normal})
    {:inputs => options}
  end

  dep :somatic_seq
  task :somatic_seq_filtered do |jobname, options, deps|
    reference = HTS.helpers[:reference_file].call(step(:somatic_seq).dependencies.flatten.select{|d| d.inputs.to_hash.include?"genome_reference"}.first.inputs.to_hash["genome_reference"])
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    args = {}
    args["reference"] = reference
    args["variant"] = step(:somatic_seq).path
    args["output"] = self.tmp_path
    args["exclude-filtered"] = true
    args["exclude-non-variants"] = true
    GATK.run_log("SelectVariants", args)
    FileUtils.rm_rf self.path
    FileUtils.ln_s self.tmp_path, self.path
  end

  dep Sample, :BAM
  dep Sample, :BAM_normal
  dep Sample, :mutect2
  dep Sample, :somatic_sniper
  dep Sample, :varscan
  #dep Sample, :muse
  dep Sample, :strelka
  dep_task :somatic_seq, HTS, :somatic_seq, :tumor_bam_file => :BAM, :normal_bam_file => :BAM_normal, :mutect2_vcf => :mutect2, :varscan_snv => :varscan, :somaticsniper_vcf => :somatic_sniper, :strelka_vcf => :strelka


  dep Sample, :BAM
  dep Sample, :BAM_normal
  dep Sample, :mutect2
  dep Sample, :somatic_sniper
  dep Sample, :varscan
  #dep Sample, :muse
  dep Sample, :strelka
  dep_task :somatic_seq_prefilter, HTS, :somatic_seq, :tumor_bam_file => :BAM, :normal_bam_file => :BAM_normal, :mutect2_vcf => :mutect2, :somaticsniper_vcf => :somatic_sniper do |jobname,options,dependencies|

    varscan = dependencies.flatten.select{|dep| dep.task_name == :varscan}.first
    options[:varscan_snv] = varscan.step(:varscan_somatic).file('output')["Default" + '.snv.vcf']
    options[:varscan_indel] = varscan.step(:varscan_somatic).file('output')["Default" + '.indel.vcf']

    strelka = dependencies.flatten.select{|dep| dep.task_name == :strelka}.first
    options[:strelka_snv] = strelka.step(:strelka_pre).path
    options[:strelka_indel] = strelka.step(:strelka_filtered_indels).step(:strelka_filtered).dependencies.first

    {:inputs => options}
  end

  input :studies, :array, "Samples to use", Sample.all_studies
  dep :mutect2_filters, :compute => :produce, :max_mnp_distance => 0 do |jobname, options|
    samples = []
    options[:studies].each do |study|
      samples += (Sample.study_samples(study))
    end
    normal_samples = samples.flatten.select{|x| x.include? "normal" }.uniq
    normal_samples.collect {|normal_sample|
      {:jobname => normal_sample, :inputs => options}
    }
  end
  extension :vcf
  task :generatePON => :text  do
    intervals = dependencies.first.recursive_inputs[:interval_list]
    orig_reference = dependencies.first.recursive_inputs[:reference]
    orig_reference = HTS.helpers[:reference_file].call(orig_reference)

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    bed_dir = file('bed')
    work_dir = file('work')
    args = {}
    args["variant"] = []

    Misc.in_dir bed_dir do
      dependencies.each do |dep|
        vcf = GATK.sort_VCF(dep.path, bed_dir)
        gvcf = HTS.prepare_BED(vcf, bed_dir)
        args["variant"].push(gvcf)
      end
    end

    Misc.in_dir work_dir do
      args["reference"] = reference
      args["genomicsdb-workspace-path"] = "pon_db"
      args["merge-input-intervals"] = "TRUE"
      args["intervals"] = intervals
      GATK.run_log("GenomicsDBImport", args)

      args = {}
      args["reference"] = reference
      args["variant"] = "gendb://pon_db"
      args["output"] = self.tmp_path
      GATK.run_log("CreateSomaticPanelOfNormals", args)
      Open.rm_rf bed_dir
      FileUtils.ln_s self.tmp_path, self.path
    end
  end

  dep :BAM
  dep_task :coverage, HTS, :genomecov do |sample,options|
    options = add_sample_options sample, options
    options[:BAM] = :BAM
    {:inputs => options, :jobname => sample}
  end

  dep :coverage
  dep_task :low_coverage, HTS, :low_coverage do |sample,options,deps|
    options = add_sample_options sample, options
    options[:coverage_file] = deps.first.path
    {:inputs => options}
  end

  #dep Sample, :BAM
  #dep_task :collect_fragment_counts, HTS, :collect_fragment_counts, :bam => :BAM
  dep Sample, :BAM
  dep Sample, :BAM_normal
  dep_task :manta, HTS, :manta_pre, :tumor => :BAM, :normal => :BAM_normal do |sample,options,deps|
    options = add_sample_options sample, options
    {:inputs => options}
  end

  dep :BAM, :compute => [:produce, :canfail]
  dep :BAM_normal, :compute => [:produce, :canfail]
  dep_task :intervals_from_BAM, HTS, :intervals_from_BAM do |sample,options,dependencies|
    options = add_sample_options sample, options
    bam_normal = dependencies.flatten.select{|dep| dep.task_name.to_s == "BAM_normal"}.first
    if bam_normal.error? && ! bam_normal.recoverable_error?
      options[:BAM] = :BAM 
    else
      options[:BAM] = :BAM_normal 
    end
    {:inputs => options, :jobname => sample}
  end

  # ToDo: Utilizar todas las mutaciones de la cohorte cuando haya mas de una
  # muestra
  #dep :BAM
  #dep :genomic_mutations
  #dep :sequenza_CNV
  #dep_task :pyclone, HTS, :pyclone, :copynumber => :sequenza_CNV, :bam => :BAM, :mutations => :genomic_mutations do |jobname,options|
  #  options[:sample_name] = jobname
  #  {:inputs => options}
  #end

  dep :BAM
  dep :genomic_mutations
  dep :sequenza_CNV
  dep :mutect2
  dep :sequenza_purity
  dep_task :pyclone, HTS, :pyclone, :copynumber => :sequenza_CNV, :bam => :BAM, :vcf_file => :mutect2, :mutations => :genomic_mutations, :tumor_content => :sequenza_purity do |jobname,options|
    options[:sample_name] = jobname
    {:inputs => options}
  end

  dep :BAM
  dep :genomic_mutations
  dep :sequenza_CNV
  dep :mutect2
  dep :sequenza_purity
  dep :organism
  dep_task :cancer_cell_fraction, HTS, :cancer_cell_fraction, :copynumber => :sequenza_CNV, :bam => :BAM, :vcf_file => :mutect2, :mutations => :genomic_mutations, :tumor_content => :sequenza_purity, :organism => :organism do |jobname,options|
    options[:sample_name] = jobname
    {:inputs => options}
  end
end
