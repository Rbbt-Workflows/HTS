module Sample

  STUDY_OPTIONS = {
    :organism => :string,
    :reference => :string,
    :type_of_sequencing => :string,
    :interval_list => :file,
    :pon => :file,
    :germline_resource => :file,
    :skip_rescore => :boolean,
    :skip_duplicates => :boolean,
    :remove_unpaired => :boolean,
    :remove_soft_clip => :boolean,
  }

  class << self
    alias original_task_info task_info

    def task_info(*args)
      info = original_task_info(*args)
      defaults = IndiferentHash.setup(info[:input_defaults])
      STUDY_OPTIONS.keys.each{|input| defaults.delete input}
      info[:input_defaults] = defaults
      info
    end
  end

  def self.load_study_files_DNA
    @@study_files_DNA ||= begin
                        study_files = {}
                        Sample.all_studies.each do |study|
                          dir = Sample.study_dir study
                          sample_files = {}

                          dir.glob('genotypes/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first

                            if File.directory?(path)
                              vcf_files = path.glob("*.vcf.*")
                              if vcf_files.any?
                                sample_files[sample] ||= {}
                                sample_files[sample]["VCF"] = vcf_files
                              end
                            else
                              if path =~ /\.vcf*/
                                sample_files[sample] ||= {}
                                sample_files[sample]["VCF"] ||= []
                                sample_files[sample]["VCF"] << path
                              end
                            end
                          end

                          dir.glob('W?S/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz") + path.glob("*.fq") + path.glob("*.fq.gz")).sort
                            if fastq_files.any?
                              fastq2_files = fastq_files.select{|f| File.basename(f) =~ /(?:\.|_)(?:2|reads?2)\.(?:fastq|fq)/ }
                              fastq2_files += fastq_files.select{|f| File.basename(f) =~ /_R2_/ }
                              fastq1_files = fastq_files - fastq2_files
                              sample_files[sample] ||= {}
                              sample_files[sample]["FASTQ"] = [fastq1_files, fastq2_files]
                            end

                            bam_files = path.glob("*.bam") + path.glob("*.BAM")
                            if bam_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["BAM"] = bam_files
                            end

                            sample_files[sample] ||= {}
                            sample_files[sample]["BAM"] = path if path =~ /.*\.bam/i

                            ubam_files = path.glob("*.ubam") + path.glob("*.uBAM")
                            if ubam_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["uBAM"] = ubam_files
                            end
                            sample_files[sample] ||= {}
                            sample_files[sample]["uBAM"] = path if path =~ /.*\.ubam/i

                            orig_files = path.glob("*.orig.bam") + path.orig.glob("*.bam") + \
                              path.glob("*.orig.BAM") + path.orig.glob("*.BAM")
                            if orig_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["orig.BAM"] = orig_files
                            end

                            sample_files[sample] ||= {}
                            sample_files[sample]["orig.BAM"] = path if path =~ /.*\.orig\.bam/i
                          end

                          dir.glob('W?S.orig/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            orig_files = path.glob("*.bam") + path.glob("*.BAM")
                            if orig_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["orig.BAM"] = orig_files
                            end
                            sample_files[sample] ||= {}
                            sample_files[sample]["orig.BAM"] = path if path =~ /.*\.bam$/i
                          end

                          sample_files.delete_if do |sample,info| info.empty? end
                          study_files[study] = sample_files
                        end
                        study_files
                      end
  end

  def self.load_study_files_RNA
    @@study_files_RNA ||= begin
                        study_files = {}
                        Sample.all_studies.each do |study|
                          dir = Sample.study_dir study

                          sample_files = {}
                          dir.glob('RNA/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first

                            fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz") + path.glob("*.fq") + path.glob("*.fq.gz")).sort
                            if fastq_files.any?
                              fastq2_files = fastq_files.select{|f| File.basename(f) =~ /(?:\.|_)(?:2|reads?2)\.(?:fastq|fq)/ }
                              fastq1_files = fastq_files - fastq2_files
                              sample_files[sample] ||= {}
                              sample_files[sample]["RNA_FASTQ"] = [fastq1_files, fastq2_files]
                            end

                            bam_files = path.glob("*.bam")
                            if bam_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["RNA_BAM"] = bam_files
                            end

                            sample_files[sample] ||= {}
                            sample_files[sample]["RNA_BAM"] = path if path =~ /.*\.bam/

                            orig_files = path.glob("*.orig.bam") + path.orig.glob("*.bam")
                            if orig_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["RNA_orig.BAM"] = orig_files
                            end

                            sample_files[sample] ||= {}
                            sample_files[sample]["RNA_orig.BAM"] = path if path =~ /.*\.orig\.bam/
                          end

                          dir.glob('RNA.orig/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            orig_files = path.glob("*.bam")
                            if orig_files.any?
                              sample_files[sample] ||= {}
                              sample_files[sample]["RNA_orig.BAM"] = orig_files
                            end
                            sample_files[sample] ||= {}
                            sample_files[sample]["RNA_orig.BAM"] = path if path =~ /.*\.bam$/
                          end

                          sample_files.delete_if do |sample,info| info.empty? end
                          study_files[study] = sample_files
                        end
                        study_files
                      end
  end

  def self.load_study_files
    @@study_files ||= Persist.persist("Study_files", :yaml, :file => Rbbt.var.cache.HTS_study_files.find, :update => true) do
                        dna = load_study_files_DNA
                        rna = load_study_files_RNA

                        study_files = {}
                        dna.each do |study,sample_files|
                          study_files[study] ||= {}
                          sample_files.each do |sample,info|
                            study_files[study][sample] ||= {}
                            study_files[study][sample].merge!(info)
                          end
                        end

                        rna.each do |study,sample_files|
                          study_files[study] ||= {}
                          sample_files.each do |sample,info|
                            study_files[study][sample] ||= {}
                            study_files[study][sample].merge!(info)
                          end
                        end

                        study_files
                      end
  end

  def self.sample_files(sample)
    study_files = load_study_files
    if sample.include?(":")
      study, _sep, ssample  = sample.partition(":")
      files = study_files[study][ssample] if study_files[study]
      return IndiferentHash.setup(files) if files
    end
    study_files.each do |study, sample_files|
      files = sample_files[sample]
      return IndiferentHash.setup(files) if files
    end
    nil
  end

  def self.study_samples(study)
    study_files = load_study_files
    samples = study_files.find {|s| s.include? study}
    samples = samples.select{|smpl| smpl.is_a? Hash}
    samples = samples.map{|hash| hash.keys}.flatten
    samples.map{|sample| study + ":" + sample}
  end

  def self.sample_study(sample)
    study_files = load_study_files
    if sample.include?(":")
      study = sample.split(":").first
      return study if study_files.include?(study)
    end
    study_files.each do |study, sample_files|
      return study if sample_files[sample]
    end
    nil
  end


  def self.sample_options(sample, sstudy = nil)
    if sample.include?(":") and sstudy.nil?
      study, _sep, ssample  = sample.partition(":")
      return sample_options(ssample, study)
    end

    load_study_files.each do |study, sample_files|
      next if sstudy && study.to_s != sstudy.to_s
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


  def self.study_options(sample, sstudy = nil)
    if sample.include?(":") and sstudy.nil?
      study, _sep, ssample  = sample.partition(":")
      return study_options(ssample, study)
    end

    @@study_options ||= {}
    @@study_options[[study,sstudy]] ||= 
      begin
        options = {}
        IndiferentHash.setup options
        load_study_files.each do |study, sample_files|
          next if sstudy && study.to_s != sstudy.to_s
          next unless sstudy || Sample.sample_study(sample) == study
          options_file = Sample.study_dir(study).options.find
          next unless options_file.exists?
          study_options = if File.directory? options_file
                            input_names = STUDY_OPTIONS.keys
                            input_types = STUDY_OPTIONS
                            options = Workflow.load_inputs(options_file, input_names, input_types)
                            Dir.glob(File.join(options_file, "*#*")).each do |od|
                              options[File.basename(od)] = od
                            end
                            options 
                          else
                            YAML.load(options_file)
                          end
          options.merge! study_options
        end
        options
      end
  end

  def self.add_sample_options(sample, options)
    @@sample_options ||= {}
    sample_options = @@sample_options[sample] ||= 
      begin
        Sample.study_options(sample).merge(Sample.sample_options(sample))
      end
    options = sample_options.merge(options)
    load_sample_workflow(sample)
    IndiferentHash.setup options
    options
  end

  def self.load_study_workflow(study)
    wf_file = Sample.study_dir(study)["workflow.rb"]
    require wf_file.find if wf_file.exists?
  end

  def self.load_sample_workflow(sample)
    study = sample_study(sample)
    load_study_workflow(study)
  end

  input :organism, :string, "Organism Code", nil
  task :organism => :string do |organism|
    sample = clean_name
    options = Sample.sample_options(sample).merge(Sample.study_options(sample))
    organism || options[:organism] || Organism.organism_for_build(options[:reference] || 'b37') || Organism.default_code("Hsa")
  end

  CALLERS = %w(strelka varscan mutect2 muse somatic_sniper delly svABA)
end

require 'HTS/tasks/sample/DNA'
require 'HTS/tasks/sample/RNA'
