module Sample

  def self.load_study_files_DNA
    @@study_files_DNA ||= begin
                        study_files = {}
                        Sample.all_studies.each do |study|
                          dir = Sample.study_dir study

                          sample_files = {}
                          dir.glob('WES/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz")).sort
                            if fastq_files.any?
                              fastq2_files = fastq_files.select{|f| File.basename(f) =~ /(?:_2|_reads2)\.fastq/ }
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

  def self.load_study_files_RNA
    @@study_files_RNA ||= begin
                        study_files = {}
                        Sample.all_studies.each do |study|
                          dir = Sample.study_dir study

                          sample_files = {}
                          dir.glob('RNA/*').each do |path|
                            file = File.basename path
                            sample = file.split(".").first
                            fastq_files = (path.glob("*.fastq") + path.glob("*.fastq.gz")).sort
                            if fastq_files.any?
                              fastq2_files = fastq_files.select{|f| File.basename(f) =~ /(?:_2|_reads2)\.fastq/ }
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
    @@study_files ||= begin
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

  STUDY_OPTIONS = {:organism => :string, :reference => :string, :interval_list => :file, 
                   :pon => :file, :germline_resource => :file}
  def self.study_options(sample, sstudy = nil)
    if sample.include?(":") and sstudy.nil?
      study, _sep, ssample  = sample.partition(":")
      return study_options(ssample, study)
    end

    options = {}
    load_study_files.each do |study, sample_files|
      next if sstudy && study.to_s != sstudy.to_s
      next unless sstudy || Sample.sample_study(sample) == study
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

  def self.add_sample_options(sample,options)
    options = Sample.sample_options(sample).merge(options)
    options = Sample.study_options(sample).merge(options)
    IndiferentHash.setup options
    options
  end

  task :organism => :string do
    sample = clean_name
    options = Sample.sample_options(sample).merge(Sample.study_options(sample))
    options[:organism] || Organism.organism_for_build(options[:reference] || 'b37') || Organism.default_code("Hsa")
  end

end

require 'HTS/tasks/sample/DNA'
require 'HTS/tasks/sample/RNA'
