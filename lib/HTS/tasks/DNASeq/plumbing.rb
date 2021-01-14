require 'tools/BAM_shard'
module HTS

  input :bam_files, :array, "BAM filenames to multiplex"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :bam
  task :BAM_multiplex => :binary do |bam_filenames,reference|
    bam_filenames = Dir.glob(File.join(bam_filenames.first, "*.bam")) if Array === bam_filenames && bam_filenames.length == 1 && ! Step === bam_filenames && File.directory?(bam_filenames.first)
    bam_filenames = bam_filenames.collect{|f| Step === f ? f.path : f}

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    unsorted = file('unsorted.bam')

    args= {}
    args["INPUT"] = bam_filenames
    args["OUTPUT"] = unsorted
    args["METRICS_FILE"] = file('metrics.txt')
    args["ASSUME_SORT_ORDER"] = 'queryname'
    gatk("MarkDuplicates", args)
    
    sort = ! info[:spark] 

    if sort
      sorted = file('sorted.bam')
      args= {}
      args["INPUT"] = unsorted
      args["OUTPUT"] = sorted
      args["SORT_ORDER"] = 'coordinate'
      gatk("SortSam", args)
      FileUtils.rm unsorted
      Open.mv sorted, self.path
    else
      Open.mv unsorted, self.path
    end
    nil
  end

  #input :ubam_files, :array, "uBAM filenames to multiplex"
  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  #extension :bam
  #task :BAM_multi_bwa => :binary do |ubam_filenames,reference|
  #  args= {}
  #  ubam_filenames = Dir.glob(File.join(ubam_filenames.first, "*.bam") + File.join(ubam_filenames.first, "*.ubam") ) if Array === ubam_filenames && ubam_filenames.length == 1 && File.directory?(ubam_filenames.first)
  #  ubam_filenames = ubam_filenames.collect{|f| Step === f ? f.path : f}
  #  output = file('out.bam')

  #  args["INPUT"] = ubam_filenames
  #  args["OUTPUT"] = output
  #  args["METRICS_FILE"] = file('metrics.txt')
  #  args["ASSUME_SORT_ORDER"] = 'queryname'

  #  FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
  #  gatk("MarkDuplicates", args)
  #  
  #  Open.mv output, self.path
  #  nil
  #end


  input :bam_file, :binary, "Bam file", nil, :nofile => true
  input :by_group, :boolean, "Separate files by read-group", false
  input :max_discard_fraction, :boolean, "Max dicard fraction", 0.05
  input :tmp_dir, :string, "Temporary directory", nil
  task :revert_BAM => :binary do |bam_file,by_group,max_discard_fraction,tmp_dir|
    args = {}
    args["INPUT"] = bam_file

    if by_group
      Open.mkdir file("uBAM").find
      args["OUTPUT"] = file("uBAM").find
      args["OUTPUT_BY_READGROUP"] = "true"
    else
      args["OUTPUT"] = self.tmp_path
      args["OUTPUT_BY_READGROUP"] = "false"
    end
    
    args["SANITIZE"] = "true"
    args["SORT_ORDER"] = "queryname"
    args["KEEP_FIRST_DUPLICATE"] = "true"

    args["MAX_DISCARD_FRACTION"] = max_discard_fraction
    args["ATTRIBUTE_TO_CLEAR"] = ["XA", "BD", "XS", "BI"]
    args["RESTORE_ORIGINAL_QUALITIES"] = "true"
    args["REMOVE_DUPLICATE_INFORMATION"] = "true"
    args["REMOVE_ALIGNMENT_INFORMATION"] = "true"

    gatk("RevertSam", args, tmp_dir)
    if by_group
      file("uBAM").glob("*")
    end
  end

  input :bam_file , :file, "BAM", :nofile => true
  task :revert_BAM_sharded do |bam_file|
    bam_file = Samtools.prepare_BAM(bam_file)
    BAMShard.revert_BAM(bam_file, self.path)
  end

  input :fastq1_files, :array, "FASTQ files for first mate"
  input :fastq2_files, :array, "FASTQ files for second mate", []
  input :uBAM_files, :array, "uBAM files for second mate", []
  dep :BAM_bwa, :compute => :bootstrap do |jobname,options|
    uBAM_files = options[:uBAM_files]
    fastq1_files = options[:fastq1_files]
    if fastq1_files
      fastq2_files = options[:fastq2_files]
      fastq1_files.zip(fastq2_files).collect do |fastq1,fastq2|
        read_group_name = File.basename(fastq1).sub(/(_{1,2})?\.fastq.*/,'')
        options = options.merge({:fastq1 => fastq1, :fastq2 => fastq2, :read_group_name => read_group_name})
        {:inputs => options, :jobname => [jobname, read_group_name] * "." }
      end
    elsif uBAM_files
      uBAM_files = Dir.glob(File.join(uBAM_files, '*.bam')) + Dir.glob(File.join(uBAM_files, '*.ubam')) if (String === uBAM_files && File.directory?(uBAM_files))
      [uBAM_files].flatten.collect do |uBAM|
        read_group_name = File.basename(uBAM).sub(/.u?bam/i,'')
        options = options.merge({"HTS#uBAM" => uBAM})
        {:inputs => options, :jobname => [jobname, read_group_name] * "." }
      end
    else
      raise "No FASTQ or uBAM files for #{ jobname }"
    end
  end
  dep :BAM_multiplex do |jobname, options,dependencies|
    bam_files = dependencies.flatten.collect{|dep| dep }
    {:jobname => jobname, :inputs => options.merge(:bam_files => bam_files)}
  end
  extension :bam
  dep_task :BAM_rescore_mutiplex, HTS, :BAM_rescore do |jobname,options, dependencies|
    mutiplex = dependencies.flatten.select{|dep| dep.task_name == :BAM_multiplex}.first
    {:inputs => options.merge("HTS#BAM_sorted" =>  mutiplex), :jobname => jobname}
  end

  dep :revert_BAM, :compute => :produce
  dep :BAM_bwa, :compute => :produce, 
    :fastq1 => :placeholder, :fastq2 => :placeholder,
    :sample_name => :placeholder,
    :read_group_name => :placeholder,
    :library_name => :placeholder,
    :platform_unit => :placeholder,
    :platform => :placeholder,
    :sequencing_center => :placeholder do |jobname, options, dependencies|

    read_groups = CMD.cmd("samtools view --no-PG -H #{options[:bam_file]}").read.split("\n").select{|line| 
      line =~ /^@RG/
    }.collect{|line| 
      line.split("\t").select{|part| part =~ /^ID:(.*)/}.first.split(":").last
    }

    read_groups.collect do |read_group|
      uBAM = Step.new dependencies.first.file('uBAM')[read_group] + ".bam"
      #uBAM.dependencies << dependencies.first
      #{:task => :BAM_bwa, :inputs => options.merge({"HTS#uBAM" => uBAM}), :jobname => [jobname, read_group] * "."}
      job = HTS.job(:BAM_bwa, [jobname, read_group] * ".",  options.merge({"HTS#uBAM" => uBAM}))
      job.step(:mark_adapters).dependencies += [dependencies.flatten.first]
      job
    end
  end
  dep :BAM_multiplex do |jobname, options,dependencies|
    bam_files = dependencies.flatten.select{|dep| dep.task_name == :BAM_bwa}
    {:jobname => jobname, :inputs => options.merge(:bam_files => bam_files)}
  end
  extension :bam
  dep_task :BAM_rescore_realign_by_group, HTS, :BAM_rescore do |jobname,options, dependencies|
    mutiplex = dependencies.flatten.select{|dep| dep.task_name == :BAM_multiplex}.first
    {:inputs => options.merge("HTS#BAM_duplicates" =>  mutiplex), :jobname => jobname}
  end

  dep :revert_BAM, :compute => :produce
  extension :bam
  dep_task :BAM_rescore_realign, HTS, :BAM_rescore do |jobname,options,dependencies|
    {:inputs => options.merge("HTS#uBAM" =>  dependencies.first), :jobname => jobname}
  end
end
