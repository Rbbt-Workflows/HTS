require 'rbbt/sources/organism'

module HTS

  input :fastq1, :file, "FASTQ 1 file"
  input :fastq2, :file, "FASTQ 2 file"
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil, :jobname => true
  input :library_name, :string, "LIBRARY_NAME BAM field", "DefaultLibraryName"
  input :platform_unit, :string, "PLATFORM_UNIT BAM field", "DefaultPlatformUnit"
  input :platform, :string, "PLATFORM BAM field", "DefaultPlatform"
  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field", "DefaultSequencingCenter"
  extension :ubam
  task :uBAM => :binary do |fastq1, fastq2, read_group_name, sample_name, library_name, platform_unit, platform, sequencing_center, run_date|
    sample_name ||= self.clean_name
    read_group_name ||= sample_name

    args = {}
    args["FASTQ"] = fastq1
    args["FASTQ2"] = fastq2
    args["OUTPUT"] = self.tmp_path
    args["READ_GROUP_NAME"] = read_group_name
    args["SAMPLE_NAME"] = sample_name
    args["LIBRARY_NAME"] = library_name
    args["PLATFORM_UNIT"] = platform_unit
    args["PLATFORM"] = platform
    args["SEQUENCING_CENTER"] = sequencing_center
    args["RUN_DATE"] = run_date

    args[:progress_bar] = gatk_read_count_monitor("FastqToSam") do |bar|
      set_info :reads, bar.ticks
    end

    gatk("FastqToSam", args)

    nil
  end

  dep :uBAM, :compute => :produce
  extension :ubam
  task :mark_adapters => :binary do
    args = {}
    args["INPUT"] = step(:uBAM).path
    args["METRICS"] = file('metrics.txt')
    args["OUTPUT"] = self.tmp_path

    max = step(:uBAM).info[:reads]
    args[:progress_bar] = gatk_read_count_monitor("MarkIlluminaAdapters", max) do |bar|
      ticks = bar.ticks
      ticks = max if ticks.to_i == 0
      set_info :reads, ticks
    end

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    gatk("MarkIlluminaAdapters", args)
  end

  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :bwa_mem_args, :string, "Arg string", "-M -p -K 1000000 -v 0"
  input :remove_unpaired, :boolean, "Remove possible unpaired reads", false
  input :skip_mark_adapters, :boolean, "Skip mark adapters", false
  dep :uBAM
  dep :mark_adapters do |jobname,options,dependencies|
    if options[:skip_mark_adapters]
      nil
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  extension :bam
  task :BAM_bwa => :binary do |reference, bwa_mem_args,remove_unpaired,skip_mark_adapters|

    orig_reference = reference_file(reference)
    #reference = BWA.prepare_FASTA orig_reference
    #reference = GATK.prepare_FASTA orig_reference
    #reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.prepare_FASTA orig_reference

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    fsam_to_fastq = file('SamToFastq.fastq')
    fbwa_bam = file('bwa.bam')
    ffilter_sam = file('FilterSam')
    Open.rm fsam_to_fastq
    Open.rm fbwa_bam
    Open.rm ffilter_sam
    Misc.with_fifo(fsam_to_fastq) do |s2f_path|
      Misc.with_fifo(fbwa_bam) do |bwa_bam|
        Misc.with_fifo(ffilter_sam) do |filter_sam|

          mark_adapters = if skip_mark_adapters
                            step(:uBAM)
                          else
                            step(:mark_adapters)
                          end

          if remove_unpaired
            io_filter = CMD.cmd(:samtools, "view -h --no-PG -f 0x1 '#{mark_adapters.path}'", :pipe => true)

            Misc.consume_stream io_filter, true, filter_sam
          else
            filter_sam = mark_adapters.path
          end

          args = {}
          args["INPUT"] = filter_sam
          args["FASTQ"] = "/dev/stdout"
          args["CLIPPING_ATTRIBUTE"] = "XT"
          args["CLIPPING_ACTION"] = "2"
          args["INTERLEAVE"] = "true"

          io = gatk_io("SamToFastq", args)

          bwa_mem_args += " -t " << (config('cpus', 'bwa', :default => 8) || "1").to_s.strip

          Thread.new do
            BWA.mem_pipe([s2f_path], reference, bwa_mem_args, io, bwa_bam)
          end

          args = {}
          args["ALIGNED_BAM"] = bwa_bam
          args["UNMAPPED_BAM"] = step('uBAM').path
          args["OUTPUT"] = self.tmp_path
          args["R"] = reference
          args["CREATE_INDEX"] = "false"
          args["ADD_MATE_CIGAR"] = "true"
          args["CLIP_ADAPTERS"] = "false"
          args["CLIP_OVERLAPPING_READS"] = "true"
          args["INCLUDE_SECONDARY_ALIGNMENTS"] = "true"
          args["MAX_INSERTIONS_OR_DELETIONS"] = "-1"
          args["PRIMARY_ALIGNMENT_STRATEGY"] = "MostDistant"
          args["ATTRIBUTES_TO_RETAIN"] = "XS"
          args["SORT_ORDER"] = "queryname"

          # FROM GOAST
          args["VALIDATION_STRINGENCY"] = "SILENT"
          args["MAX_RECORDS_IN_RAM"]=200_000_000
          args["EXPECTED_ORIENTATIONS"] = "FR"
          args["ATTRIBUTES_TO_RETAIN"] = "X0"
          args["ATTRIBUTES_TO_REMOVE"] = %w(NM MD)
          args["PAIRED_RUN"] = true
          args["IS_BISULFITE_SEQUENCE"] = false
          args["UNMAPPED_READ_STRATEGY"]="COPY_TO_TAG"
          args["ALIGNER_PROPER_PAIR_FLAGS"]=true 
          args["UNMAP_CONTAMINANT_READS"]=true 
          args["ADD_PG_TAG_TO_READS"]=false

          max = mark_adapters.info[:reads]
          args[:progress_bar] = gatk_read_count_monitor("BWA", max) do |bar|
            ticks = bar.ticks
            ticks = max if ticks.to_i == 0
            set_info :reads, ticks
          end

          gatk("MergeBamAlignment", args)

          FileUtils.rm_rf filter_sam if remove_unpaired
        end
        FileUtils.rm_rf bwa_bam
      end
      FileUtils.rm_rf s2f_path
    end
    nil
  end

  dep :BAM_bwa
  extension :bam
  task :BAM_duplicates => :binary do 
    Open.mkdir files_dir 
    output = file('out.bam')

    args = {}
    args["INPUT"] = step(:BAM_bwa).path
    args["OUTPUT"] = output
    args["METRICS_FILE"] = file('metrics.txt') 
    args["ASSUME_SORT_ORDER"] = 'queryname'
    args["CREATE_INDEX"] = 'false'

    max = step(:BAM_bwa).info[:reads]
    args[:progress_bar] = gatk_read_count_monitor("BAM_duplicates", max) do |bar|
      ticks = bar.ticks
      ticks = max if ticks.to_i == 0
      set_info :reads, ticks
    end

    gatk("MarkDuplicates", args)

    Open.mv output, self.path
    nil
  end

  input :skip_duplicates, :boolean, "Skip MarkDuplicates", false
  input :type_of_sequencing, :select, "Whole genome or whole exome", nil, :select_options => %w(WGS WES panel)
  dep :BAM_duplicates do |jobname,options|
    if options[:skip_duplicates] || (options[:type_of_sequencing].to_s == "panel" && ! options[:skip_duplicates] == false)
      task = :BAM_bwa
    else
      task = :BAM_duplicates
    end
    {:task => task, :jobname => jobname, :inputs => options}
  end
  extension :bam
  task :BAM_sorted => :binary do 
    split_sort = config :split_sort, :bam, :gatk, :GATK, :default => false
    samtools_sort = config :samtools_sort, :bam_sort, :bam, :sort, :GATK, :default => false
    dep = dependencies.first

    job = if dep.info[:spark]
            dep
          else
            job = if samtools_sort
                    HTS.job(:sort_BAM_samtools, self.clean_name, :BAM => dep)
                  elsif split_sort
                    HTS.job(:sort_BAM_split, self.clean_name, :BAM => dep)
                  else
                    HTS.job(:sort_BAM, self.clean_name, :BAM => dep)
                  end
            self.dependencies = self.dependencies + [job]
            job
          end

    job.produce
    Open.link job.path, self.tmp_path
    nil
  end

  dep :BAM_sorted
  extension :bam
  input :interval_list, :file, "Interval list", nil, :nofile => true, :noload => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :known_sites, :array, "List of population resources for this reference"
  task :BAM_rescore => :binary do |interval_list,reference,known_sites|

    interval_list = nil if interval_list == "none"

    args = {}

    if known_sites.nil?
      known_sites = [] 
      known_site_codes = if reference.include?('mm10') || reference.include?('GRCm38')
                           ["mm10_variation", "mm10_structural"]
                         elsif reference.include?('rn6') || reference.include?('Rnor_6.0')
                           ["rn6_variation"]
                         else
                           ["miller_indels", "dbsnp", "1000g_snps"]
                         end
      known_site_codes.each do |file|
        vcf = vcf_file reference, file
        next if vcf.nil?
        vcf = GATK.prepare_VCF vcf 
        known_sites << vcf
      end
    end

    args["known-sites"] = known_sites

    orig_reference = reference_file(reference)
    reference = HTS.prepare_FASTA orig_reference
    #reference = BWA.prepare_FASTA orig_reference
    #reference = GATK.prepare_FASTA orig_reference
    #reference = Samtools.prepare_FASTA orig_reference

    args["reference"] = reference
    args["intervals"] = interval_list if interval_list
    args["output"] = file('recal_data.table')

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    input_bam_job = step(:BAM_sorted)

    #{{{ Recalibration
    shard = config('shard', :BaseRecalibrator, :baserecalibrator, :rescore, :gatk)

    if shard.to_s == 'true'
      break_interval = Rbbt::Config.get('break_interval', 'shard', 'GATK', 'gatk', 'interval', :default => false) if break_interval.nil?
      break_interval = false if break_interval.to_s.downcase == 'false'

      contigs = Samtools.bam_contigs(input_bam_job)
      bam_file = Samtools.prepare_BAM(input_bam_job)
      args["input"] = bam_file

      cpus = config('cpus', :BaseRecalibrator, :baserecalibrator, :rescore, :shard, :gatk)
      args["intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing BaseRecalibrator sharded")

      outfiles = Path.setup(file('outfiles_recall'))
      GATKShard.cmd("BaseRecalibrator", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar, break_interval) do |ioutfile|
        bar.tick
        Open.mv ioutfile, outfiles[File.basename(ioutfile + '.report')]
        nil
      end
      bar.remove 

      args = {}
      args["I"] = outfiles.glob("*.report")
      args["O"] = file('recal_data.table')
      gatk("GatherBQSRReports", args)
      Open.rm_rf outfiles
      nil
    else
      bam_file = interval_list ? Samtools.prepare_BAM(input_bam_job) : input_bam_job.path

      args["input"] = bam_file
      gatk("BaseRecalibrator", args)
    end

    #{{{ Apply
    output = file('out.bam')
    args = {}
    args["input"] = bam_file
    args["output"] = output
    args["bqsr-recal-file"] = file('recal_data.table')

    shard = config('shard', :gatk, :rescore, :apply_rescore, :apply_bqsr, :ApplyBQSR)

    if shard.to_s == 'true'
      contigs = Samtools.bam_contigs(input_bam_job)
      args["input"] = bam_file

      cpus = config('cpus', :ApplyBQSR, :apply_bqsr, :apply_rescore, :shard, :gatk)
      #args["intervals"] ||= nil
      #args["interval-padding"] ||= GATKShard::GAP_SIZE 
      #intervals = (interval_list || intervals_for_reference(reference))
      
      args["intervals"] = nil
      args["interval-padding"] = GATKShard::GAP_SIZE 
      intervals = intervals_for_reference(reference)

      outfiles = Path.setup(file('outfiles'))
      Open.mkdir outfiles

      bar = self.progress_bar("Processing ApplyBQSR sharded")
      GATKShard.cmd("ApplyBQSR", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar, false) do |ioutfile|
        bar.tick
        chr, pos = Samtools.BAM_start ioutfile
        target = outfiles[File.basename(ioutfile) + '.bam']
        Open.mv ioutfile, target if ! (chr.nil? or chr.empty?)
        nil
      end
      bar.remove 

      sorted_parts = outfiles.glob("*.bam").sort{|a,b| Misc.genomic_location_cmp_contigs(File.basename(a).split(",").first, File.basename(b).split(",").first, contigs, '__')}

      args = {}
      args["I"] = sorted_parts
      args["O"] = output
      args["CREATE_INDEX"] = 'false'
      args["CREATE_MD5_FILE"] = 'false'
      gatk("GatherBamFiles", args)

      Open.rm_rf outfiles
    else
      bam_file = interval_list ? Samtools.prepare_BAM(input_bam_job) : input_bam_job.path

      args["intervals"] = nil
      args["input"] = bam_file
      gatk("ApplyBQSR", args)
    end

    FileUtils.mv output, self.path
    nil
  end

  extension :bam
  input :skip_rescore, :boolean, "Skip BAM rescore", false
  dep_task :BAM, HTS, :BAM_rescore do |jobname,options|
    if options[:skip_rescore]
      task = :BAM_sorted
    else
      task = :BAM_rescore
    end
    {:task => task, :jobname => jobname, :inputs => options}
  end
end

