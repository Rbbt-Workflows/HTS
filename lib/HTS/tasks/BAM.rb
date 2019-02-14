require 'rbbt/sources/organism'
module HTS
  input :fastq1, :file, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :file, "FASTQ 2 file", nil, :nofile => true
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
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

    GATK.run_log("FastqToSam", args)
  end

  dep :uBAM, :compute => :produce
  extension :ubam
  task :mark_adapters => :binary do
    args = {}
    args["INPUT"] = step(:uBAM).path
    args["METRICS"] = file('metrics.txt')
    args["OUTPUT"] = self.tmp_path

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkIlluminaAdapters", args)
  end

  dep :mark_adapters
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :bwa_mem_args, :string, "Arg string", "-M -p"
  extension :bam
  task :BAM => :binary do |reference, bwa_mem_args|

    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    Open.rm file('SamToFastq')
    Misc.with_fifo(file('SamToFastq')) do |s2f_path|

      args = {}
      args["INPUT"] = step(:mark_adapters).path
      args["FASTQ"] = s2f_path
      args["CLIPPING_ATTRIBUTE"] = "XT"
      args["CLIPPING_ACTION"] = "2"
      args["INTERLEAVE"] = "true"
      args["NON_PF"] = "true"

      io_s2f = GATK.run("SamToFastq", args)
      t_s2f = Thread.new do
        while line = io_s2f.gets
          Log.debug line
        end
      end

      bwa_mem_args += " -t " << config('cpus', 'bwa', :default => 8) 
      io_bwa = BWA.mem([s2f_path], reference, bwa_mem_args)

      args = {}
      args["ALIGNED_BAM"] = "/dev/stdin"
      args["UNMAPPED_BAM"] = step('uBAM').path
      args["OUTPUT"] = self.tmp_path
      args["R"] = reference
      args["CREATE_INDEX"] = "true"
      args["ADD_MATE_CIGAR"] = "true"
      args["CLIP_ADAPTERS"] = "false"
      args["CLIP_OVERLAPPING_READS"] = "true"
      args["INCLUDE_SECONDARY_ALIGNMENTS"] = "true"
      args["MAX_INSERTIONS_OR_DELETIONS"] = "-1"
      args["PRIMARY_ALIGNMENT_STRATEGY"] = "MostDistant"
      args["ATTRIBUTES_TO_RETAIN"] = "XS"
      GATK.run_log("MergeBamAlignment", args, io_bwa)

      FileUtils.rm_rf s2f_path
    end
    nil
  end

  dep :BAM
  extension :bam
  task :BAM_duplicates => :binary do
    args = {}
    args["INPUT"] = step(:BAM).path
    args["OUTPUT"] = self.tmp_path
    args["METRICS_FILE"] = file('metrics.txt')

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkDuplicates", args)
  end

  input :bam_files, :array, "BAM filenames to multiplex"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension :bam
  task :BAM_multiplex => :binary do |bam_filenames|
    args= {}
    bam_filenames = Dir.glob(File.join(bam_filenames.first, "*.bam")) if Array === bam_filenames && bam_filenames.length == 1 && File.directory?(bam_filenames.first)
    args["INPUT"] = bam_filenames
    args["OUTPUT"] = self.tmp_path
    args["METRICS_FILE"] = file('metrics.txt')

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("MarkDuplicates", args)
  end


  dep :BAM_duplicates
  extension :bam
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_rescore => :binary do |interval_list|

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference
    reference_code = self.recursive_inputs[:reference]

    DbSNP["All.vcf.gz.tbi"].produce.find
    db_SNP = DbSNP["All.vcf.gz"].produce.find

    bam_file = interval_list ? Samtools.prepare_BAM(step(:BAM_duplicates)) : step(:BAM_duplicates).path

    args = {}
    args["input"] = bam_file
    args["reference"] = reference
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    args["output"] = file('recal_data.table')

    known_sites = [] 
    ["dbsnp_138.vcf.gz","Miller_1000G_indels.vcf.gz", "1000G_phase1.indels.vcf.gz"].each do |file|
      known_sites << GATK.known_sites[reference_code][file].produce.find
    end
    args["known-sites"] = known_sites

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)
    GATK.run_log("BaseRecalibrator", args)

    args = {}
    args["input"] = bam_file
    args["output"] = self.tmp_path
    args["bqsr-recal-file"] = file('recal_data.table')
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    GATK.run_log("ApplyBQSR", args)
  end

  input :bam_file, :binary, "Bam file"
  task :revert_BAM => :binary do |bam_file|
    args = {}
    args["INPUT"] = bam_file
    args["OUTPUT"] = file("uBAM").find
    args["OUTPUT_BY_READGROUP"] = "true"
    
    args["SANITIZE"] = "true"
    args["MAX_DISCARD_FRACTION"] = "0.005"
    args["ATTRIBUTE_TO_CLEAR"] = "XA"
    args["ATTRIBUTE_TO_CLEAR"] = "BD"
    args["ATTRIBUTE_TO_CLEAR"] = "BI"
    args["SORT_ORDER"] = "queryname"
    args["RESTORE_ORIGINAL_QUALITIES"] = "true"
    args["REMOVE_DUPLICATE_INFORMATION"] = "true"
    args["REMOVE_ALIGNMENT_INFORMATION"] = "true"

    Open.mkdir file("uBAM").find
    GATK.run_log("RevertSam", args)
    file("uBAM").glob("*")
  end

  input :fastq1_files, :array, "FASTQ files for first mate"
  input :fastq2_files, :array, "FASTQ files for second mate", []
  input :uBAM_files, :array, "uBAM files for second mate", []
  dep :BAM, :compute => :produce do |jobname,options,uBAM_files|
    fastq1_files = options[:fastq1_files]
    if fastq1_files
      fastq2_files = options[:fastq2_files]
      fastq1_files.zip(fastq2_files).collect do |fastq1,fastq2|
        read_group_name = File.basename(fastq1).sub(/(_{1,2})?\.fastq.*/,'')
        options = options.merge({:fastq1 => fastq1, :fastq2 => fastq2, :read_group_name => read_group_name})
        {:inputs => options, :jobname => [jobname, read_group_name] * "." }
      end
    elsif uBAM_files
      uBAM_files = Dir.glob(File.join(uBAM_files, '*.bam')) + Dir.glob(File.join(uBAM_files, '*.ubam')) if String === uBAM_files && File.directory?(uBAM_files)
      [uBAM_files].flatten.collect do |uBAM|
        read_group_name = File.basename(uBAM).sub(/.u?bam/i,'')
        options = options.merge({"HTS#uBAM" => uBAM})
        {:inputs => options, :jobname => [jobname, read_group_name] * "." }
      end
    else
      raise "No FASTQ or uBAM files for #{ jobname }"
    end
  end
  dep :BAM_multiplex, :compute => :produce do |jobname, options,dependencies|
    bam_files = dependencies.flatten.collect{|dep| dep.path}
    {:jobname => jobname, :inputs => options.merge(:bam_files => bam_files)}
  end
  dep_task :BAM_rescore_mutiplex, HTS, :BAM_rescore do |jobname,options, dependencies|
    mutiplex = dependencies.flatten.select{|dep| dep.task_name == :BAM_multiplex}.first
    {:inputs => options.merge("HTS#BAM_duplicates" =>  mutiplex), :jobname => jobname}
  end




  dep :revert_BAM, :compute => :produce
  dep :BAM, :compute => :produce do |jobname, options, dependencies|
    read_groups = CMD.cmd("samtools view -H #{options[:bam_file]}").read.split("\n").select{|line| 
      line =~ /^@RG/
    }.collect{|line| 
      line.split("\t").select{|part| part =~ /^ID:(.*)/}.first.split(":").last
    }

    read_groups.collect do |read_group|
      uBAM = dependencies.first.file('uBAM')[read_group] + ".bam"
      {:task => :BAM, :inputs => options.merge({"HTS#uBAM" => uBAM}), :jobname => [jobname, read_group] * "."}
    end
  end
  dep :BAM_multiplex, :compute => :produce do |jobname, options,dependencies|
    bam_files = dependencies.flatten.select{|dep| dep.task_name == :BAM}.collect{|dep| dep.path}
    {:jobname => jobname, :inputs => options.merge(:bam_files => bam_files)}
  end
  dep_task :BAM_rescore_realign, HTS, :BAM_rescore do |jobname,options, dependencies|
    mutiplex = dependencies.flatten.select{|dep| dep.task_name == :BAM_multiplex}.first
    {:inputs => options.merge("HTS#BAM_duplicates" =>  mutiplex), :jobname => jobname}
  end





  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 hs37d5), :nofile => true
  extension :vcf
  task :BAM_pileup_sumaries_known_biallelic => :tsv do |reference|
    variants_file = case reference
                    when 'b37', 'hg19', 'hg38', 'hs37d5'
                      GATK.known_sites[reference]["1000G_phase1.snps.high_confidence.vcf.gz"].produce.find
                    else 
                      if m = Pathname.new(reference).realpath.to_s.match(/(b37|hg19|hg38|GRCh38|hs37d5)/)
                        code = m[1]
                        code = 'hg38' if code == 'GRCh38'
                        GATK.known_sites[code]["1000G_phase1.snps.high_confidence.vcf.gz"].produce.find
                      else
                        raise ParameterException.new "Cannot file a suitable variant file for reference: #{Misc.fingerprint reference}"
                      end
                    end

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference

    args = {}
    args["reference"] = reference
    args["variant"] = variants_file
    args["restrict-alleles-to"] = 'BIALLELIC'
    args["output"] = tmp_path
    GATK.run_log("SelectVariants", args)
    nil
  end

  dep :BAM_pileup_sumaries_known_biallelic, :jobname => "Default"
  input :BAM, :file, "BAM file", nil, :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_pileup_sumaries => :text do |bam,interval_list|

    variants_file = step(:BAM_pileup_sumaries_known_biallelic).path

    args = {}
    args["feature-file"] = variants_file
    GATK.run_log("IndexFeatureFile", args)

    args = {}
    args["input"] = Samtools.prepare_BAM bam 
    args["variant"] = variants_file
    args["output"] = self.tmp_path
    args["intervals"] = interval_list ? interval_list : variants_file
    args["interval-padding"] = 100 if interval_list
    GATK.run_log("GetPileupSummaries", args)
    nil
  end

  dep :BAM_pileup_sumaries
  task :contamination => :text do
    args = {}
    args["input"] = step(:BAM_pileup_sumaries).path
    args["output"] = self.tmp_path
    GATK.run_log("CalculateContamination", args)
  end


  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  extension :bam
  task :razers3_BAM => :text do |reference,fastq1,fastq2|

    reference = reference_file reference
    reference = BWA.prepare_FASTA reference
    cpus = Rbbt::Config.get(:cpus, :razers)
    TmpFile.with_file :extension => :bam do |tmp_bam|
      if fastq2
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}' '#{fastq2}'")
      else
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}'")
      end
      CMD.cmd("samtools sort '#{tmp_bam}' > '#{self.tmp_path}'")
    end
    nil
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  input :bowtie_args, :string, "Bowtie2 arguments", "--end-to-end"
  extension :bam
  task :bowtie_BAM => :text do |reference,fastq1,fastq2,bowtie_args|

    reference_index_dir = file('reference_index') if Step === reference

    reference = reference_file reference
    reference = Bowtie.prepare_FASTA reference, reference_index_dir
    cpus = Rbbt::Config.get(:cpus, :bowtie) || 1

    original_fastq1 = fastq1
    original_fastq2 = fastq2 if fastq2

    fastq1 = file('f1.fastq')
    fastq2 = file('f2.fastq') if fastq2

    fastq1 += '.gz' if original_fastq1 =~ /\.gz$/
    fastq2 += '.gz' if fastq2 && original_fastq2 =~ /\.gz$/

    Open.ln_s original_fastq1, fastq1
    Open.ln_s original_fastq2, fastq2 if fastq2

    TmpFile.with_file :extension => :bam do |tmp_bam|
      if fastq2
        CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -1 '#{fastq1}' -2 '#{fastq2}' -S '#{tmp_bam}'  ")
      else
        CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -U '#{fastq1}' -S '#{tmp_bam}'")
      end
      CMD.cmd("samtools sort '#{tmp_bam}' > '#{self.tmp_path}'")
    end
    nil
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  input :novoalign_args, :string, "NovoAlign arguments", "-F STDFQ -R 0 -r All 9999 -o SAM -o FullNW"
  extension "bam"
  task :novoalign_BAM => :binary do |reference,fastq1,fastq2,novoalign_args|
    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    reference = reference_file reference
    cpus = Rbbt::Config.get(:cpus, :novoalign)
    begin
      original_fastq1 = fastq1
      original_fastq2 = fastq2 if fastq2
      original_reference = reference

      reference = file('reference.fasta').find if original_reference =~ /\.gz/
      if original_reference != reference
        CMD.cmd("gunzip '#{original_reference}' -c > '#{reference}'")

        CMD.cmd("'#{Rbbt.software.opt.NovoAlign.novoindex.find}' '#{ reference }.nix' '#{reference}'")
      else
        reference = NovoAlign.prepare_FASTA reference 
      end

      fastq1 = file('file.1.fastq').find if original_fastq1 =~ /\.gz/
      fastq2 = file('file.2.fastq').find if fastq2 && original_fastq2 =~ /\.gz/

      CMD.cmd("gunzip '#{original_fastq1}' -c > '#{fastq1}'") if original_fastq1 != fastq1
      CMD.cmd("gunzip '#{original_fastq2}' -c > '#{fastq2}'") if fastq2 && original_fastq2 != fastq2


      TmpFile.with_file :extension => :sam do |tmp_sam|
        if fastq2
          CMD.cmd_log("#{Rbbt.software.opt.NovoAlign.novoalign.find} -d '#{reference}.nix' -f '#{fastq1}' '#{fastq2}' #{novoalign_args} 1> #{tmp_sam}")
        else
          CMD.cmd_log("#{Rbbt.software.opt.NovoAlign.novoalign.find} -d '#{reference}.nix' -f '#{fastq1}' #{novoalign_args} 1> #{tmp_sam}")
        end
        CMD.cmd("samtools sort -O BAM  '#{tmp_sam}' > '#{self.tmp_path}'")
      end
    ensure
      FileUtils.rm fastq1 if File.exists?(fastq1) && original_fastq1 != fastq1
      FileUtils.rm fastq2 if fastq2 && File.exists?(fastq2) && original_fastq2 != fastq2
      FileUtils.rm reference if File.exists?(reference) && original_reference != reference
      FileUtils.rm reference + '.nix' if File.exists?(reference + '.nix') && original_reference != reference
    end
    nil
  end

  input :bam, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension 'pileup.gz'
  task :pileup => :text do |bam,reference|
    orig_reference = reference_file(reference)
    reference = Samtools.prepare_FASTA orig_reference

    monitor_cmd_genome "samtools mpileup -f '#{reference}' -Q 20 '#{bam}'", false, true
  end

end
