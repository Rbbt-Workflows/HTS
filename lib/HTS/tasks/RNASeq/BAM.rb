module HTS

  input :fastq1, :array, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :array, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil
  input :library_name, :string, "LIBRARY_NAME BAM field", "DefaultLibraryName"
  input :platform_unit, :string, "PLATFORM_UNIT BAM field", "DefaultPlatformUnit"
  input :platform, :string, "PLATFORM BAM field", "DefaultPlatform"
  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field", "DefaultSequencingCenter"
  input :read_length, :integer, "Read Length for index creation", 100
  extension :bam
  task :STAR => :binary do |fastq1, fastq2,organism,reference,read_group_name,sample_name,library_name,platform_unit,platform,sequencing_center,read_length|
    cpus = config("cpus", :STAR_build, :STAR, :default => 6) || 1
    samtools_cpus = config("cpus", :samtools_index, :samtools, :index, :default => nil)
    output = file('output')
    Open.mkdir output

    reference = 'b37' if reference.nil? && organism.nil?
    reference = 'hg38_noalt' if reference == 'hg38'
    organism = Organism.organism_for_build(reference) if organism.nil? and reference

    index = STAR.build_gft_index(organism, read_length, cpus)

    sample_name = clean_name if sample_name.nil?
    rg_str = '"ID:' << sample_name << '"'
    rg_str << ' "SM:' << sample_name << '"'
    rg_str << ' "LB:' << library_name << '"' if library_name
    rg_str << ' "PL:' << platform << '"' if platform
    rg_str << ' "PU:' << platform_unit << '"' if platform_unit
    rg_str << ' "CN:' << sequencing_center << '"' if sequencing_center

    fastqs_str = fastq2 ? "'#{fastq1*","}' '#{fastq2*","}'" : "'#{fastq1}'"
    CMD.cmd_log("STAR --runThreadN #{cpus} --genomeDir #{index} --readFilesIn #{fastqs_str} --readFilesCommand zcat --outFileNamePrefix #{output}/ --twopassMode Basic --outSAMattrRGline #{rg_str} --outSAMtype BAM SortedByCoordinate")
    #sam = output["Aligned.out.sam"]
    bam = output["Aligned.sortedByCoord.out.bam"]
    #CMD.cmd_log("samtools sort -@ #{samtools_cpus || 1} -O BAM -o #{bam} #{sam} ")
    #Open.rm sam
    Open.mv bam, self.tmp_path
    nil
  end

  input :fastq1, :array, "FASTQ 1 file", nil, :nofile => true
  input :fastq2, :array, "FASTQ 2 file", nil, :nofile => true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :phred, :select, "Phred Qualities", 'phred33', :select_options => %w(phred33 phred64)
  input :rna_strandness, :select, "RNA Strandness", 'FR', :select_options => %w(FR)
  input :read_group_name, :string, "READ_GROUP_NAME BAM field", nil
  input :sample_name, :string, "SAMPLE_NAME BAM field", nil
  input :library_name, :string, "LIBRARY_NAME BAM field", "DefaultLibraryName"
  input :platform_unit, :string, "PLATFORM_UNIT BAM field", "DefaultPlatformUnit"
  input :platform, :string, "PLATFORM BAM field", "DefaultPlatform"
  input :sequencing_center, :string, "SEQUENCING_CENTER BAM field", "DefaultSequencingCenter"
  extension :bam
  task :hisat => :binary do |fastq1, fastq2,organism,reference,phred,rna_strandness,read_group_name,sample_name,library_name,platform_unit,platform,sequencing_center|
    cpus = config("cpus", :hisat_build, :hisat, :default => 1)
    samtools_cpus = config("cpus", :samtools_index, :samtools, :index, :default => nil)

    reference = 'b37' if reference.nil? && organism.nil?
    organism = Organism.organism_for_build(reference) if organism.nil? and reference

    index = HISAT.build_gft_index(organism, cpus)

    set_info :organism, organism

    splice_in = Path.setup(index).replace_extension('ss')

    Open.mkdir files_dir
    sam = file('out.sam')
    splice_out = file('novel_splicesites.ss')
    summary = file('alignment_summary.txt')

    sample_name = clean_name if sample_name.nil?
    rg_str = "--rg-id " << sample_name
    rg_str << " --rg SM:" << sample_name 
    rg_str << " --rg LB:" << library_name if library_name
    rg_str << " --rg PL:" << platform if platform
    rg_str << " --rg PU:" << platform_unit if platform_unit
    rg_str << " --rg CN:" << sequencing_center if sequencing_center

    Open.mkdir files_dir
    CMD.cmd_log("hisat2 -p #{cpus || 1} --dta --#{phred} #{rg_str} --summary-file #{summary} --known-splicesite-infile #{splice_in} --novel-splicesite-outfile #{splice_out} -x #{index} -1 #{fastq1*","} -2 #{fastq2*""} -S #{sam}")
    CMD.cmd_log("samtools sort --no-PG -@ #{samtools_cpus || 1} -O BAM -o #{self.tmp_path} #{sam} ")
    Open.rm sam
    nil
  end

  input :aligner, :select, "RNA aligner", :STAR, :select_options => %w(hisat STAR)
  dep :STAR do |jobname,options|
    if options[:aligner].to_s == "STAR"
      {:jobname => jobname, :inputs => options}
    end
  end
  dep :hisat do |jobname,options|
    if options[:aligner].to_s == "hisat"
      {:jobname => jobname, :inputs => options}
    end
  end
  extension :bam
  task :RNA_BAM_duplicates => :binary do
    Open.mkdir files_dir 
    output = file('out.bam')

    args = {}
    args["INPUT"] = dependencies.first.path
    args["OUTPUT"] = output
    args["METRICS_FILE"] = file('metrics.txt') 
    args["ASSUME_SORT_ORDER"] = 'queryname'
    args["CREATE_INDEX"] = 'false'

    gatk("MarkDuplicates", args)

    Open.mv output, self.path
    nil
  end

  dep :RNA_BAM_duplicates
  extension :bam
  task :RNA_BAM_sorted => :binary do
    Open.mkdir files_dir 
    sorted = file('sorted.bam')
    
    args = {}
    args["INPUT"] = step(:RNA_BAM_duplicates).path
    args["OUTPUT"] = sorted
    args["SORT_ORDER"] = 'coordinate'
    gatk("SortSam", args)
    Open.mv sorted, self.path
    nil
  end

  dep :RNA_BAM_sorted
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :RNA_BAM_cigars => :binary do |organism,reference|
    Open.mkdir files_dir 
    fixed = file('fixed.bam')

    if reference.nil? && organism
      reference = Organism.hg_build(organism) 
      reference = 'b37' if reference == 'hg19'
    end
    
    orig_reference = reference_file(reference)
    reference = BWA.prepare_FASTA orig_reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    args = {}
    args["input"] = step(:RNA_BAM_sorted).path
    args["output"] = fixed
    args["reference"] = reference
    gatk("SplitNCigarReads", args)
    Open.mv fixed, self.path
    nil
  end

  dep :RNA_BAM_cigars
  extension :bam
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :RNA_BAM_rescore => :binary do |interval_list|

    interval_list = nil if interval_list == "none"

    reference = reference_file self.recursive_inputs[:reference]
    reference = GATK.prepare_FASTA reference

    args = {}
    args["reference"] = reference
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list
    args["output"] = file('recal_data.table')

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
      vcf = GATK.prepare_VCF_AF_only vcf 
      known_sites << vcf
    end

    args["known-sites"] = known_sites

    FileUtils.mkdir_p files_dir unless Open.exists?(files_dir)

    #{{{ Recalibration
    shard = config('shard', :gatk, :rescore, :baserecalibrator, :BaseRecalibrator)
    if shard == 'true'
      contigs = Samtools.bam_contigs(step(:RNA_BAM_cigars))
      bam_file = Samtools.prepare_BAM(step(:RNA_BAM_cigars))
      args["input"] = bam_file

      cpus = config('cpus', :shard, :rescore, :baserecalibrator, :BaseRecalibrator)
      args["intervals"] ||= nil
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing BaseRecalibrator sharded")

      outfiles = Path.setup(file('outfiles_recall'))
      GATKShard.cmd("BaseRecalibrator", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
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
      bam_file = interval_list ? Samtools.prepare_BAM(step(:RNA_BAM_cigars)) : step(:RNA_BAM_cigars).path

      args["input"] = bam_file
      gatk("BaseRecalibrator", args)
    end

    #{{{ Apply
    output = file('out.bam')
    args = {}
    args["input"] = bam_file
    args["output"] = output
    args["bqsr-recal-file"] = file('recal_data.table')
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list

    shard = config('shard', :gatk, :rescore, :apply_rescore, :apply_bqsr, :ApplyBQSR)

    if shard == 'true'
      contigs = Samtools.bam_contigs(step(:RNA_BAM_cigars))
      args["input"] = bam_file

      cpus = config('cpus', :shard, :rescore, :apply_rescore, :apply_bqsr, :ApplyBQSR)
      args["intervals"] ||= nil
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      intervals = (interval_list || intervals_for_reference(reference))

      outfiles = Path.setup(file('outfiles'))
      Open.mkdir outfiles

      bar = self.progress_bar("Processing ApplyBQSR sharded")
      GATKShard.cmd("ApplyBQSR", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
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
      bam_file = interval_list ? Samtools.prepare_BAM(step(:RNA_BAM_cigars)) : step(:RNA_BAM_cigars).path

      args["input"] = bam_file
      gatk("ApplyBQSR", args)
    end

    FileUtils.mv output, self.path
    nil
  end

  extension :bam
  dep_task :RNA_BAM, HTS, :RNA_BAM_rescore
end
