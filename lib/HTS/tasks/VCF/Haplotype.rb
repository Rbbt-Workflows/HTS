module HTS
  input :BAM, :file, "BAM file", nil, :nofile => true
  input :BAM_list, :array, "List of BAM files"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :reference_confidence_mode, :select, "Mode for emitting reference confidence scores", "NONE", :select_options => %w(NONE BP_RESOLUTION GVCF)
  extension :vcf
  task :haplotype => :text do |bam,bam_list,reference,interval_list,reference_confidence_mode|

    interval_list = nil if interval_list == "none"

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    prepared_bams = []

    if bam_list.nil?
      bam = bam.path if Step === bam
      bam = Samtools.prepare_BAM(bam)
    else
      bam_list.each{ |file| prepared_bams.push(Samtools.prepare_BAM(file))}
      bam = prepared_bams
    end

    FileUtils.mkdir_p files_dir unless File.exists? files_dir
    output = file('calls.vcf')

    args = {}
    args["-I"] = bam
    args["--output"] = output
    args["-R"] = reference
    args["--intervals"] = interval_list if interval_list
    args["-ip"] = 500 if interval_list
    args["-ERC"] = reference_confidence_mode
    args["--max-alternate-alleles"] = 3
    args["--create-output-variant-index"] = false
    shard = config('shard', :HaplotypeCaller, :haplotype,  :gatk)

    if shard == 'true'
      contigs = Samtools.reference_contigs reference
      cpus = config('cpus', :HaplotypeCaller, :haplotype,  :gatk, :default => 2)
      headervcf = file('tmp.header')
      contentvcf = file('tmp.content')
      args["--intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing HaplotypeCaller sharded")

      outfiles = file('output')
      GATKShard.cmd("HaplotypeCaller", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
        bar.tick
        Open.mv ioutfile, outfiles[File.basename(ioutfile) + '.vcf']
      end

      contigs = Samtools.reference_contigs reference
      sorted_parts = outfiles.glob("*.vcf").sort{|a,b| Misc.genomic_location_cmp_contigs(File.basename(a), File.basename(b), contigs, '__')}

      args = {}
      args["INPUT"] = sorted_parts
      args["OUTPUT"] = output
      gatk("MergeVcfs", args)

      Open.rm_rf outfiles
      nil
    else
      gatk("HaplotypeCaller", args)
    end

    Open.mv output, self.path
    nil
  end
end
