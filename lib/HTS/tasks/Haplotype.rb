module HTS
  input :BAM, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  extension :vcf
  task :haplotype => :text do |bam,reference,interval_list|

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    bam = bam.path if Step === bam
    
    bam = Samtools.prepare_BAM(bam)

    FileUtils.mkdir_p files_dir unless File.exists? files_dir
    output = file('calls.vcf')

    args = {}
    args["-I"] = bam
    args["--output"] = output
    args["-R"] = reference
    args["--intervals"] = interval_list if interval_list
    args["-ip"] = 500 if interval_list
    args["-ERC"] = "GVCF"
    args["--max-alternate-alleles"] = 3
    args["--create-output-variant-index"] = false

    shard = config('shard', :gatk, :haplotype, :HaplotypeCaller)

    if shard == 'true'
      contigs = Samtools.reference_contigs reference
      cpus = config('cpus', :shard, :haplotype, :HaplotypeCaller)
      headervcf = file('tmp.header')
      contentvcf = file('tmp.content')
      args["--intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing HaplotypeCaller sharded")

      outfiles = file('output')
      GATKShard.cmd("HaplotypeCaller", args, intervals, 10_000_000, cpus, contigs, bar) do |ioutfile|
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
