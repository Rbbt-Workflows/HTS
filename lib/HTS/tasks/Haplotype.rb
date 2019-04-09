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
    args["-ip"] = 100 if interval_list
    args["-ERC"] = "GVCF"
    args["--max-alternate-alleles"] = 3

    shard = Rbbt::Config.get('shard', :gatk, :haplotype, :HaplotypeCaller)

    if shard == 'true'
      headervcf = file('tmp.header')
      contentvcf = file('tmp.content')
      args["--intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing HaplotypeCaller sharded")

      bar.init
      GATKShard.cmd("HaplotypeCaller", args, intervals, 30_000_000) do |ioutfile|
        bar.tick
        `grep "#" "#{ioutfile}" > "#{headervcf}"` unless File.exists? headervcf
        `grep -v "#" #{ioutfile} >> #{contentvcf}` 
      end
      bar.done

      `cat "#{headervcf}" > "#{output}"`
      `sort -k 1,2 "#{contentvcf}" | uniq >> "#{output}"`

      FileUtils.rm Path.setup(files_dir).glob("tmp.*")
      nil
    else
      GATK.run_log("HaplotypeCaller", args)
    end

    Open.cp output, self.path
    nil
  end
end
