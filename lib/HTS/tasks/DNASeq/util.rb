require 'tools/Qualimap'

module HTS

  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 hs37d5), :nofile => true
  #extension :vcf
  #task :BAM_pileup_sumaries_known_biallelic => :tsv do |reference|
  #  case reference.sub('_noalt','')
  #  when 'b37', 'hg19', 'hg38', 'GRCh38'
  #    variants_file = vcf_file reference, "small_exac"
  #  when 'mm10', 'GRCm38'
  #    raise RbbtException, "No Pileup Summaries for mouse"
  #    variants_file = vcf_file reference, "mm10_variation"
  #  end

  #  variants_file = GATK.prepare_VCF variants_file

  #  reference = reference_file self.recursive_inputs[:reference]
  #  reference = GATK.prepare_FASTA reference

  #  args = {}
  #  args["reference"] = reference
  #  args["variant"] = variants_file
  #  args["restrict-alleles-to"] = 'BIALLELIC'
  #  args["output"] = tmp_path
  #  GATK.run_log("SelectVariants", args)
  #  nil
  #end

  #dep :BAM_pileup_sumaries_known_biallelic, :jobname => "Default"
  
  input :bam_file, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10)
  input :pileup_germline_resource, :file, "Germline resource for BAM_pileup_sumaries", :small_exac, :nofile => true
  task :BAM_pileup_sumaries => :text do |bam,reference,pileup_germline_resource|

    vcf = vcf_file(reference, pileup_germline_resource)

    if reference.include?('mm10') || reference.include?('GRCm38')
      Log.warn "Germline resource #{pileup_germline_resource} not found for #{reference} using mm10_variation instead"
      vcf = vcf_file(reference, "mm10_variation")
    elsif reference.include?('rn6') || reference.include?('Rnor_6.0')
      Log.warn "Germline resource #{pileup_germline_resource} not found for #{reference} using rn6_variation instead"
      vcf = vcf_file(reference, "rn6_variation")
    end if vcf.nil?

    raise ParameterException, "No population VCF for pileup BAM pileup summaries" if vcf.nil?

    orig_reference = reference_file(reference)
    reference = HTS.prepare_FASTA orig_reference

    variants_file = GATK.prepare_VCF vcf
    args = {}
    args["input"] = Samtools.prepare_BAM bam 
    args["variant"] = variants_file
    args["reference"] = reference
    args["output"] = self.tmp_path
    args["intervals"] = variants_file

    begin
      GATK.run_log("GetPileupSummaries", args)
    rescue ProcessFailed
      raise RbbtException, "GetPileupSummaries failed"
    end
    nil
  end

  dep :BAM_pileup_sumaries, :bam_file => :placeholder do |jobname,options|
    jobs = []
    jobs << {:inputs => options.merge(:bam_file => options[:tumor_bam]), :jobname => jobname}
    jobs << {:inputs => options.merge(:bam_file => options[:normal_bam]), :jobname => jobname + "_normal"} if options[:normal_bam]
    jobs
  end
  input :tumor_bam, :file, "Main BAM", nil, :nofile => true
  input :normal_bam, :file, "Matched BAM", nil, :nofile => true
  task :contamination => :text do 
    args = {}

    tumor_pileup = dependencies[0]
    matched = dependencies[1]

    args["input"] = tumor_pileup.path
    args["output"] = self.tmp_path

    args["matched"] = matched.path if matched

    Open.mkdir files_dir
    args["segments"] = file('segments.tsv')
    GATK.run_log("CalculateContamination", args)
  end


  input :bam, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10)
  extension 'pileup.gz'
  task :pileup => :text do |bam,reference|
    orig_reference = reference_file(reference)
    reference = HTS.prepare_FASTA orig_reference
    bam = bam.find if Path === bam

    monitor_cmd_genome "samtools mpileup -f '#{reference}' -Q 20 '#{bam}'", false, true
  end


  #input :tumor, :file, "Tumor BAM", nil, :nofile => true
  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  #task :BAM_artifact_metrics => :text do |tumor,reference|
  #  
  #  reference = reference_file reference
  #  orig_reference = reference

  #  reference = GATK.prepare_FASTA orig_reference
  #  reference = Samtools.prepare_FASTA orig_reference

  #  tumor = tumor.path if Step === tumor

  #  tumor = Samtools.prepare_BAM(tumor)

  #  FileUtils.mkdir_p files_dir
  #  args = {}

  #  args["-I"] = tumor
  #  args["-O"] = file('output')
  #  args["-R"] = reference
  #  args["--FILE_EXTENSION"] = '.txt'
  #  gatk("CollectSequencingArtifactMetrics", args)
  #  FileUtils.cp file('output' + '.pre_adapter_detail_metrics.txt'), self.path
  #  nil
  #end

  input :bam, :file, "Tumor BAM", nil, :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_coverage => :text do |bam,interval_list|
    args = {}

    bam = Samtools.prepare_BAM(bam)

    args["--input"] = bam
    args["--intervals"] = interval_list if interval_list
    args["--interval-merging-rule"] = "OVERLAPPING_ONLY"
    args["--format"] = "TSV"
    args["-O"] = self.tmp_path
    gatk("CollectReadCounts", args)
    nil
  end

  #input :BAM, :file, "Tumor BAM", nil, :nofile => true
  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38), :nofile => true
  #extension "tar.gz"
  #task :BAM_F1R2 => :binary do |bam,reference|
  #  args = {}

  #  reference = reference_file reference
  #  orig_reference = reference

  #  reference = GATK.prepare_FASTA orig_reference
  #  reference = Samtools.prepare_FASTA orig_reference
  #  bam = Samtools.prepare_BAM(bam)

  #  args["--input"] = bam
  #  args["--reference"] = reference
  #  args["--output"] = self.tmp_path
  #  gatk("CollectF1R2Counts", args)
  #  nil
  #end

  dep :mutect2_pre
  extension "tar.gz"
  task :BAM_orientation_model => :text do
    args = {}

    f1r2_path = step(:mutect2_pre).file('f1r2.tar.gz')
    if File.directory?(f1r2_path)
      args["--input"] = f1r2_path.glob("*.tar.gz")
    else
      args["--input"] = f1r2_path
    end
    args["--output"] = self.tmp_path
    gatk("LearnReadOrientationModel", args)
    nil
  end


  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :bam_1, :file, "BAM file", nil, :nofile => true
  input :bam_2, :file, "BAM file", nil, :nofile => true
  task :compare_BAM => :tsv do |bam_1,bam_2|
    reads_1 = file('reads_1')
    Open.write(reads_1, %w(read chr pos qual) * "\t" + "\n")
    CMD.cmd("samtools view --no-PG #{bam_1} | cut -f 1,2,3,4,11 | sed 's/\\t/:/' >> #{reads_1}")

    reads_2 = file('reads_2')
    Open.write(reads_2, %w(read chr pos qual) * "\t" + "\n")
    CMD.cmd("samtools view --no-PG #{bam_2} | cut -f 1,2,3,4,11 | sed 's/\\t/:/' >> #{reads_2}")

    first = []
    last = []
    common = []
    TSV.traverse TSV.paste_streams([reads_1, reads_2], :header_hash => "", :same_fields => true, :sort => true), :type => :array, :bar => self.progress_bar("Comparing BAM files") do |line|
      next if line =~ /^read/
      next if line =~ /^#/
      read, chr, pos, qual, *rest = line.split("\t", -1)
      aln = [read, chr, pos, qual] * "_"
      case
      when chr[0] == "|"
        first << aln
      when chr[-1] == "|"
        last << aln
      else
        common << [aln, rest]
      end
    end

    tsv = TSV.setup({}, "Statistic~Value#:type=:single")
    tsv["Missing"] = first.length
    tsv["Extra"] = last.length
    tsv["Common"] = common.length
    tsv["Common but different"] = common.select{|mutation,parts| parts.select{|p| p.split("|").uniq.length != 1}.any?}.length
    tsv
  end


  input :bam, :file, "BAM file", nil, :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  task :BAM_qualimap => :text do |bam,interval_list|
    bam = Samtools.prepare_BAM(bam)
    outdir = files_dir
    Open.mkdir outdir
    if interval_list && interval_list != :placeholder && interval_list.include?('.bed')
      TmpFile.with_file(:extension => 'bed') do |fint|
        Open.write(fint) do |f|
          TSV.traverse interval_list, :type => :array do |line|
            next if line =~ /^(@|browser|track)/
            chr, start, eend, name, score, strand, *rest = line.split("\t")
            parts = ["-"] * 6
            parts[0] = chr
            parts[1] = start
            parts[2] = eend
            parts[3] = name || "-"
            parts[4] = score || "1000"
            parts[5] = strand || "."
            f.puts parts * "\t"
          end
        end
        CMD.cmd_log(:qualimap, "bamqc -bam '#{bam}' --java-mem-size=8G -gff '#{fint}' -outdir '#{outdir}' ", :xvfb => true)
      end
    else
      CMD.cmd_log(:qualimap, "bamqc -bam '#{bam}' --java-mem-size=8G -outdir '#{outdir}' ", :xvfb => true)
    end

    Open.cp File.join(outdir, "genome_results.txt"), self.tmp_path
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :positions, :array, "Genomic position"
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :BAM_position_pileup => :tsv do |bam,positions,reference|
    bam = Samtools.prepare_BAM(bam)
    reference = reference_file reference
    reference = HTS.prepare_FASTA reference

    Open.mkdir files_dir
    output = file('output')
    TmpFile.with_file nil, true, :extension => 'bed' do |intervals|
      Open.write(intervals) do |file|
        positions.each do |position|
          chr, pos = position.split(":")
          file.puts [chr, pos.to_i - 1, pos.to_i] * "\t"
        end
      end

      args = {}
      args["input"] = bam
      args["reference"] = reference
      args["intervals"] = intervals
      args["output"] = output
      GATK.run_log("CollectAllelicCounts", args)
    end

    tsv = TSV.setup({}, "Genomic Position~Alt count,Ref count,Alt,Ref#:type=:list")
    TSV.traverse output, :type => :array, :into => tsv do |line|
      next if line =~ /^(?:@|CONTIG)/ 
      contig, pos, r_count, a_count, r, a = line.split("\t")

      position = [contig, pos] * ":"
      [position, [a_count, r_count, a, r]]
    end
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :sort_order, :select, "Sort order", :coordinate, :select_options => %w(coordinate queryname)
  extension :bam
  task :sort_BAM => :binary do |bam,sort_order|
    Open.mkdir files_dir 
    sorted = file('sorted.bam')

    args = {}
    args["INPUT"] = bam
    args["OUTPUT"] = sorted
    args["SORT_ORDER"] = sort_order
    gatk("SortSam", args)
    Open.mv sorted, self.path
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  extension :bam
  task :sort_BAM_samtools => :binary do |bam|
    Open.mkdir files_dir 
    sorted = file('sorted.bam')
    max_mem = config :max_mem, :samtools_sort_max_mem, :samtools_max_mem, :sort_samtools, :samtools, :sort
    cpus = config :cpus, :samtools_sort_cpus, :samtools_cpus, :sort_samtools, :samtools, :sort
    cpus ||= config :threads, :samtools_sort_threads, :samtools_threads, :sort_samtools, :samtools, :sort

    tmpdir ||= config :tmpdir, :samtools_sort_threads, :samtools_threads, :sort_samtools, :samtools, :sort
    if tmpdir 
      user = ENV["USER"] || `whoami`.strip
      tmpdir = tmpdir.sub("[USER]", user) 
    end
    tmpdir ||= files_dir
    Open.mkdir tmpdir
    CMD.cmd(:samtools, "sort '#{bam}' -O BAM -o '#{self.tmp_path}' -T #{tmpdir}", "-m" => max_mem, "--threads" => cpus)
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :reads, :integer, "Number of reads in chunk", 10_000_000
  input :sort_order, :select, "Sort order", :coordinate, :select_options => %w(coordinate queryname)
  extension :bam
  task :sort_BAM_split => :binary do |bam,reads,sort_order|
    Open.mkdir files_dir 

    header = file('header')
    CMD.cmd_log('samtools', "view -H --no-PG '#{bam}' > '#{header}'")
    txt = Open.read(header)

    tmp = file('tmp.sam')
    sorted = file('sorted.bam')
    chunks_dir = file('chunk')

    count = 0
    chunk = 1
    tmp_io = Open.open(tmp, :mode => 'w') 
    tmp_io.puts txt
    TSV.traverse CMD.cmd('samtools', "view --no-PG '#{bam}'", :pipe => true), :type => :array do |line|
      tmp_io.puts line
      count += 1

      if count == reads
        tmp_io.close

        args = {}
        args["INPUT"] = tmp
        args["OUTPUT"] = sorted
        args["SORT_ORDER"] = sort_order
        gatk("SortSam", args)

        Open.mv sorted, chunks_dir["chunk_#{chunk}.bam"]

        tmp_io = Open.open(tmp, :mode => 'w') 
        tmp_io.puts txt

        chunk += 1
        count = 0
      end
    end

    if count != 0
      tmp_io.close

      args = {}
      args["INPUT"] = tmp
      args["OUTPUT"] = sorted
      args["SORT_ORDER"] = sort_order
      gatk("SortSam", args)

      Open.mv sorted, chunks_dir["chunk_#{chunk}.bam"]
    end

    args = {}
    args["INPUT"] = chunks_dir.glob("*.bam").sort
    args["OUTPUT"] = self.tmp_path
    args["ASSUME_SORTED"] = true
    args["SORT_ORDER"] = sort_order
    gatk("MergeSamFiles", args)

    Open.rm_rf self.files_dir
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :min_cov, :integer, "Min. coverage", 20
  extension :bed
  task :genomecov => :text do |bam,min_cov|
    CMD.cmd(:bedtools, "genomecov -bg -max #{min_cov} -ibam '#{bam}'", :pipe => true)
  end

  dep :genomecov
  extension :bed
  task :intervals_from_BAM => :text do
    last_chr = nil
    last_start = nil
    last_eend = nil
    min = recursive_inputs[:min_cov].to_s
    Misc.open_pipe do |sin|
      begin
        io = TSV.traverse step(:genomecov), :into => :stream, :type => :array do |line|
          chr, start, eend, count = line.split("\t")
          next unless count == min
          if (last_chr && last_chr != chr) || (last_eend && last_eend != start)
            res = [last_chr, last_start, last_eend] * "\t" 
            last_start = start
          else
            res = nil
          end
          last_chr = chr 
          last_eend = eend
          last_start = start if last_start.nil?
          next unless res
          res  + "\n"
        end

        Misc.consume_stream(io, false, sin, false)

        sin << [last_chr, last_start, last_eend] * "\t" << "\n"
      rescue
        sin.stream_raise_exception $!
      end
    end
  end

  input :bam, :file, "BAM file", nil, :nofile => true
  input :bed_file, :file, "BED file", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension :bam
  task :extract_BAM_region_with_mates_samtools => :binary do |bam,bed_file,reference|
    cpus = config :cpus, :samtools_bam, :samtools_view, :samtools, :default => 2


    TmpFile.with_file do |reads|
      log :samtools, "Extracting read list with samtools"
      cpus = config :cpus, :samtools_view, :samtools, :view, :default => 1
      if reference
        orig_reference = reference_file(reference)
        reference = GATK.prepare_FASTA orig_reference
        CMD.cmd_log(:samtools, "view --threads #{cpus} -L '#{bed_file}' -T #{reference} '#{bam}' |cut -f 1 > '#{reads}'")
      else
        CMD.cmd_log(:samtools, "view --threads #{cpus} -L '#{bed_file}' '#{bam}' |cut -f 1 > '#{reads}'")
      end

      log :FilterSamReads, "Extracting reads with FilterSamReads"
      args = {}
      args["INPUT"] = bam
      args["OUTPUT"] = self.tmp_path
      args["READ_LIST_FILE"] = reads
      args["REFERENCE_SEQUENCE"] = reference if reference
      args["FILTER"] = 'includeReadList'

      gatk("FilterSamReads", args)

    end
    nil
  end

  input :bam, :file, "BAM file", nil, :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  extension :bam
  task :extract_BAM_region_with_mates => :binary do |bam,intervals|
    if intervals =~ /\.bed$/i
      fintervals = file('intervals.list')
      bed_to_intervals(intervals, fintervals, bam)
      intervals = fintervals
    end
    args = {}
    args["INPUT"] = bam
    args["OUTPUT"] = self.tmp_path
    args["INTERVAL_LIST"] = intervals
    args["FILTER"] = 'includePairedIntervals'

    gatk("FilterSamReads", args)
    nil
  end

  input :bam, :file, "BAM file", nil, :nofile => true
  extension :bam
  task :remove_spillover => :binary do |bam|
    io_in = CMD.cmd(:samtools, "view -h #{bam}", :pipe => true)

    last = nil
    io_out = TSV.traverse io_in, :type => :array, :into => :stream, :bar => true do |line|
      next line if line[0] == "@"
      pos = line.split("\t")[3].to_i
      next if last && pos <= last
      last = pos
      line
    end

    CMD.cmd(:samtools, "view -h -b - > #{self.tmp_path}", :in => io_out)
    nil
  end
end
