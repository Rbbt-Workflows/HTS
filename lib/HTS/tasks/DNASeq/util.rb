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
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :BAM_pileup_sumaries => :text do |bam,reference|

    vcf = vcf_file(reference, "small_exac")

    raise ParameterException, "No population VCF for pileup BAM pileup summaries" if vcf.nil?

    variants_file = GATK.prepare_VCF vcf
    args = {}
    args["input"] = Samtools.prepare_BAM bam 
    args["variant"] = variants_file
    args["output"] = self.tmp_path
    args["intervals"] = variants_file

    begin
      GATK.run_log("GetPileupSummaries", args)
    rescue ProcessFailed
      raise RbbtException, "GetPileupSummaries failed"
    end
    nil
  end

  dep :BAM_pileup_sumaries do |jobname,options|
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
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  extension 'pileup.gz'
  task :pileup => :text do |bam,reference|
    orig_reference = reference_file(reference)
    reference = Samtools.prepare_FASTA orig_reference

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


  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
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

    iif first
    iif last
    iif common.select{|mutation,parts| parts.select{|p| p.split("|").uniq.length != 1}.any? }

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
        CMD.cmd_log(:qualimap, "bamqc -bam '#{bam}' --java-mem-size=8G -gff '#{fint}' -outdir '#{outdir}' ")
      end
    else
      CMD.cmd_log(:qualimap, "bamqc -bam '#{bam}' --java-mem-size=8G -outdir '#{outdir}' ")
    end

    Open.cp File.join(outdir, "genome_results.txt"), self.tmp_path
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :positions, :array, "Genomic position"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :BAM_position_pileup => :tsv do |bam,positions,reference|
    bam = Samtools.prepare_BAM(bam)
    reference = reference_file reference
    reference = GATK.prepare_FASTA reference
    reference = Samtools.prepare_FASTA reference

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
  input :min_cov, :integer, "Min. coverage", 5
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
    out = ""
    TSV.traverse step(:genomecov), :into => out, :type => :array do |line|
      chr, start, eend, count = line.split("\t")
      last_chr = chr if last_chr.nil?
      if last_chr != chr || last_eend != start
        res = [last_chr, last_start, last_eend] * "\t"
        last_start = start
      else
        res = nil
      end
      last_chr = chr 
      last_eend = eend
      next unless res
      res  + "\n"
    end

    out << [last_chr, last_start, last_eend] * "\t"
    out
  end

end
