module HTS

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 hs37d5), :nofile => true
  extension :vcf
  task :BAM_pileup_sumaries_known_biallelic => :tsv do |reference|
    case reference.sub('_noalt','')
    when 'b37', 'hg19', 'hg38', 'GRCh38'
      variants_file = vcf_file reference, "1000g_snps"
    when 'mm10', 'GRCm38'
      raise RbbtException, "No Pileup Summaries for mouse"
      variants_file = vcf_file reference, "mm10_variation"
    end

    variants_file = GATK.prepare_VCF_AF_only variants_file

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

    variants_file = GATK.prepare_VCF step(:BAM_pileup_sumaries_known_biallelic).path, file('index')

    args = {}
    args["input"] = Samtools.prepare_BAM bam 
    args["variant"] = variants_file
    args["output"] = self.tmp_path
    args["intervals"] = interval_list ? interval_list : variants_file
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list
    begin
      GATK.run_log("GetPileupSummaries", args)
    rescue ProcessFailed
      raise RbbtException, "GetPileupSummaries failed"
    end
    nil
  end

  dep :BAM_pileup_sumaries
  input :matched, :file, "Matched contamiation table"
  task :contamination => :text do |matched|
    args = {}
    args["input"] = step(:BAM_pileup_sumaries).path
    args["output"] = self.tmp_path
    args["matched"] = matched if matched
    Open.mkdir files_dir
    args["segments"] = file('segments.tsv')
    GATK.run_log("CalculateContamination", args)
  end


  input :bam, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38), :nofile => true
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

  input :BAM, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38), :nofile => true
  extension "tar.gz"
  task :BAM_F1R2 => :binary do |bam,reference|
    args = {}

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    bam = Samtools.prepare_BAM(bam)

    args["--input"] = bam
    args["--reference"] = reference
    args["--output"] = self.tmp_path
    gatk("CollectF1R2Counts", args)
    nil
  end

  dep :BAM_F1R2
  extension "tar.gz"
  task :BAM_orientation_model => :text do
    args = {}

    args["--input"] = step(:BAM_F1R2).path
    args["--output"] = self.tmp_path
    gatk("LearnReadOrientationModel", args)
    nil
  end


  input :bam, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  task :BAM_summary => :text do |bam,reference|
    args = {}

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    bam = Samtools.prepare_BAM(bam)

    args["--INPUT"] = bam
    args["--REFERENCE_SEQUENCE"] = reference
    args["-O"] = self.tmp_path
    gatk("CollectAlignmentSummaryMetrics", args)
    nil
  end


  input :bam_1, :file, "BAM file", nil, :nofile => true
  input :bam_2, :file, "BAM file", nil, :nofile => true
  task :compare_BAM => :tsv do |bam_1,bam_2|
    reads_1 = file('reads_1')
    Open.write(reads_1, %w(read chr pos qual) * "\t" + "\n")
    CMD.cmd("samtools view #{bam_1} | cut -f 1,2,3,4,11 | sed 's/\\t/:/' >> #{reads_1}")

    reads_2 = file('reads_2')
    Open.write(reads_2, %w(read chr pos qual) * "\t" + "\n")
    CMD.cmd("samtools view #{bam_2} | cut -f 1,2,3,4,11 | sed 's/\\t/:/' >> #{reads_2}")
    
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
    if interval_list and interval_list.include? '.bed'
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
        CMD.cmd_log("qualimap bamqc -bam '#{bam}' --java-mem-size=4G -gff '#{fint}' -outdir '#{outdir}' ")
      end
    else
      CMD.cmd_log("qualimap bamqc -bam '#{bam}' --java-mem-size=4G -outdir '#{outdir}' ")
    end

    Open.cp File.join(outdir, "genome_results.txt"), self.tmp_path
    nil
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :positions, :array, "Genomic position"
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38), :nofile => true
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
end
