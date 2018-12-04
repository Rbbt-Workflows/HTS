module HTS

  
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension 'interval_list'
  task :preprocess_intervals => :tsv do |interval_list,reference|

    reference = reference_file reference
    reference = GATK.prepare_FASTA reference

    args = {}
    args["intervals"] =  interval_list
    args["reference"] =  reference
    args["bin-length"] =  0
    args["interval-merging-rule"] =  "OVERLAPPING_ONLY"
    args["output"] = self.tmp_path
    GATK.run_log("PreprocessIntervals", args)
    nil
  end

  dep :preprocess_intervals
  input :bam, :file, "BAM file", nil, :nofile => true
  extension 'hdf5'
  task :collect_fragment_counts => :binary do |bam|

    bam = Samtools.prepare_BAM(bam)
    args = {}
    args["input"] =  bam
    args["intervals"] =  step(:preprocess_intervals).path
    args["interval-merging-rule"] =  "OVERLAPPING_ONLY"
    args["output"] = self.tmp_path
    GATK.run_log("CollectFragmentCounts", args)
    nil
  end

  dep :preprocess_intervals
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension 'bed'
  task :germline_intervals => :tsv do |reference|
    file = GATK.known_sites[reference]["af-only-gnomad.vcf.gz"].produce.find

    TmpFile.with_file do |tmpfile|
      TmpFile.with_file nil, :extension => 'interval_list' do |tmp_interval_file|
        txt = step(:preprocess_intervals).path.read.split("\n").reject{|line| line =~ /^(Y|MT|GL|hs|NC)/} * "\n" + "\n"
        Open.write(tmp_interval_file, txt)
        args = {}
        args["variant"] =  file
        args["intervals"] =  tmp_interval_file
        args["output"] = tmpfile
        GATK.run_log("SelectVariants", args)
        TSV.traverse tmpfile, :type => :array, :into => :stream do |line|
          next if line =~ /^#/
          chr, pos, id, *rest = line.split("\t")
          [chr, pos, pos, "+", id] * "\t"
        end
      end
    end
  end


  dep :germline_intervals
  input :bam, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  extension 'tsv'
  task :collect_allelic_counts => :tsv do |bam,reference|

    reference = reference_file reference
    reference = GATK.prepare_FASTA reference

    bam = Samtools.prepare_BAM(bam)

    args = {}
    args["input"] =  bam
    args["reference"] =  reference
    args["intervals"] = step(:germline_intervals).path
    args["output"] = self.tmp_path
    GATK.run_log("CollectAllelicCounts ", args)
    nil
  end


end
