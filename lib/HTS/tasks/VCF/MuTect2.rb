require 'tools/shard'
module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :pon, :file, "Panel of normals", nil, :nofile => true
  input :germline_resource, :file, "Germline resource", :gnomad, :nofile => true
  input :af_not_in_resource, :float, "Allele frequency of alleles not in resource", nil
  extension :vcf
  task :mutect2_pre => :text do |tumor,normal,reference,interval_list,pon,germline_resource,af_not_in_resource|

    interval_list = nil if interval_list == "none"

    af_not_in_resource = germline_min_af germline_resource if af_not_in_resource.nil? and germline_resource
    germline_resource = vcf_file reference, germline_resource if germline_resource
    germline_resource = GATK.prepare_VCF_AF_only germline_resource if germline_resource

    pon = GATK.prepare_VCF_AF_only pon if pon

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal) if normal

    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    output = file('calls.vcf')

    args["input"] = [tumor, normal].compact
    args["output"] = output
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample if normal_sample
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = GATKShard::GAP_SIZE if interval_list
    args["panel-of-normals"] = pon if pon
    args["bam-output"] = file('haplotype.bam')
    args["germline-resource"] = germline_resource

    # UPDATE FOR GATK 4.1.2
    #args["af-of-alleles-not-in-resource"] = "%.10f" % af_not_in_resource.to_s if af_not_in_resource

    shard = config('shard', :gatk, :mutect, :mutect2)

    if shard == 'true'
      contigs = Samtools.reference_contigs reference
      cpus = config('cpus', :shard, :mutect, :mutect2)
      headervcf = file('tmp.header')
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      contentvcf = file('tmp.content')
      headervcf_stats = file('tmp.header.stats')
      contentvcf_stats = file('tmp.content.stats')
      tmp_stats = file('tmp.stats')
      args["intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing Mutect2 sharded")

      GATKShard.cmd("Mutect2", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
        bar.tick
        `grep "#" "#{ioutfile}" > "#{headervcf}"` unless File.exists? headervcf
        `grep -v "#" #{ioutfile} >> #{contentvcf}` 
        `head -n 1 "#{ioutfile}.stats" > "#{headervcf_stats}"` unless File.exists? headervcf_stats
        `tail -n 1 #{ioutfile}.stats >> #{contentvcf_stats}` 
      end
      bar.remove 

      `cat "#{headervcf}" > "#{output}"`
      `sort -k 1,2 "#{contentvcf}" | uniq >> "#{output}"`

      `cat "#{headervcf_stats}" > "#{tmp_stats}"`
      `cat "#{contentvcf_stats}" >> "#{tmp_stats}"`
      stats = TSV.open( tmp_stats, :type => :double, :header_hash => '', :merge => true)
      stats = stats.to_list{|l| Misc.sum(l.collect{|e| (e.nil? || e.empty?) ? nil : e.to_f }.compact) }
      Open.write(output + '.stats') do |f|
        f.puts ["statistic", "value"] * "\t"
        f.puts ["callable", stats["callable"]] * "\t"
      end
      FileUtils.rm Path.setup(files_dir).glob("tmp.*")
      nil
    else
      gatk("Mutect2", args)
    end

    Open.cp output, self.path
    nil
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  dep :mutect2_pre
  dep :contamination, :BAM => :normal, :compute => true do |jobname,options,dependencies|
    if options[:normal]
      options[:BAM] = options[:normal]
      {:inputs => options}
    end
  end
  dep :contamination, :BAM => :tumor do |jobname,options,dependencies|
    matched_dep = dependencies.flatten.select{|dep| dep.task_name.to_sym == :contamination}.first
    options[:matched] = matched_dep.step(:BAM_pileup_sumaries).path if matched_dep
    options[:BAM] = options[:tumor]
    {:inputs => options}
  end
  dep :BAM_orientation_model, :BAM => :tumor
  extension :vcf
  task :mutect2_filtered => :text do |reference|
    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2_pre).path, tmp
    FileUtils.ln_s step(:mutect2_pre).file('calls.vcf.stats'), tmp + '.stats'

    FileUtils.mkdir_p files_dir
    args = {}

    contamination = dependencies.select{|dep| dep.task_name.to_sym == :contamination}.last

    args["variant"] = tmp
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["--orientation-bias-artifact-priors"] = step(:BAM_orientation_model).path
    if step(:contamination).path.read.include? "NaN"
      set_info :missing_contamination, true
      Log.warn "NaN in contamination file: #{Log.color :blue, self.path}"
    else
      args["tumor-segmentation"] = contamination.file('segments.tsv')
      args["contamination-table"] = contamination.path
    end
    gatk("FilterMutectCalls", args)
    nil
  end

  dep :mutect2_filtered
  extension :vcf
  task :mutect2_clean => :text do
    TSV.traverse step(:mutect2_filtered), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end

  #dep :BAM_artifact_metrics, :compute => :bootstrap
  #dep :mutect2_clean, :compute => :bootstrap
  #extension :vcf
  #task :mutect2_orientation_bias => :text do
  #  FileUtils.mkdir_p files_dir
  #  args = {}

  #  args["-P"] = step(:BAM_artifact_metrics).path
  #  args["-V"] = step(:mutect2_clean).join.path
  #  args["-O"] = self.tmp_path
  #  gatk("FilterByOrientationBias", args)
  #  nil
  #end

  extension :vcf
  dep_task :mutect2, HTS, :mutect2_clean

end
