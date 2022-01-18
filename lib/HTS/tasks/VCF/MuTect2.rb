require 'tools/shard'

module HTS

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :interval_padding, :integer, "Interval padding", 50
  input :pon, :file, "Panel of normals", 'default', :nofile => true
  input :germline_resource, :file, "Germline resource", :gnomad, :nofile => true
  input :af_not_in_resource, :float, "Allele frequency of alleles not in resource", nil
  input :remove_soft_clip, :boolean, "Don't consider soft clip bases", false
  input :max_mnp_distance, :integer, "Max distance for mnp merge", 1
  extension :vcf
  task :mutect2_pre => :text do |tumor,normal,reference,interval_list,interval_padding,pon,germline_resource,af_not_in_resource,remove_soft_clip,max_mnp_distance|

    interval_list = nil if interval_list == "none"

    af_not_in_resource = germline_min_af germline_resource if af_not_in_resource.nil? and germline_resource
    germline_resource = vcf_file reference, germline_resource if germline_resource
    germline_resource = GATK.prepare_VCF germline_resource if germline_resource

    pon = nil if pon == 'none'
    pon = Open.read(pon).strip if pon && Open.exists?(pon) && File.size(pon) < 10

    if pon == 'default' || pon == 'Broad'
      pon = case reference.to_s
            when 'b37'
              Organism["Hsa"].b37.known_sites["panel_of_normals.vcf"].produce.find
            when 'hg38'
              Organism["Hsa"].hg38.known_sites["panel_of_normals.vcf"].produce.find
              raise "No default or Broad panel of normals for reference '#{reference}'"
            end
    end
    pon = GATK.prepare_VCF pon if pon

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    args = {}
    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    tumor_sample = GATK.BAM_sample_name(tumor, reference)
    normal_sample = GATK.BAM_sample_name(normal, reference) if normal

    raise "No normal sample name" if normal and normal_sample.nil?

    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    output = file('calls.vcf')
    args["input"] = [tumor, normal].compact
    args["output"] = output
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample if normal_sample
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = interval_padding || GATKShard::GAP_SIZE if interval_list
    args["panel-of-normals"] = pon if pon
    #args["bam-output"] = file('haplotype.bam')
    args["germline-resource"] = germline_resource
    args["f1r2-tar-gz"] = file('f1r2.tar.gz')
    args["dont-use-soft-clipped-bases"] = remove_soft_clip if remove_soft_clip
    args["max-mnp-distance"] = max_mnp_distance

    # UPDATE FOR GATK 4.1.2
    #args["af-of-alleles-not-in-resource"] = "%.10f" % af_not_in_resource.to_s if af_not_in_resource

    shard = config('shard', :mutect2, :mutect, :gatk)

    if shard == 'true'
      contigs = Samtools.reference_contigs reference

      cpus = config('cpus', :mutect2, :mutect, :shard)

      headervcf = file('tmp.header')
      args["interval-padding"] ||= GATKShard::GAP_SIZE 
      args["intervals"] ||= nil
      args["bam-output"] = nil
      args["f1r2-tar-gz"] = '[OUTPUT]-f1r2.tar.gz'
      bar = self.progress_bar("Processing Mutect2 sharded")

      intervals = (interval_list || intervals_for_reference(reference))
      contentvcf = file('tmp.content')
      headervcf_stats = file('tmp.header.stats')
      contentvcf_stats = file('tmp.content.stats')
      tmp_stats = file('tmp.stats')

      Open.mkdir file('f1r2.tar.gz')
      GATKShard.cmd("Mutect2", args, intervals, GATKShard::CHUNK_SIZE, cpus, contigs, bar) do |ioutfile|
        bar.tick
        `grep "#" "#{ioutfile}" > "#{headervcf}"` unless File.exists? headervcf
        `grep -v "#" #{ioutfile} >> #{contentvcf}` 
        `head -n 1 "#{ioutfile}.stats" > "#{headervcf_stats}"` unless File.exists? headervcf_stats
        `tail -n 1 #{ioutfile}.stats >> #{contentvcf_stats}` 
        `mv "#{ioutfile}-f1r2.tar.gz" "#{file('f1r2.tar.gz')}"`
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

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  dep :contamination, :tumor_bam => :tumor, :normal_bam => :normal, :compute => :canfail
  dep :BAM_orientation_model
  dep :mutect2_pre
  extension :vcf
  task :mutect2_filters => :text do |reference|
    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2_pre).path, tmp
    FileUtils.ln_s step(:mutect2_pre).file('calls.vcf.stats'), tmp + '.stats'

    FileUtils.mkdir_p files_dir
    args = {}

    contamination = step(:contamination)

    output = file('output')
    args["variant"] = tmp
    args["output"] = output
    args["reference"] = reference
    args["--orientation-bias-artifact-priors"] = step(:BAM_orientation_model).path
    if contamination.error? || contamination.path.read.include?("NaN")
      set_info :missing_contamination, true
      Log.warn "NaN in contamination file or file missing: #{Log.color :blue, self.path}"
    else
      args["tumor-segmentation"] = contamination.file('segments.tsv')
      args["contamination-table"] = contamination.path
    end
    gatk("FilterMutectCalls", args)
    Open.link output, self.tmp_path
    nil
  end

  dep :mutect2_filters
  extension :vcf
  task :mutect2_clean => :text do
    TSV.traverse step(:mutect2_filters), :into => :stream, :type => :array do |line|
      next line if line[0] =~ /^#/
      
      chr = line.split("\t").first
      next unless chr =~ /^(chr)?[0-9MTXY]+$/
      next unless line =~ /PASS/

      line
    end
  end

  extension :vcf
  dep_task :mutect2, HTS, :mutect2_clean

end
