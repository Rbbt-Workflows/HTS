require 'tools/shard'
module HTS

  helper :intervals_for_reference do |reference|
    fai = reference + '.fai'

    intervals = StringIO.new 
    TSV.traverse fai, :type => :array do |line|
      chr, size, rest = line.split("\t")
      intervals << ([chr, "1", size] * "\t") << "\n"
    end

    intervals.rewind
    intervals
  end

  helper :germline_resource_file do |germline_resource,reference|
    germline_resource = germline_resource.to_s if Symbol === germline_resource
    case germline_resource
    when 'gnomad'
      GATK.known_sites[reference]["af-only-gnomad.vcf.gz"].produce.find
    when String
      return germline_resource if Misc.is_filename?(germline_resource) && Open.exists?(germline_resource)
      gr = GATK.known_sites[reference].glob("*" << germline_resource << "*").first
      if gr
        gr
      else
        germline_resource
      end
    else
      germline_resource
    end
  end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :interval_list, :file, "Interval list", nil, :nofile => true
  input :pon, :file, "Panel of normals", nil, :nofile => true
  input :germline_resource, :file, "Germline resource", :gnomad, :nofile => true
  input :af_not_in_resource, :float, "Allele frequency of alleles not in resource", nil
  extension :vcf
  task :mutect2 => :text do |tumor,normal,reference,interval_list,pon,germline_resource,af_not_in_resource|

    germline_resource = germline_resource_file germline_resource, reference

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

    args["input"] = [tumor, normal].compact
    args["output"] = self.tmp_path
    args["reference"] = reference
    args["tumor-sample"] = tumor_sample
    args["normal-sample"] = normal_sample if normal_sample
    args["intervals"] = interval_list if interval_list
    args["interval-padding"] = 100 if interval_list
    args["panel-of-normals"] = pon if pon
    args["bam-output"] = file('haplotype.bam')
    args["germline-resource"] = germline_resource
    args["af-of-alleles-not-in-resource"] = af_not_in_resource.to_s if af_not_in_resource

    shard = Rbbt::Config.get('shard', :gatk, :mutect, :mutect2)

    if shard == 'true'
      headervcf = file('tmp.header')
      contentvcf = file('tmp.content')
      args["intervals"] ||= nil
      intervals = (interval_list || intervals_for_reference(reference))
      bar = self.progress_bar("Processing Mutect2 sharded")

      bar.init
      GATKShard.cmd("Mutect2", args, intervals, 30_000_000) do |ioutfile|
        bar.tick
        `grep "#" "#{ioutfile}" > "#{headervcf}"` unless File.exists? headervcf
        `grep -v "#" #{ioutfile} >> #{contentvcf}` 
      end
      bar.done

      `cat "#{headervcf}" > "#{self.tmp_path}"`
      `sort -k 1,2 "#{contentvcf}" | uniq >> "#{self.tmp_path}"`
      nil
    else
      GATK.run_log("Mutect2", args)
    end
  end

  dep :mutect2
  dep :contamination, :BAM => :tumor
  extension :vcf
  task :mutect2_filtered => :text do
    args = {}
    FileUtils.mkdir_p files_dir

    tmp = TmpFile.tmp_file
    FileUtils.ln_s step(:mutect2).path, tmp

    args["variant"] = tmp
    args["output"] = self.tmp_path
    args["contamination-table"] = step(:contamination).path
    GATK.run_log("FilterMutectCalls", args)
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
end
