module Study

  dep Sample, :normal_mutect2_for_panel do |jobname,options|
    jobname.split(",").collect do |study|
      samples = Sample.study_samples(study.strip).uniq
      samples.select{|sample| ! sample.include? 'normal'}.
        collect do |sample|
        {:jobname => sample.split(":").last}
      end
    end.flatten
  end
  input :interval_list, :file, "Interval list", nil, :nofile => true
  extension :vcf
  task :panel_of_normals => :text do |intervals|

    orig_reference = dependencies.first.recursive_inputs[:reference]
    type_of_sequencing = dependencies.first.recursive_inputs[:type_of_sequencing] if dependencies.first.recursive_inputs.fields.include?(:type_of_sequencing)
    type_of_sequencing = 'WGS' if type_of_sequencing.nil?
    intervals = HTS.helpers[:capture_intervals].call orig_reference, type_of_sequencing if intervals.nil?

    orig_reference = HTS.helpers[:reference_file].call(orig_reference)


    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference
    reference = HTS.uncompress_FASTA orig_reference

    bed_dir = file('bed')
    work_dir = file('work')
    args = {}
    args["variant"] = []

    Misc.in_dir bed_dir do
      dependencies.each do |dep|
        vcf = GATK.sort_VCF(dep.path, bed_dir)
        gvcf = HTS.prepare_BED(vcf, bed_dir)
        args["variant"].push(gvcf)
      end
    end

    Misc.in_dir work_dir do
      args["reference"] = reference
      args["genomicsdb-workspace-path"] = "./pon_db"
      args["merge-input-intervals"] = "TRUE"
      args["intervals"] = intervals
      GATK.run_log("GenomicsDBImport", args)

      args = {}
      args["reference"] = reference
      args["variant"] = "gendb://pon_db"
      args["create-output-variant-index"] = "false"
      args["output"] = self.tmp_path
      GATK.run_log("CreateSomaticPanelOfNormals", args)
      Open.rm_rf bed_dir
    end
    nil
  end
end
