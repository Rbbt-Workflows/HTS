
module HTS
  Rbbt.claim Rbbt.software.opt.LoFreq, :install do
    commands =<<-EOF
get_git $name $url
    EOF
    {:git => "https://github.com/CSB5/lofreq.git", :commands => commands}
  end

  CMD.tool :lofreq, Rbbt.software.opt.LoFreq

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :lofreq_pre => :text do |tumor,normal,reference|
    output = file('output')
    orig_reference = reference

    reference_code = File.basename(reference).sub('.fa','').sub('gz','')
    reference = reference_file reference

    normal = Samtools.prepare_BAM(normal) if normal
    tumor = Samtools.prepare_BAM(tumor) if tumor

    reference = Samtools.prepare_FASTA(reference)

    cpus = config :cpus, :lofreq, :LoFreq, :default => 1

    output = file('output')
    Open.mkdir output

    build = Organism.GRC_build reference_code
    if build
      dbsnp = DbSNP[build]["common.vcf.gz"].find
      dbsnp = GATK.prepare_VCF dbsnp
      CMD.cmd_log(:lofreq, "somatic -n #{normal} -t #{tumor} -f #{reference} --threads #{cpus} -o #{output}/ -d #{dbsnp}")
    else
      CMD.cmd_log(:lofreq, "somatic -n #{normal} -t #{tumor} -f #{reference} --threads #{cpus} -o #{output}/")
    end
    
    Dir.glob(File.join(output,'*'))
  end

  dep :lofreq_pre, :compute => :produce
  extension :vcf
  dep_task :lofreq_joined, HTS, :join_vcfs, :vcf1 => :placeholder, :vcf2 => :placeholder do |jobname,options,dependencies|
    lofreq_pre = dependencies.flatten.first
    output = lofreq_pre.file('output')
    options[:vcf1] = output["somatic_final.snvs.vcf.gz"]
    options[:vcf2] = output["somatic_final.indels.vcf.gz"]
    {:inputs => options}
  end

  dep :lofreq_joined
  extension :vcf
  task :lofreq => :text do 
    TSV.traverse step(:lofreq_joined), :type => :array, :into => :stream do |line|
      next line if line =~ /^##/
      next line + "\tFORMAT\tTumor" if line =~ /^#/

      parts = line.split("\t")
      info = parts.pop
      parts << ""
      format = []
      values = []
      info.split(";").each do |e|
        key, value = e.split("=")

        value = key if value.nil?

        format << key
        values << value 
      end

      (parts + [format*":", values*":"]) * "\t"
    end
  end
end
