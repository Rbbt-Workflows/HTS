
module HTS
  CMD.tool :lofreq do
    'conda install -c bioconda lofreq'
  end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  task :lofreq => :text do |tumor,normal,reference|
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
      CMD.cmd_log(:lofreq, "somatic -n #{normal} -t #{tumor} -f #{reference} --threads #{cpus} -o #{output} -d #{dbsnp}")
    else
      CMD.cmd_log(:lofreq, "somatic -n #{normal} -t #{tumor} -f #{reference} --threads #{cpus} -o #{output}")
    end
    
    Dir.glob(File.join(output,'*'))
  end
end
