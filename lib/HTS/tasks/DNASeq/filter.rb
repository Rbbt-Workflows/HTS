module HTS

  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :fastq1, :file, "FASTQ file"
  input :fastq2, :file, "FASTQ 2 file", nil
  extension "fastq.gz"
  task :bwa_filter => :binary do |reference,fastq1,fastq2|

    reference = reference_file reference
    #reference = BWA.prepare_FASTA reference
    reference = HTS.prepare_FASTA reference
    bwa_mem_args = " -t " << config('cpus', 'bwa')
    if fastq2
      CMD.cmd_log("bwa mem #{bwa_mem_args} '#{reference}' '#{fastq1}' '#{fastq2}' | samtools fastq -F4 - | gzip > #{self.tmp_path}")
    else
      CMD.cmd_log("bwa mem #{bwa_mem_args} '#{reference}' '#{fastq1}' | samtools fastq -F4 - | gzip > #{self.tmp_path}")
    end
    nil
  end

  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :fastq1, :file, "FASTQ file"
  input :fastq2, :file, "FASTQ 2 file", nil
  input :max_missmatches, :integer, "Maximum number of allowed missmatches", 2
  extension "fastq.gz"
  task :bowtie_filter => :binary do |reference,fastq1,fastq2,max_missmatches|

    reference_index_dir = file('reference_index') if Step === reference

    reference = reference_file reference
    reference = Bowtie.prepare_FASTA reference, reference_index_dir
    cpus = Rbbt::Config.get(:cpus, :bowtie) || 1

    original_fastq1 = fastq1
    original_fastq2 = fastq2 if fastq2

    fastq1 = file('f1.fastq.gz')
    fastq2 = file('f2.fastq.gz') if fastq2

    Open.ln_s original_fastq1, fastq1
    Open.ln_s original_fastq2, fastq2 if fastq2

    bowtie_args = '--end-to-end'
    missmatches = (0..max_missmatches).to_a.collect{|i| i.to_s}
    if fastq2
      CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -1 '#{fastq1}' -2 '#{fastq2}' | grep '^@\\|XM:i:[#{missmatches *","}]\\s' | samtools view -h -u -F4 - | samtools bam2fq - |gzip > '#{self.tmp_path}'  ")
    else
      CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -U '#{fastq1}' | grep '^@\\|XM:i:[#{missmatches *","}]\\s' | samtools view -h -u -F4 - | samtools bam2fq - | gzip > '#{self.tmp_path}'  ")
    end
    nil
  end


  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :fastq1, :file, "FASTQ file"
  input :fastq2, :file, "FASTQ file 2", nil
  extension "fastq.gz"
  task :razers3_filter => :binary do |reference,fastq1,fastq2|
    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    reference = reference_file reference
    reference = HTS.prepare_FASTA reference
    cpus = Rbbt::Config.get(:cpus, :razers)
    TmpFile.with_file :extension => :bam do |tmp_bam|
      if fastq2
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}' '#{fastq2}'")
      else
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}'")
      end
      CMD.cmd("samtools bam2fq -F 4 '#{tmp_bam}' | gzip > '#{self.tmp_path}'")
    end
    nil
  end

  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :fastq1, :file, "FASTQ file"
  input :fastq2, :file, "FASTQ file 2", nil
  input :novoalign_args, :string, "NovoAlign arguments", "-F STDFQ -R 0 -r All 9999 -o SAM -o FullNW"
  extension "fastq.gz"
  task :novoalign_filter => :binary do |reference,fastq1,fastq2,novoalign_args|
    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    reference = reference_file reference
    reference = NovoAlign.prepare_FASTA reference
    cpus = Rbbt::Config.get(:cpus, :novoalign)
    TmpFile.with_file :extension => :sam do |tmp_sam|
      if fastq2
        CMD.cmd_log(:novoalign, "-d '#{reference}.nix' -f '#{fastq1}' '#{fastq2}' #{novoalign_args} 1> #{tmp_sam}")
      else
        CMD.cmd_log(:novoalign, "-d '#{reference}.nix' -f '#{fastq1}' #{novoalign_args} 1> #{tmp_sam}")
      end
      CMD.cmd("samtools sam2fq -F 4 '#{tmp_sam}' | gzip > '#{self.tmp_path}'")
    end
    nil
  end
end
