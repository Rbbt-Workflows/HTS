module HTS

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  extension :bam
  task :razers3_BAM => :text do |reference,fastq1,fastq2|

    reference = reference_file reference
    reference = BWA.prepare_FASTA reference
    cpus = Rbbt::Config.get(:cpus, :razers)
    TmpFile.with_file :extension => :bam do |tmp_bam|
      if fastq2
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}' '#{fastq2}'")
      else
        CMD.cmd_log("#{Rbbt.software.opt.seqan.bin.razers3.find} -tc #{cpus} -i 95 -m 1 -dr 0 -o '#{tmp_bam}' '#{reference}' '#{fastq1}'")
      end
      CMD.cmd("samtools sort '#{tmp_bam}' > '#{self.tmp_path}'")
    end
    nil
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  input :bowtie_args, :string, "Bowtie2 arguments", "--end-to-end"
  extension :bam
  task :bowtie_BAM => :text do |reference,fastq1,fastq2,bowtie_args|

    reference_index_dir = file('reference_index') if Step === reference

    reference = reference_file reference
    reference = Bowtie.prepare_FASTA reference, reference_index_dir
    cpus = Rbbt::Config.get(:cpus, :bowtie) || 1

    original_fastq1 = fastq1
    original_fastq2 = fastq2 if fastq2

    fastq1 = file('f1.fastq')
    fastq2 = file('f2.fastq') if fastq2

    fastq1 += '.gz' if original_fastq1 =~ /\.gz$/
    fastq2 += '.gz' if fastq2 && original_fastq2 =~ /\.gz$/

    Open.ln_s original_fastq1, fastq1
    Open.ln_s original_fastq2, fastq2 if fastq2

    TmpFile.with_file :extension => :bam do |tmp_bam|
      if fastq2
        CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -1 '#{fastq1}' -2 '#{fastq2}' -S '#{tmp_bam}'  ")
      else
        CMD.cmd_log("#{Rbbt.software.opt.Bowtie2.bowtie2.find} #{bowtie_args} -p #{cpus} -x '#{reference}' -U '#{fastq1}' -S '#{tmp_bam}'")
      end
      CMD.cmd("samtools sort '#{tmp_bam}' > '#{self.tmp_path}'")
    end
    nil
  end

  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  input :novoalign_args, :string, "NovoAlign arguments", "-F STDFQ -R 0 -r All 9999 -o SAM -o FullNW"
  extension "bam"
  task :novoalign_BAM => :binary do |reference,fastq1,fastq2,novoalign_args|
    FileUtils.mkdir_p files_dir unless File.exists? files_dir

    reference = reference_file reference
    cpus = Rbbt::Config.get(:cpus, :novoalign)
    begin
      original_fastq1 = fastq1
      original_fastq2 = fastq2 if fastq2
      original_reference = reference

      reference = file('reference.fasta').find if original_reference =~ /\.gz/
      if original_reference != reference
        CMD.cmd("gunzip '#{original_reference}' -c > '#{reference}'")

        CMD.cmd("'#{Rbbt.software.opt.NovoAlign.novoindex.find}' '#{ reference }.nix' '#{reference}'")
      else
        reference = NovoAlign.prepare_FASTA reference 
      end

      fastq1 = file('file.1.fastq').find if original_fastq1 =~ /\.gz/
      fastq2 = file('file.2.fastq').find if fastq2 && original_fastq2 =~ /\.gz/

      CMD.cmd("gunzip '#{original_fastq1}' -c > '#{fastq1}'") if original_fastq1 != fastq1
      CMD.cmd("gunzip '#{original_fastq2}' -c > '#{fastq2}'") if fastq2 && original_fastq2 != fastq2


      TmpFile.with_file :extension => :sam do |tmp_sam|
        if fastq2
          CMD.cmd_log("#{Rbbt.software.opt.NovoAlign.novoalign.find} -d '#{reference}.nix' -f '#{fastq1}' '#{fastq2}' #{novoalign_args} 1> #{tmp_sam}")
        else
          CMD.cmd_log("#{Rbbt.software.opt.NovoAlign.novoalign.find} -d '#{reference}.nix' -f '#{fastq1}' #{novoalign_args} 1> #{tmp_sam}")
        end
        CMD.cmd("samtools sort -O BAM  '#{tmp_sam}' > '#{self.tmp_path}'")
      end
    ensure
      FileUtils.rm fastq1 if File.exists?(fastq1) && original_fastq1 != fastq1
      FileUtils.rm fastq2 if fastq2 && File.exists?(fastq2) && original_fastq2 != fastq2
      FileUtils.rm reference if File.exists?(reference) && original_reference != reference
      FileUtils.rm reference + '.nix' if File.exists?(reference + '.nix') && original_reference != reference
    end
    nil
  end

end
