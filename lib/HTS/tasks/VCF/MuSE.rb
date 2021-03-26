module HTS
  CMD.tool :MuSE, nil, "bash -c 'type MuSE'" do
    CMD.cmd('conda install muse -c bioconda')
  end

  #input :tumor, :file, "Tumor BAM", nil, :nofile => true
  #input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  #input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  #extension :vcf
  #task :muse => :text do |tumor,normal, reference|
  #  reference = reference_file reference
  #  orig_reference = reference

  #  reference = GATK.prepare_FASTA orig_reference
  #  reference = Samtools.prepare_FASTA orig_reference

  #  reference = HTS.uncompress_FASTA orig_reference

  #  tumor = tumor.path if Step === tumor
  #  normal = normal.path if Step === normal

  #  tumor = Samtools.prepare_BAM(tumor)
  #  normal = Samtools.prepare_BAM(normal) if normal

  #  Open.mkdir files_dir
  #  int_file = file('intermediate')
  #  CMD.cmd_log(:MuSE, "call -O #{int_file} -f #{reference} #{tumor} #{normal}")

  #  CMD.cmd_log(:MuSE, "sump -I #{int_file}.MuSE.txt -E -O #{self.tmp_path}")
  #  nil
  #end

  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Normal BAM (optional)", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :type_of_sequencing, :select, "Whole genome or whole exome", nil, :select_options => %w(WGS WES panel)
  extension :vcf
  task :muse => :text do |tumor,normal, reference,type_of_sequencing|

    reference = reference_file reference
    orig_reference = reference

    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    reference = HTS.uncompress_FASTA orig_reference

    tumor = tumor.path if Step === tumor
    normal = normal.path if Step === normal

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    Open.mkdir files_dir
    int_file = file('intermediate')

    cpus = config :cpus, :muse, :default => 3
    cpus = 1 if cpus.nil?

    if cpus.to_s == "1"
      CMD.cmd_log(:MuSE, "call -O #{int_file} -f #{reference} #{tumor} #{normal}")
    else

      interval_list = intervals_for_reference(reference)
      chunk_size = 10_000_000
      chunks = GATKShard.chunk_intervals(interval_list, chunk_size)
      chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 5) if chunks.length < (cpus.to_i * 2) || chunks.length < 25
      chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 20) if chunks.length < (cpus.to_i * 2) || chunks.length < 25

      q = RbbtProcessQueue.new cpus

      q.callback do |output|
        Open.open(int_file + '.MuSE.txt', :mode => 'a') do |file|
          file.puts Open.read(output + '.MuSE.txt')
        end
      end

      bar = self.progress_bar("Processing MuSE intervals")

      bar.max = chunks.length if bar
      bar.init if bar

      TmpFile.with_file do |workdir|
        Open.mkdir workdir

        q.init do |intervals|
          Log.low "MuSE processing intervals #{Misc.fingerprint intervals}"
          name = [intervals.first * "__", intervals.last * "__"] * ","
          output = File.join(workdir, name)
          interval_file = File.join(workdir, name + '.intervals')

          Open.write(interval_file, intervals.collect{|chr,start,eend| "#{chr}:#{start}-#{eend.to_i + 1}" } * "\n")

          CMD.cmd_log(:MuSE, "call -O #{output} -f #{reference} #{tumor} #{normal} -l #{interval_file}")

          output
        end

        TSV.traverse chunks, :type =>:array do |intervals|
          q.process intervals
        end

        q.join
      end

      bar.done if bar
    end

    CMD.cmd_log(:MuSE, "sump -I #{int_file}.MuSE.txt #{type_of_sequencing.to_s == "WGS" ? "-G" : "-E" } -O #{self.tmp_path}")
    nil
  end
end
