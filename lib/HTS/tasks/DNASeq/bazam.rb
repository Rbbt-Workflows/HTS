module HTS
  input :bam, :file, "BAM file", nil, :nofile => true
  input :regions, :file, "BED file", nil, :nofile => true
  input :region_pad, :integer, "BED region padding", 0
  task :bazam_revert_BAM => :array do |bam,regions,pad|
    bam = Samtools.prepare_BAM(bam)
    f1 = file('fastq_1.fq.gz')
    f2 = file('fastq_2.fq.gz')
    Open.mkdir files_dir
    if regions
      CMD.cmd_log(:bazam, "-bam '#{bam}' -L '#{regions}' -pad #{pad} -r1 '#{f1}' -r2 '#{f2}'")
    else
      CMD.cmd_log(:bazam, "-bam '#{bam}' -r1 '#{f1}' -r2 '#{f2}'")
    end
    Dir.glob(File.join(files_dir, '*'))
  end

  input :bam, :file, "BAM file", nil, :nofile => true
  input :regions, :file, "BED file", nil, :nofile => true
  input :region_pad, :integer, "BED region padding", 0
  task :bazam_revert_BAM_RG => :array do |bam,regions,pad|
    bam = Samtools.prepare_BAM(bam)

    Open.mkdir files_dir
    io = if regions
           CMD.cmd(:bazam, "-bam '#{bam}' -L '#{regions}' -pad #{pad} ", :pipe => true)
         else
           CMD.cmd(:bazam, "-bam '#{bam}' ", :pipe => true)
         end

    # ADAPTED FROM
    # https://raw.githubusercontent.com/bsmn/fastq-split-readgroup-tool/master/bin/split-read-group
    awk_script =<<-'EOF'
{
  header=$1;
  sub(/^@/,"",header);

  x = split(header, arr1, " ");
  part1 = arr1[1];
  part2 = arr1[2];

  l = split(part1, arr2, ":");
  MACHINE = arr2[1];
  RUN = arr2[2];
  FLOWCELL = arr2[3];
  LANE = arr2[4];

  m = split(part2, arr3, ":")
  READ = arr3[1];

  print $1 "\n" $2 "\n" $3 "\n" $4 > OUTPREFIX "_" MACHINE "_" RUN "_R" READ ".fastq"
}
      EOF

    TmpFile.with_file(awk_script) do |script|
      CMD.cmd_log("paste - - - - | awk -F\"\\t\" -v OUTPREFIX='#{files_dir}/RG' -f #{script} ", :in => io)
    end
    CMD.cmd_log("gzip #{Dir.glob(File.join(files_dir, '*')) * " "}")
    Dir.glob(File.join(files_dir, '*.gz'))
  end
end
