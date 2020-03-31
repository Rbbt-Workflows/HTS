require 'tools/FastQC'

module HTS

  input :fastq_file, :file, "List of FASTQ files to process"
  task :fastqc => :string do |fastq_file|

    Open.mkdir files_dir

    CMD.cmd_log(:fastqc, " --outdir '#{files_dir}' '#{fastq_file}'")

    zipfile = Dir.glob(files_dir + '/*').first
    name = File.basename(zipfile).sub('.zip', '')

    CMD.cmd("unzip -p #{zipfile} #{name}/summary.txt").read
  end
end
