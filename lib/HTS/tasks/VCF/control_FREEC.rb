require 'tools/control_FREEC'
require 'pathname'
require 'fileutils'
module HTS

	input :sample_mateFile, :file, "File with mapped reads from tumor sample", nil, :nofile => true
	input :control_mateFile, :file, "File with mapped reads from control sample", nil, :nofile => true
	input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38s GRCh38 hs37d5), :nofile => true
	input :snpFile, :file, "File containing known SNP", nil, :nofile => true
	input :intervals, :file, "List of Intervals for the capture regions", nil, :nofile => true
    input :mate_orientation,:string, "0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)" "FR" 
  dep :BAM_pileup_sumaries_known_biallelic, :jobname => "Default"
  extension :vcf
  task :control_freeC => :tsv do |sample_mateFile, control_mateFile, reference, snpFile, intervals, mate_orientation|
    output = file('output')
    Open.mkdir output


    snpFile = vcf_file reference, snpFile if snpFile
    referencePath = BWA.prepare_FASTA(reference_file(reference))
    makePileup = step(:BAM_pileup_sumaries_known_biallelic).path

    chrFilesPath = HTS.unfold_FASTA(referencePath)

    if intervals
      new = file('reference.fai')
      contigs = Open.read(intervals).split("\n").collect{|l| l.split("\t").first}.uniq
      Open.open(new, :mode => 'w') do |fout|
        TSV.traverse referencePath + ".fai",  :type => :array do |line|
          parts = line.split("\t")
          next unless contigs.include? parts.first
          fout.puts line
        end
      end
      chrLenFilePath = new
    else
      chrLenFilePath = referencePath + ".fai" 
    end

    sex = 'XX'

    script = <<-EOF
[general]
BedGraphOutput = TRUE
chrFiles = #{chrFilesPath}
chrLenFile = #{chrLenFilePath}
window = 0
ploidy = 2
outputDir = #{output}
maxThreads = 48
minimalSubclonePresence = 30
sex=#{sex}
breakPointType=2
breakPointThreshold=0.8
noisyData=FALSE
printNA=FALSE
readCountThreshold=50

[sample]

mateFile = #{sample_mateFile} 
inputFormat = BAM
mateOrientation = #{mate_orientation}

[control]

mateFile = #{control_mateFile} 
inputFormat = BAM
mateOrientation = #{mate_orientation}

[BAF]

SNPfile = #{snpFile}
makePileup = #{makePileup}
minimalCoveragePerPosition = 0
fastaFile = #{referencePath} 

[target]

captureRegions = #{intervals}
		EOF
		
    IO::write(output + "/config_file.txt", script)
    ControlFREEC.run(output + "/config_file.txt")
    result = output.glob("*.bam_ratio.txt").first
    TSV.open result, :type => :list, :header_hash => '', :key_field => 'Gene'
  end

  dep :control_freeC
  extension :png
  task :plot_freec_results => :string do
    output = file('output')
    input = step(:control_freeC).file('output')
    ControlFREEC.makegraphs(input, output)
    Open.cp output.glob("*.bam_ratio.txt.png").first, self.tmp_path
    nil
  end
end
