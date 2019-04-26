require 'tools/control_FREEC'
require 'pathname'
require 'fileutils'
module HTS

	input :outputDir, :file, "output folder path. If relative, it will be stored in job results folder path", nil, :nofile => true
	input :sample_mateFile, :file, "File with mapped reads from tumor sample", nil, :nofile => true
	input :control_mateFile, :file, "File with mapped reads from control sample", nil, :nofile => true
	input :snpFile, :file, "File containing known SNP", nil, :nofile => true
	input :makePileup, :file, "SNP positions to create a mini pileup file from the initial BAM file provided in mateFileFolder", nil, :nofile => true
	input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38s GRCh38 hs37d5), :nofile => true
	input :captureRegions, :file, "List of Intervals for the capture regions", nil, :nofile => true
  extension :vcf
  task :control_freeC => :text do |outputDir, sample_mateFile, control_mateFile, snpFile, makePileup, reference, captureRegions|
		output = (Pathname.new outputDir).relative? ? file(outputDir): outPutDir
		FileUtils.mkdir_p output

		referencePath = reference_file(reference)	
                chrFilesPath = File.dirname(reference_file(reference)) + "s"
		chrLenFilePath = referencePath + ".fai" 

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
sex=XY
breakPointType=2
breakPointThreshold=0.8
noisyData=TRUE
printNA=FALSE
readCountThreshold=50

[sample]

mateFile = #{sample_mateFile} 
inputFormat = BAM
mateOrientation = FR

[control]

mateFile = #{control_mateFile} 
inputFormat = BAM
mateOrientation = FR

[BAF]

SNPfile = #{snpFile}
makePileup = #{makePileup}
minimalCoveragePerPosition = 0
fastaFile = #{referencePath} 

[target]

captureRegions = #{captureRegions}
		EOF
		
    IO::write(output + "/config_file.txt", script)
    ControlFREEC.run(output + "/config_file.txt")
    "DONE"
  end

  dep :control_freeC
  extension :png
  task :plot_freec_results => :binary do
    output = step(:control_freeC)
    ControlFREEC.makegraphs(output)
  end
end
