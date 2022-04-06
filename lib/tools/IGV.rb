require 'rbbt-util' 
require 'rbbt/resource'

module IGV

  def self.run(script, reference = nil, height = 5000, width = 2000)
    Log.debug "IGV script:\n" << script
    TmpFile.with_file(script) do |script_file|
      #CMD.cmd_log("xvfb-run --server-args='-screen 0 #{height}x#{width}x24' --auto-servernum --server-num=1 java -Xmx4000m -jar '#{Rbbt.software.opt.IGV.produce["igv.jar"].find}' -b '#{script_file}'")
      #CMD.cmd_log("xvfb-run --server-args='-screen 0 #{height}x#{width}x24' --auto-servernum --server-num=1 \
      #            java -Xmx4000m -jar '#{Rbbt.software.opt.IGV.produce["lib/igv.jar"].find}' -b '#{script_file}'")
      #CMD.cmd_log("java -Xmx4000m -jar '#{Rbbt.software.opt.IGV.glob("**/igv.jar").first}' -b '#{script_file}' -g '#{reference}'", :xvfb => true)
      CMD.cmd_log('igv.sh', "-b '#{script_file}' -g '#{reference}'", :xvfb => true)
    end
  end

  Rbbt.claim Rbbt.software.opt.IGV, :install, Rbbt.share.install.software.IGV.find
end

if __FILE__ == $0
  iif Rbbt.software.opt.IGV.produce.find

  Log.severity = 0

  IGV.run <<-EOF
new
genome hg19
snapshotDirectory /home/mvazque2/tmp/
load ~/software/scm/IGV-snapshot-automator/test_data/test_alignments.bam
load ~/software/scm/IGV-snapshot-automator/test_data/test_alignments2.bam
maxPanelHeight 500
goto chr1:713167-714758
snapshot chr1_713167_714758_h500.png
goto chr1:713500-714900
snapshot chr1_713500_714900_h500.png
goto chr1:714000-715000
snapshot chr1_714000_715000_h500.png
goto chr1:714001-714001
snapshot chr1_714001_714001_h500.png
exit
  EOF
end

