$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rbbt-util'
require 'test/unit'
require 'tools/shard'
require 'tools/samtools'
require 'tools/GATK'

class TestShard < Test::Unit::TestCase
  def __test_interval_file
    interval_file = Path.setup("/data/Patients/Bellmunt/data/Padded.b37.bed")
    GATKShard.process_intervals interval_file, 1000 do |intervals|
      assert Array === intervals
      break
    end
  end
  
  def __test_interval_string
    interval_file = StringIO.new <<-EOF
1	1	10000000
1	1	10000000
2	1	10000000
3	1	10000000
EOF
    GATKShard.process_intervals interval_file, 100_000 do |intervals|
      iii intervals
    end
  end

  def __test_Mutect2
    Log.severity = 0

    tumor = "/home/mvazque2/.rbbt/var/jobs/HTSBenchmark/BAM/tumor.bam"
    normal = nil
    af_not_in_resource, pon, germline_resource = nil

    tumor = Samtools.prepare_BAM(tumor)
    normal = Samtools.prepare_BAM(normal) if normal

    tumor_sample = GATK.BAM_sample_name(tumor)
    normal_sample = GATK.BAM_sample_name(normal) if normal

    TmpFile.with_file(nil, false) do |output|
      Open.write(output) do |foutput|

        args = {}
        args["input"] = [tumor, normal].compact
        args["reference"] = '/home/mvazque2/.rbbt/var/fasta_indices/f03986b15fd8f1a60e3a8997c91e0a84/b37.fa'
        args["output"] = nil
        args["tumor-sample"] = tumor_sample
        args["normal-sample"] = normal_sample if normal_sample
        args["intervals"] = nil
        args["interval-padding"] = 100 
        args["panel-of-normals"] = pon if pon
        args["germline-resource"] = germline_resource if germline_resource
        args["af-of-alleles-not-in-resource"] = af_not_in_resource.to_s if af_not_in_resource

        GATKShatter.cmd 'Mutect2', args, interval_file, 100_000, 10 do |ioutput|
          `cat #{ioutput} >> #{output}`
        end
      end
      iii output
    end
  end
  
end

