$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '../..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rbbt-util'
require 'test/unit'
require 'tools/shard'
require 'tools/samtools'
require 'tools/GATK'

class TestShard < Test::Unit::TestCase

  def _test_intervals
    intervals = ["chr1\t1\t13000", "chr1\t13001\t15000", "chr1\t16000\t18000", "chrM\t1\t7000\n", "chrM\t10000\t14000\n"]
    Log.with_severity 0 do
      iii intervals
      iii GATKShard.chunk_intervals(intervals, 5000, nil, false)
      iii GATKShard.chunk_intervals(intervals, 5000, nil, true)
      iii GATKShard.chunk_intervals(intervals, 10000, nil, true)
    end
  end
  
  def test_hg38
    Workflow.require_workflow "HTS"
    reference = HTS.helper :reference_file, 'hg38'
    reference = HTS.prepare_FASTA reference
    intervals = HTS.helper :intervals_for_reference, reference
    chunks = GATKShard.chunk_intervals intervals, GATKShard::CHUNK_SIZE / 10, nil, false
    assert !GATKShard.chunks_overlap?(chunks)
    chunks = Log.with_severity 0 do
      GATKShard.chunk_intervals intervals, GATKShard::CHUNK_SIZE / 5, nil, true
    end
    assert GATKShard.chunks_overlap?(chunks)
    GATKShard.chunk_sizes(chunks)
  end
  


  def __test_interval_file
    interval_file = Path.setup("/data/Patients/Bellmunt/data/Padded.b37.bed")
    intervals = GATKShard.chunk_intervals interval_file
    intervals.each do |i|
      iif i
    end
    iii intervals.length
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

