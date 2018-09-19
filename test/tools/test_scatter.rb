require File.join(File.expand_path(File.dirname(__FILE__)), '..', 'test_helper.rb')
require 'rbbt-util'
require 'tools/scatter'
require 'rbbt/workflow'

class TestClass < Test::Unit::TestCase
  def setup
    Log.severity = 0

    @tumor = '/home/mvazque2/git/sandbox/HTS_test/subset.large/tumor.bam'
    @normal = '/home/mvazque2/git/sandbox/HTS_test/subset.large/normal.bam'

    @tumor = '/home/mvazque2/.rbbt/var/jobs/CFDTL/BAM/T16-2.bam'
    @normal = '/home/mvazque2/.rbbt/var/jobs/CFDTL/BAM/Blood_control.bam'

    @tumor = '/home/mvazque2/git/sandbox/HTS_test/subset.all/tumor.bam'
    @normal = '/home/mvazque2/git/sandbox/HTS_test/subset.all/normal.bam'


  end

  def _test_read_head
    TmpFile.with_file do |directory|
      s = Scatter.new @tumor, @normal, directory
      s.read_head
      assert s.instance_variable_get("@head2").include? "@RG"
    end
  end

  def _test_read_stream
    TmpFile.with_file do |directory|
      Path.setup(directory)
      Open.write(directory.foo) do |out|
        tumor = CMD.cmd("samtools view -h '#{@tumor}'", :pipe => true, :no_fail => true)
        current_chr, current_pos = Scatter.read_stream_chunk(tumor, out, 10_000_000)
        padding = Scatter.read_stream_padding(tumor, out, current_chr, current_pos + 10_000)
        padding_start = padding.split("\n").first.split("\t")[3]
        padding_end = padding.split("\n").last.split("\t")[3]
        assert padding_end.to_i - padding_start.to_i > 10_000
        assert File.size(directory.foo.find) > 10_000_000
      end
    end
  end


  def _test_prepare_files
    TmpFile.with_file do |directory|
      s = Scatter.new @tumor, @normal, directory
      s.read_head
      f11, f12 = s.prepare_files("run1")
      f21, f22 = s.prepare_files("run2")
      assert File.exists? f11
      assert File.exists? f12
      assert File.exists? f21
      assert File.exists? f22
    end
  end

  def _test_iterate
    TmpFile.with_file do |directory|
      begin
        s = Scatter.new @tumor, @normal, directory
        ls1,le1,ls2,le2 = nil,nil,nil,nil
        s.iterate do |f1,f2|
          s1 = CMD.cmd("grep -v ^@ '#{f1}' |head -n 1 | cut -f 4").read
          s2 = CMD.cmd("grep -v ^@ '#{f2}' |head -n 1 | cut -f 4").read
          e1 = CMD.cmd("tail -n 1 '#{f1}' | cut -f 4").read
          e2 = CMD.cmd("tail -n 1 '#{f2}' | cut -f 4").read
          if ls1.nil? 
            ls1,le1,ls2,le2 = s1,e1,s2,e2
          else
            assert le1 > s1 
            assert le1 > s2
            assert le2 > s1  
            assert le2 > s2
            raise Interrupt
          end
          Open.rm f1
          Open.rm f2
          sleep 1
        end
      rescue Interrupt
      end
    end
  end

  def test_process
    Workflow.require_workflow 'HTS'
    TmpFile.with_file(nil, false) do |directory|
      chromosome_sizes = Persist.persist('chromosome_sizes', :yaml, :organism => nil) do
        Organism.chromosome_sizes
      end

      total = Misc.sum(chromosome_sizes.collect{|k,v| v}).to_i

      chromosomes = chromosome_sizes.keys.sort{|a,b| Misc.chr_cmp_strict(a,b)}

      max_chr_pos = chromosome_sizes

      gspos = 1_000_000
      gepos = 1_200_000

      begin
        s = Scatter.new @tumor, @normal, directory, 20_000_000, 5_000
        s.process(15) do |tumor, normal|
          output = tumor + '.result.vcf'

          schr1, spos1 = CMD.cmd("samtools view '#{tumor}' | head -n 1 | cut -f 3,4").read.split("\t")
          echr1, epos1 = CMD.cmd("samtools view '#{tumor}' | tail -n 1 | cut -f 3,4").read.split("\t")
          iii [schr1, spos1]
          iii [echr1, epos1]

          schr2, spos2 = CMD.cmd("samtools view '#{tumor}' | head -n 1 | cut -f 3,4").read.split("\t")
          echr2, epos2 = CMD.cmd("samtools view '#{tumor}' | tail -n 1 | cut -f 3,4").read.split("\t")
          iii [schr2, spos2]
          iii [echr2, epos2]

          min_chr = [schr1, schr2].sort.first
          min_pos = [gspos, spos1.to_i, spos2.to_i].sort.first
          
          max_chr = [echr1, echr2].sort.last
          max_pos = [epos1.to_i, epos2.to_i, gspos].sort.last

          intervals = []
          if min_chr != max_chr
            chromosomes.select{|c| Misc.chr_cmp_strict(c, min_chr) == 1 && Misc.chr_cmp_strict(c, max_chr) == -1}.each do |c|
              intervals << [c, [gspos, gepos] * "-"] * ":"
            end
            #intervals << [min_chr, [min_pos, max_chr_pos[min_chr]] * "-"] * ":"
            intervals << [min_chr, [min_pos, gepos] * "-"] * ":"
            intervals << [max_chr, [gspos, max_pos] * "-"] * ":"
          else
            intervals << [min_chr, [min_pos, max_pos] * "-"] * ":"
          end

          intervals.select!{|i| i.split(":").first != "MT" }

          tumor_bam = tumor + '.bam'
          normal_bam = normal + '.bam'

          CMD.cmd("samtools view -h -b '#{tumor}' > '#{tumor_bam}'")

          CMD.cmd("samtools view -h -b '#{normal}' > '#{normal_bam}'")
          #Open.rm tumor
          #Open.rm normal

          tumor_bam = Samtools.prepare_BAM(tumor_bam, directory)
          normal_bam = Samtools.prepare_BAM(normal_bam, directory)

          tumor_sample = GATK.BAM_sample_name(tumor_bam)
          normal_sample = GATK.BAM_sample_name(normal_bam)

          reference = BWA.references.b37["reference.fa"].produce.find

          args = {}
          args["input"] = [tumor_bam, normal_bam]
          args["output"] = output
          args["reference"] = "/home/mvazque2/git/sandbox/HTS_test/reference/reference.fa"
          args["tumor-sample"] = tumor_sample
          args["bam-output"] = tumor + '.haplotype.bam'
          args["-L"] = intervals

          GATK.run_log("Mutect2", args)

          output
        end
      end
    end
  end
end

