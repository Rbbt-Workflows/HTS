require 'HTS/tasks/VCF/MuTect2'
require 'HTS/tasks/VCF/Haplotype'
require 'HTS/tasks/VCF/VarScan'
require 'HTS/tasks/VCF/Strelka'
require 'HTS/tasks/VCF/Delly'
require 'HTS/tasks/VCF/svABA'
require 'HTS/tasks/VCF/control_FREEC'
require 'HTS/tasks/VCF/MuSE'
require 'HTS/tasks/VCF/somaticsniper'

module HTS

  input :vcf_1, :file, "VCF file"
  input :vcf_2, :file, "VCF file"
  input :sort, :boolean, "Sort files", true
  task :compare_vcf => :tsv do |vcf_1,vcf_2,sort|
    first = []
    last = []
    common = []

    Misc.with_fifo do |pipe1|
      Misc.with_fifo do |pipe2|
        th1 = Thread.new do
          Open.write(pipe1) do |pipe1|
            TSV.traverse vcf_1, :type => :array do |line|
              next if line =~ /^#/ 
                chr, pos, id, ref, alt, *rest = line.split("\t")
              pipe1.puts [chr, pos, ref, alt] * ":" + "\t" + line.gsub("|", "---")
            end
          end
        end

        th2 = Thread.new do
          Open.write(pipe2) do |pipe2|
            TSV.traverse vcf_2, :type => :array do |line|
              next if line =~ /^#/ 
                chr, pos, id, ref, alt, *rest = line.split("\t")
              pipe2.puts [chr, pos, ref, alt] * ":" + "\t" + line.gsub("|", "---")
            end
          end
        end

        TSV.traverse TSV.paste_streams([pipe1, pipe2], :header_hash => "", :same_fields => true, :sort => sort), :type => :array do |line|
          next if line =~ /^(?:#|CHR)/
          mutation, chr, pos, id, *parts = line.split("\t", -1)
          case
          when id[0] == "|"
            first << mutation
          when id[-1] == "|"
            last << mutation
          else
            common << [mutation, parts]
          end
        end
        th1.join
        th2.join
      end
    end
    tsv = TSV.setup({}, "Statistic~Value#:type=:single")

    iif first
    iif last
    iii common.select{|mutation,parts| parts.select{|p| p.split("|").uniq.length != 1}.any? }
    tsv["Missing"] = first.length
    tsv["Extra"] = last.length
    tsv["Common"] = common.length
    tsv["Common but different"] = common.select{|mutation,parts| parts.select{|p| p.split("|").uniq.length != 1}.any?}.length
    tsv
  end
end
