require 'tools/FuSeq'
require 'HTS/tasks/RNASeq/BAM'
require 'HTS/tasks/RNASeq/salmon'
require 'HTS/tasks/RNASeq/kallisto'
require 'HTS/tasks/RNASeq/FuSeq'
require 'HTS/tasks/RNASeq/splicing'
require 'HTS/tasks/RNASeq/tximport'

module HTS

  dep :stringtie
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :htseq_count => :tsv do |organism|
    CMD.cmd_log("htseq-count -s no -f bam '#{step(:RNA_BAM).path}' '#{step(:stringtie).path}'  > '#{self.tmp_path}'")
    nil
  end
end
