require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/HTS'

require 'tools/BWA'
require 'tools/GATK'
require 'tools/samtools'
require 'tools/NovoAlign'
require 'tools/Bowtie'

module HTS
  extend Workflow
  
  helper :reference_file do |reference|
    case reference
    when 'hg19'
      BWA.references.hg19["reference.fa"].produce.find
    when 'hg38'
      BWA.references.hg38["reference.fa"].produce.find
    when 'b37'
      BWA.references.b37["reference.fa"].produce.find
    when 'GRCh38'
      BWA.references.GRCh38["reference.fa"].produce.find
    else
      reference
    end
  end

  helper :monitor_genome do |stream,bgzip=true|
    chromosome_sizes = Persist.persist('chromosome_sizes', :yaml, :organism => nil) do
      Organism.chromosome_sizes
    end

    total = Misc.sum(chromosome_sizes.collect{|k,v| v}).to_i

    chromosomes = chromosome_sizes.keys.sort{|a,b| Misc.chr_cmp_strict(a,b)}

    chromosome_offset = {}
    chromosomes.each_with_index do |chr,i|
      offset = 0
      chromosomes[0..i-1].each{|c| offset += chromosome_sizes[c] } if i > 0
      chromosome_offset[chr] = offset
    end

    self.monitor_stream stream, :bar => {:desc => task_name.to_s, :max => total}, :bgzip => bgzip do |line,bar|
      chr, pos = line.split(/[\t:]/)
      pos = pos.to_i
      if pos > 0 and chromosome_offset.include?(chr)
        position = chromosome_offset[chr] + pos
        bar.pos(position.to_i)
      end
    end
  end

  helper :monitor_cmd_genome do |cmd,fifo=false,bgzip=false|
    cmd, args = cmd if Array === cmd
    args = {} if args.nil?
    if fifo
      fifo = file('fifo') if TrueClass === fifo
      if cmd.include? "{FIFO}"
        cmd = cmd.sub('{FIFO}',"'#{fifo}'")

        Misc.with_fifo(fifo) do |fifo|
          io = CMD.cmd(cmd, args.merge(:pipe => true))
          res = self.monitor_genome io, bgzip 
        end
      else
        Misc.with_fifo(fifo) do |fifo|
          CMD.cmd(cmd, args.merge(:pipe => true))
          io = Open.open(fifo)
          res = self.monitor_genome io, bgzip 
          ConcurrentStream.setup res do
            Open.rm fifo if File.exists? fifo
          end
        end
      end
    else
      io = CMD.cmd(cmd, args.merge(:pipe => true))
      self.monitor_genome io, bgzip 
    end
  end

  helper :monitor_cmd do |cmd,bar,fifo=false,bgzip=false|
    cmd, args = cmd if Array === cmd
    args = {} if args.nil?
    if fifo
      if cmd.include? "{FIFO}"
        fifo = file('fifo') if TrueClass === fifo
        Misc.with_fifo(fifo) do |fifo|
          cmd = cmd.replace('{FIFO}',"'#{fifo}'")
          io = CMD.cmd(cmd, args.merge(:pipe => true))
          self.monitor_stream io, :bar => bar, :bgzip => true 
        end
      else
        cmd += " & while [ ! -f '#{fifo}' ]; do sleep 0.5; done; tail -f '#{fifo}'"
        io = CMD.cmd(cmd, args.merge(:pipe => true))
        self.monitor_genome io, bgzip 
      end
    else
      io = CMD.cmd(cmd, args.merge(:pipe => true))
      self.monitor_stream io, :bar => bar, :bgzip => true 
    end
  end

end

require 'HTS/tasks/BAM'
require 'HTS/tasks/filter'
require 'HTS/tasks/sequenza'
require 'HTS/tasks/MuTect2'
require 'HTS/tasks/VarScan'
require 'HTS/tasks/Strelka'
require 'HTS/tasks/test'

#require 'rbbt/knowledge_base/HTS'
#require 'rbbt/entity/HTS'

