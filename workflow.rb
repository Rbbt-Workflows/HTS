require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/HTS'

require 'tools/BWA'
require 'tools/GATK'
require 'tools/samtools'
require 'tools/NovoAlign'
require 'tools/control_FREEC'
require 'tools/svABA'
require 'tools/stringTie'
require 'tools/HISAT'
#require 'tools/pfff'
require 'sources/HTS_reference'

module HTS
  extend Workflow
  
  #helper :reference_file do |reference|
  #  case reference
  #  when 'hg19'
  #    BWA.references.hg19[reference + ".fa"].produce.find
  #  when 'hg38'
  #    BWA.references.hg38[reference + ".fa"].produce.find
  #  when 'b37'
  #    BWA.references.b37[reference + ".fa"].produce.find
  #  when 'GRCh38'
  #    BWA.references.GRCh38[reference + ".fa"].produce.find
  #  when 'hs37d5'
  #    BWA.references.hs37d5[reference + ".fa"].produce.find
  #  else
  #    reference
  #  end
  #end

  helper :reference_file do |reference|
    reference = Open.read(reference).strip if String === reference && File.exists?(reference) && File.size(reference) < 2000
    case reference
    when 'hg19', 'hg38', 'b37', 'hs37d5'
      Organism["Hsa"][reference][reference + ".fa"].produce.find
    else
      reference
    end
  end

  helper :vcf_file do |reference, file|
    reference = reference.scan(/(?:b37|hg19|hg38)/).first 
    file = Open.read(file).strip if String === file && File.exists?(file) && File.size(file) < 2000

    case file.to_s.downcase
    when 'miller_indels', 'miller'
      Organism["Hsa"][reference].known_sites["Miller_1000G_indels.vcf.gz"].produce.find
    when '1000g_indels'
      Organism["Hsa"][reference].known_sites["1000G_phase1.indels.vcf.gz"].produce.find
    when '1000g_snps_hc', '1000g_snps', '1000g'
      Organism["Hsa"][reference].known_sites["1000G_phase1.snps.high_confidence.vcf.gz"].produce.find
    when 'gnomad'
      Organism["Hsa"][reference].known_sites["af-only-gnomad.vcf.gz"].produce.find
    when 'dbsnp'
      if reference == "hg38"
        Organism["Hsa"][reference].known_sites["dbsnp_146.vcf.gz"].produce.find
      else
        Organism["Hsa"][reference].known_sites["dbsnp_138.vcf.gz"].produce.find
      end
    when 'none'
      nil
    else
      file
    end
  end

  helper :germline_min_af do |file|
    file = Open.read(file).strip if String === file && File.exists?(file) && File.size(file) < 2000
    case file.to_s.downcase
    when 'gnomad'
      0.0000025
    when 'miller_indels', 'miller'
      nil
    when '1000g_indels'
      nil
    when '1000g_snps_hc', '1000g_snps', '1000g'
      nil
    when 'dbsnp'
      nil
    when 'none'
      nil
    else
      nil
    end
  end

  #helper :intervals_for_reference do |reference|
  #  fai = reference + '.fai'

  #  intervals = StringIO.new 
  #  TSV.traverse fai, :type => :array do |line|
  #    chr, size, rest = line.split("\t")
  #    next if BLACKLISTED_CONTIGS.select{|c| Regexp === c ? c.match(chr) : c === chr }.any?
  #    intervals << ([chr, "1", size] * "\t") << "\n"
  #  end

  #  intervals.rewind
  #  intervals
  #end

  BLACKLISTED_CONTIGS = ['hs37d5', /GL/, /NC_/, /KI/, /decoy/, /chrUn/, /random/]
  helper :intervals_for_reference do |reference|
    output = reference + '.byNS.interval_list' 
    unless File.exists?(output)
      args = {}
      args["REFERENCE"] = reference
      args["OUTPUT"] = output
      args["MAX_TO_MERGE"] = GATKShard::GAP_SIZE
      gatk("ScatterIntervalsByNs", args)
    end

    StringIO.new(Open.read(output).split("\n").reject{|line|
      if line =~ /^@SQ/ 
        true
      else
        chr = line.split("\t").first
        if BLACKLISTED_CONTIGS.select{|c| Regexp === c ? c.match(chr) : c === chr }.any?
          true
        else
          line =~ /Nmer$/
        end
      end
    }.collect{|line| line.split("\t").values_at(0,1,2) * "\t"} * "\n")
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

  helper :fix_spark_args do |command,args|
    if Hash === args
      args_new = {}
      args.each do |k,v|
        k = k.downcase if k.length > 1
        k = k.gsub('_', '-')
        args_new[k] = v
      end
      args = args_new
    end

    case command.to_s
    when "BaseRecalibrator"
      args.delete_if{|k,v| k.include? 'interval'}
    when "SortSam", "ApplyBQSR"
      args['--create-output-bam-index'] = false
    when "MarkDuplicates"
      args['--create-output-bam-index'] = false
      args.delete_if{|k,v| k.include? "sort-order" }
      args.delete_if{|k,v| k.include? "metrics" }
      args.delete_if{|k,v| k.include? "create-index" }
    end

    args
  end

  helper :gatk_io do |command,args,sin=nil|

    if GATK::SPARK_COMMANDS.include?(command) and config('spark', :gatk, command) 
      args = fix_spark_args command, args
      command += "Spark"
    end

    GATK.run(command, args, sin)
  end

  helper :gatk do |command,args,sin=nil|

    if GATK::SPARK_COMMANDS.include?(command) and config('spark', :gatk, command) 
      args = fix_spark_args command, args
      command += "Spark"
    end

    GATK.run_log(command, args, sin)
  end

end

require 'HTS/tasks/BAM'
require 'HTS/tasks/sequenza'
require 'HTS/tasks/CNV'
require 'HTS/tasks/VCF'
require 'HTS/tasks/test'
require 'HTS/tasks/RNASeq'
require 'HTS/tasks/IGV'

require 'HTS/tasks/sample' if defined? Sample
