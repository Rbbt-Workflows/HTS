class GATKShard

  GAP_SIZE = 1_000
  CHUNK_SIZE = 100_000_000

  def self.chunk_intervals(interval_list=nil, chunk_size=CHUNK_SIZE, contigs = nil, break_interval = false)
    intervals = TSV.traverse interval_list, :type => :array, :into => [] do |line|
      next if line =~ /^@/ 
      chr, start, eend, *rest = line.split("\t")
      next unless start =~ /^\d+$/
      next unless eend =~ /^\d+$/
      [chr, start.to_i, eend.to_i]
    end

    if contigs
      # Deal with contigs with ":" characters
      intervals = intervals.select{|i| contigs.include?(i.first) }
      contigs = contigs.collect{|s| s.gsub(":","__")}
      intervals = intervals.sort do |a,b| 
        astr = a.values_at(0,1).collect{|s| s.to_s.gsub(":","__")} * ":"
        bstr = b.values_at(0,1).collect{|s| s.to_s.gsub(":","__")} * ":"
        Misc.genomic_location_cmp_contigs(astr, bstr, contigs)
      end
    end

    gap = GATKShard::GAP_SIZE * 4
    chunks = []

    current = []
    current_size = 0
    
    last_eend = 0
    last_chr = nil
    if ! break_interval
      TSV.traverse intervals, :type => :array do |chr,start,eend|
        interval_size = eend - start + 1
        last_chr = chr if last_chr.nil?

        if current_size >= chunk_size && (chr != last_chr || (start > last_eend + gap))
          chunks << current
          current = []
          current_size = 0
        end
        current << [chr, start, eend]
        current_size += interval_size

        last_eend = eend
        last_chr = chr
      end
    else
      breaks = 0

      TSV.traverse intervals, :type => :array do |chr,start,eend|
        interval_size = eend - start + 1
        last_chr = chr if last_chr.nil?

        if current_size + interval_size > chunk_size * 1.1
          if current.any? && (current_size + interval_size > chunk_size * 1.1) && current_size > chunk_size * 0.9
            chunks << current
            current = []
            current_size = 0
          end

          remaining = interval_size

          while remaining > 0
            partial_size = [remaining, chunk_size - current_size].min
            breaks += 1
            Log.low "Breaking interval #{[chr, start, eend]} at #{start + partial_size - 1}" if partial_size != remaining
            current << [chr, start, start + partial_size - 1]
            current_size += partial_size
            remaining = remaining - partial_size
            start = start + partial_size
            if current_size > chunk_size * 0.8
              chunks << current
              current = []
              current_size = 0
            end
          end

        elsif current_size + interval_size < chunk_size || (chr == last_chr && (start < last_eend + gap))
          current << [chr, start, eend]
          current_size += interval_size
          next
        else
          current << [chr, start, eend]
          chunks << current
          current = []
          current_size = 0
        end
      end
      Log.low "Broke #{breaks} intervals"
    end

    chunks << current if current.any?

    chunks
  end

  def self.chunks_overlap?(chunks, gap = GATKShard::GAP_SIZE)
    gap = gap.to_i
    last_chr = nil
    last_eend = nil

    chunks.each do |list|
      return true if last_eend && last_chr == list.first[0] && list.first[1].to_i - gap <= last_eend  + gap
      last_chr = list.last[0]
      last_eend = list.last[2].to_i
    end

    return false
  end

  def self.chunk_sizes(chunks)
    sizes = chunks.collect do |list|
      list.inject(0){|acc,e| acc += e[2].to_i - e[1].to_i}
    end
    TmpFile.with_file(sizes.sort * "\n") do  |file|
      ppp CMD.cmd(:gnuplot, :in => <<-EOF).read
set terminal dumb 
plot '#{file}'
      EOF
    end
    [Misc.min(sizes), Misc.mean(sizes), Misc.max(sizes)]
  end

  def self.cmd(command, args, interval_list, chunk_size = CHUNK_SIZE, cpus = nil, contigs = nil, bar = nil, break_interval = nil, &callback)
    interval_file_field = args.keys.select{|f| f =~ /interval/i and f !~ /padding/ }.first
    padding_field = args.keys.select{|f| f =~ /interval/i and f =~ /padding/ }.first
    output_field = args.keys.select{|f| f =~ /output/i }.first
    cpus ||= Rbbt::Config.get('cpus', 'shard', :default => 3) 

    is_bedfile = Path.is_filename?(interval_list) && interval_list =~ /\.bed$/i

    chunks = GATKShard.chunk_intervals(interval_list, chunk_size, contigs)
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size, contigs, break_interval)      if break_interval && chunks.length < (cpus.to_i * 2) || chunks.length < 20
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 2, contigs, break_interval)  if chunks.length < (cpus.to_i * 2) || chunks.length < 20
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 5, contigs, break_interval)  if chunks.length < (cpus.to_i * 2) || chunks.length < 20
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 10, contigs, break_interval) if chunks.length < (cpus.to_i * 2) || chunks.length < 20
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 25, contigs, break_interval) if chunks.length < (cpus.to_i * 2) || chunks.length < 20
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 50, contigs, break_interval) if chunks.length < (cpus.to_i * 2) || chunks.length < 20

    if padding_field && chunks_overlap?(chunks, args[padding_field] || GATKShard::GAP_SIZE)
      args[padding_field] = "0"
    end 

    bar.max = chunks.length if bar
    bar.init if bar

    cpus = [cpus.to_i, chunks.length].min if cpus

    bar.max_history = cpus * 2
    q = RbbtProcessQueue.new cpus

    q.callback &callback

    TmpFile.with_file do |workdir|
      Open.mkdir workdir

      q.init do |intervals|
        Log.low "GATKShard processing intervals #{Misc.fingerprint intervals}"

        if ! is_bedfile
          intervals = intervals.collect{|chr,start,eend| [chr, start.to_i - 1, eend] }
        end

        name = [intervals.first * "__", intervals.last * "__"] * ","
        output = File.join(workdir, name)

        iargs = {} 
        args.each do |k,v|
          if Array === v
            iv = v.collect{|_v| String === _v ? _v.sub("[OUTPUT]", output) : _v }
          else
            iv = String === v ? v.sub("[OUTPUT]", output) : v
          end
          iargs[k] = iv
        end

        if ENV["_JAVA_OPTIONS"]
          java_options = ENV["_JAVA_OPTIONS"]
          max = java_options.scan(/-Xmx(\d+)m/).first[0]
          java_options = java_options.sub("-Xmx#{max}m", "-Xmx#{(max.to_i / cpus.to_i).to_s}m")
        else
          java_options = ""
        end

        #TmpFile.with_file(intervals.collect{|e| [e[0], e[1], e[2].to_i + 1] * "\t"} * "\n", :extension => 'bed') do |interval_file|
        TmpFile.with_file(intervals.collect{|e| [e[0], e[1], e[2].to_i] * "\t"} * "\n", :extension => 'bed') do |interval_file|
          iargs[interval_file_field] = interval_file 
          iargs[output_field] = output 
          Misc.with_env "_JAVA_OPTIONS", java_options do
              GATK.run_log(command, iargs)
          end
        end
        output
      end

      TSV.traverse chunks, :type =>:array do |intervals|
        q.process intervals
      end

      q.join
    end

    bar.done if bar
  end
end
