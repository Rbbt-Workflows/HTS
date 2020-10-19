class GATKShard

  GAP_SIZE = 1_000
  CHUNK_SIZE = 10_000_000

  def self.chunk_intervals(interval_list=nil, chunk_size=10_000_000, contigs = nil, break_interval = false)
    intervals = TSV.traverse interval_list, :type => :array, :into => [] do |line|
      next if line =~ /^@/ 
      chr, start, eend, *rest = line.split("\t")
      next unless start =~ /^\d+$/
      next unless eend =~ /^\d+$/
      [chr, start.to_i, eend.to_i]
    end

    if contigs
      # Deal with contigs with ":" characters
      contigs = contigs.collect{|s| s.gsub(":","__")}
      intervals = intervals.sort do |a,b| 
        astr = a[0,1].collect{|s| s.gsub(":","__")} * ":"
        bstr = b[0,1].collect{|s| s.gsub(":","__")} * ":"
        Misc.genomic_location_cmp_contigs(astr, bstr, contigs)
      end
    end

    chunks = []
    current = []

    current_size = 0
    last_eend = 0
    last_chr = nil
    gap = GATKShard::GAP_SIZE * 4
    TSV.traverse intervals, :type => :array do |chr,start,eend|
      remaining = eend - start + 1
      last_chr = chr if last_chr.nil?

      while remaining > 0
        if break_interval
          size = [chunk_size, remaining].min
        else
          size = remaining
        end

        if current_size >= chunk_size && (chr != last_chr || (start > last_eend + gap))
          chunks << current
          current = []
          current_size = 0
        end

        current << [chr, start.to_s, (start + size - 1).to_s ]
        current_size += size

        remaining = remaining - size
        start += size
      end

      last_eend = eend
      last_chr = chr
    end
    chunks << current if current.any?

    chunks
  end


  def self.cmd(command, args, interval_list, chunk_size = 10_000_000, cpus = nil, contigs = nil, bar = nil, break_interval = false, &callback)
    interval_file_field = args.keys.select{|f| f =~ /interval/i and f !~ /padding/ }.first
    output_field = args.keys.select{|f| f =~ /output/i }.first
    cpus ||= Rbbt::Config.get('cpus', 'shard', :default => 3) 

    q = RbbtProcessQueue.new cpus

    q.callback &callback

    chunks = GATKShard.chunk_intervals(interval_list, chunk_size, contigs)
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 5, contigs, break_interval) if chunks.length < (cpus.to_i * 2) || chunks.length < 25
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 20, contigs, break_interval) if chunks.length < (cpus.to_i * 2) || chunks.length < 25

    bar.max = chunks.length if bar
    bar.init if bar

    TmpFile.with_file do |workdir|
      Open.mkdir workdir

      q.init do |intervals|
        Log.low "GATKShard processing intervals #{Misc.fingerprint intervals}"
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

        #TmpFile.with_file(intervals.collect{|e| [e[0], e[1], e[2].to_i + 1] * "\t"} * "\n", :extension => 'bed') do |interval_file|
        TmpFile.with_file(intervals.collect{|e| [e[0], e[1], e[2].to_i] * "\t"} * "\n", :extension => 'bed') do |interval_file|
          iargs[interval_file_field] = interval_file 
          iargs[output_field] = output 
          GATK.run_log(command, iargs)
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
