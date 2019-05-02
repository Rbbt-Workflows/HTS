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
      intervals = intervals.sort{|a,b| Misc.genomic_location_cmp_contigs(a[0,1] * ":", b[0,1] * ":", contigs)}
    end

    chunks = []
    current = []

    current_size = 0
    last_eend = 0
    last_chr = nil
    gap = GATKShard::GAP_SIZE * 4
    TSV.traverse intervals, :type => :array do |chr,start,eend|
      remaining = eend - start
      last_chr = chr if last_chr.nil?

      while remaining > 0
        if break_interval
          size = [chunk_size, remaining].min
        else
          size = remaining
        end

        if chr != last_chr || (current_size >= chunk_size && start > last_eend + gap)
          chunks << current
          current = []
          current_size = 0
        end

        current << [chr, start.to_s, (start + size).to_s ]
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


  def self.cmd(command, args, interval_list, chunk_size = 10_000_000, cpus = nil, contigs = nil, bar = nil, &callback)
    interval_file_field = args.keys.select{|f| f =~ /interval/i and f !~ /padding/ }.first
    output_field = args.keys.select{|f| f =~ /output/i }.first
    cpus ||= Rbbt::Config.get('cpus', 'shard', :default => 3) 

    q = RbbtProcessQueue.new cpus

    q.callback &callback

    chunks = GATKShard.chunk_intervals(interval_list, chunk_size, contigs)
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 10, contigs) if chunks.length < (cpus.to_i * 2) || chunks.length < 25
    chunks = GATKShard.chunk_intervals(interval_list, chunk_size / 100, contigs) if chunks.length < (cpus.to_i * 2) || chunks.length < 25
    bar.max = chunks.length if bar
    bar.init if bar

    TmpFile.with_file do |workdir|
      Open.mkdir workdir

      q.init do |intervals|
        Log.low "GATKShard processing intervals #{Misc.fingerprint intervals}"
        iargs = args.dup
        name = [intervals.first * "__", intervals.last * "__"] * ","
        output = File.join(workdir, name)
        TmpFile.with_file(intervals.collect{|e| e * "\t"} * "\n", :extension => 'bed') do |interval_file|
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
