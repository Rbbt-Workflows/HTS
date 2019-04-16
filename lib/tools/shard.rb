class GATKShard

  def self.chunk_intervals(interval_list=nil, chunk_size=10_000_000)
    current = []
    current_size = 0
    chunks = []
    TSV.traverse interval_list, :type => :array do |line|
      next if line =~ /^@/
      chr, start, eend, *rest = line.split("\t")

      start = start.to_i
      eend = eend.to_i
      remaining = eend - start

      while remaining > 0
        size = [chunk_size, remaining].min

        current << [chr, start.to_s, (start + size).to_s ]
        current_size += size

        if current_size >= chunk_size
          chunks << current
          current = []
          current_size = 0
        end

        remaining = remaining - size
        start += size
      end
    end
    chunks << current if current.any?

    chunks
  end


  def self.cmd(command, args, interval_list, chunk_size = 10_000_000, cpus = nil, bar = nil, &callback)
    interval_file_field = args.keys.select{|f| f =~ /interval/i and f !~ /padding/ }.first
    output_field = args.keys.select{|f| f =~ /output/i }.first
    cpus ||= Rbbt::Config.get('cpus', 'shard', :default => 3) 

    q = RbbtProcessQueue.new cpus

    q.callback &callback

    chunks = GATKShard.chunk_intervals(interval_list, chunk_size)
    bar.max = chunks.length if bar
    bar.init if bar

    TmpFile.with_file do |workdir|
      Open.mkdir workdir

      q.init do |intervals|
        Log.low "GATKShard processing intervals #{Misc.fingerprint intervals}"
        iargs = args.dup
        output = File.join(workdir, intervals.first * "__")
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
