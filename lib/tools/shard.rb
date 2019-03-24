class GATKShard

  def self.process_intervals(interval_list=nil, chunk_size=10_000_000, &block)
    current = []
    current_size = 0
    TSV.traverse interval_list, :type => :array do |line|
      chr, start, eend, *rest = line.split("\t")

      start = start.to_i
      eend = eend.to_i
      remaining = eend - start

      while remaining > 0
        size = [chunk_size, remaining].min

        current << [chr, start.to_s, (start + size).to_s ]
        current_size += size

        if current_size >= chunk_size
          yield current
          current = []
          current_size = 0
        end

        remaining = remaining - size
        start += size
      end
    end
    yield current if current.any?
  end


  def self.cmd(command, args, interval_list, chunk_size = 10_000_000, cpus = nil, &callback)
    interval_file_field = args.keys.select{|f| f =~ /interval/i and f !~ /padding/ }.first
    output_field = args.keys.select{|f| f =~ /output/i }.first
    cpus = Rbbt::Config.get('cpus', 'shard', :default => 3) 

    q = RbbtProcessQueue.new cpus

    q.callback &callback

    TmpFile.with_file do |workdir|
      Open.mkdir workdir

      q.init do |intervals|
        Log.low "GATKShard processing intervals #{Misc.fingerprint intervals}"
        iargs = args.dup
        output = File.join(workdir, Misc.obj2digest(intervals))
        TmpFile.with_file(intervals.collect{|e| e * "\t"} * "\n", :extension => 'bed') do |interval_file|
          iargs[interval_file_field] = interval_file 
          iargs[output_field] = output 
          GATK.run_log(command, iargs)
        end
        output
      end

      GATKShard.process_intervals(interval_list, chunk_size) do |intervals|
        q.process intervals
      end
      q.join
    end
  end

end
