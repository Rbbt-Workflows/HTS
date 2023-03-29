class Scatter
  def initialize(bam1, bam2, directory, chunk_size=10_000_000, padding=10_000, &block)
    @directory = directory.dup

    Path.setup @directory unless Path === @directory

    if Misc.is_filename? bam1
      @stream1 = CMD.cmd("samtools view --no-PG -h '#{bam1}'", :pipe => true, :no_fail => true)
    else
      @stream1 = TSV.get_stream bam1
    end

    if Misc.is_filename? bam2
      @stream2 = CMD.cmd("samtools view --no-PG -h '#{bam2}'", :pipe => true, :no_fail => true)
    else
      @stream2 = TSV.get_stream bam2
    end

    @chunk_size = chunk_size
    @padding = padding
    @block = block

    @padding_str1 = ""
    @padding_str2 = ""
  end

  def read_head(comment_char = "@")
    @head1 = ""
    while @line1 = @stream1.gets
      break unless @line1[0] == comment_char
      @head1 << @line1
    end

    @head2 = ""
    while @line2 = @stream2.gets
      break unless @line2[0] == comment_char
      @head2 << @line2
    end
  end

  def self.read_stream_chunk(stream, out, chunk_size)
    total = 0
    closed = false
    while total < chunk_size && ! closed
      chunk = nil
      begin
        chunk = stream.readpartial(Misc::BLOCK_SIZE)
      rescue EOFError
        closed = true
      end
      total += chunk.length unless chunk.nil?
      out.write chunk unless chunk.nil?
    end
    return :done if closed

    line = stream.gets
    return :done if line.nil?
    out.write line

    line = stream.gets
    return :done if line.nil?
    out.write line

    current_chr, current_pos = line.split("\t").values_at 2, 3
    current_pos = current_pos.to_i

    [current_chr, current_pos]
  end

  def self.read_stream_padding(stream, out, min_chr, min_pos)
    do_padding = true
    padding_str = ""
    while do_padding
      line = stream.gets
      break if line.nil?
      out.write line
      padding_str << line

      chr, pos = line.split("\t").values_at 2, 3
      pos = pos.to_i
      do_padding = false if Misc.chr_cmp_strict(chr, min_chr) == 1 || (Misc.chr_cmp_strict(chr, min_chr) == 0 && (min_pos != 0 && pos > min_pos))
    end

    padding_str
  end

  def prepare_files(name)
    file1 = @directory[name + ".1"]
    file2 = @directory[name + ".2"]

    read_head if @head1.nil?
    Open.write(file1) do |out1|
      Open.write(file2) do |out2|
        out1.write @head1
        out2.write @head2

        out1.write @line1  if @line1
        out2.write @line2 if @line2

        @line1, @line2 = nil, nil

        out1.write @padding_str1
        out2.write @padding_str2

        current_chr1, current_pos1 = Scatter.read_stream_chunk(@stream1, out1, @chunk_size)
        current_chr2, current_pos2 = Scatter.read_stream_chunk(@stream2, out2, @chunk_size)

        current_chr1, current_pos1 = "ZZZ", 0 if current_chr1 == :done
        current_chr2, current_pos2 = "ZZZ", 0 if current_chr2 == :done

        case Misc.chr_cmp_strict(current_chr1, current_chr2)
        when 0
          min_chr = current_chr1
          min_pos = [current_pos1, current_pos2].max
        when -1
          min_chr = current_chr2
          min_pos = current_pos2
        when 1
          min_chr = current_chr1
          min_pos = current_pos1
        end

        @padding_str1 = ""
        @padding_str2 = ""

        @padding_str1 = Scatter.read_stream_padding(@stream1, out1, min_chr, min_pos + @padding) unless current_chr1 == :done
        @padding_str2 = Scatter.read_stream_padding(@stream2, out2, min_chr, min_pos + @padding) unless current_chr2 == :done

      end
    end

    [file1, file2]
  end

  def iterate(&block)
    @chunk = 0
    while ! @stream1.eof? && ! @stream2.eof?
      name = "chunk-" + @chunk.to_s
      file1, file2 = prepare_files(name)
      yield file1, file2
      @chunk += 1
    end
  end

  def process(cpus, &block)
    q = RbbtProcessQueue.new cpus

    lock = @directory.lock.find

    res = []
    q.callback do |res_file|
      res << res_file
    end

    q.init do |file1, file2, lock|
      Open.rm lock if File.exist? lock
      block.call file1, file2
    end

    times = 500

    iterate do |file1, file2|
      sleep 2 while File.exist?(lock)
      FileUtils.touch lock
      q.process file1, file2, lock
    end

    q.join

    res
  end
end
