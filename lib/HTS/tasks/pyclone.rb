module HTS
  CMD.tool "pyclone-vi" do
    'pip install git+https://github.com/Roth-Lab/pyclone-vi.git'
  end


  input :copynumber, :tsv, "Copy number segments by sequenza"
  input :sample_name, :string, "Sample name"
  dep HTS, :mutation_pileup
  task :pyclone_variant_info => :tsv do |copynumber,sample_name|
    dumper = TSV::Dumper.new :key_field => "mutation_id", :fields => ["sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn"], :type => :list
    dumper.init(:preamble => "", :header_hash => '')

    sample_name ||= clean_name

    copynumber.add_field "Start" do |k,v|
      k.split(":")[1].to_i
    end

    copynumber.add_field "End" do |k,v|
      k.split(":")[2].to_i
    end

    cn_index = copynumber.range_index("Start", "End")
    TSV.traverse step(:mutation_pileup), :type => :array, :into => dumper do |line|
      chr, pos, ref, count, bases, *rest = line.split("\t")

      mutation = rest.pop

      base_counts = Misc.counts(bases.upcase.chars.select{|c| %w(A C T G).include? c })
      chr, pos, alt = mutation.split(":")

      alt_counts = base_counts[alt.upcase] 
      ref_counts = (base_counts.values - [alt_counts]).max
      cn_key = cn_index[pos.to_i].select{|key| key.split(":").first.sub('chr', '') == chr.sub('chr', '')}.first

      if cn_key 
        values = copynumber[cn_key]
        major_cn = values[7]
        minor_cn = values[8]
      else
        major_cn = 1
        minor_cn = 1
      end

      normal_cn = 2

      [mutation, [sample_name, ref_counts || 1, alt_counts || 1, major_cn, minor_cn, normal_cn]]
    end
  end


  dep :pyclone_variant_info
  input :clusters, :integer, "Number of clusters", 10
  input :restarts, :integer, "Number of restarts", 10
  input :distribution, :select, "Distribution to use", "binomial", :select_options => %w(binomial beta-binomial)
  input :min_alt_reads, :integer, "Minimum number of alt reads in at least a sample to consider mutations", 5
  task :pyclone => :tsv do |clusters,restarts,distribution,min_alt_reads|
    info = file('pyclone_info.tsv')
    good_mutations = step(:pyclone_variant_info).join.path.tsv(:header_hash => '').select("alt_counts"){|list| list.select{|v| v.to_i >= min_alt_reads}.any? }.keys
    io = TSV.traverse step(:pyclone_variant_info), :into => :stream, :type => :array do |line|
      next line if line =~ /mutation_id/
      next unless good_mutations.include? line.partition("\t").first
      line
    end
    Open.write(info, io)

    h5 = file('result.h5')
    result = file('result.tsv')
    CMD.cmd_log("pyclone-vi", "fit -i #{info} -o #{h5} -c #{clusters} -d #{distribution} -r #{restarts}")
    CMD.cmd_log("pyclone-vi", "write-results-file -i #{h5} -o #{result}")
    TSV.open(result, :header_hash => '', :type => :list)
  end
end
