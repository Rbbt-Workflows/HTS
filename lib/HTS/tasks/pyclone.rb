Workflow.require_workflow "Sequence"

module HTS
  CMD.tool "pyclone-vi" do
    'pip install git+https://github.com/Roth-Lab/pyclone-vi.git'
  end

  input :copynumber, :tsv, "Copy number segments by sequenza"
  input :sample_name, :string, "Sample name"
  input :tumor_content, :float, "Tumor content", 1
  dep HTS, :mutation_pileup
  dep Sequence, :expanded_vcf, :compute => :produce
  task :pyclone_variant_info => :tsv do |copynumber,sample_name, tumor_content|
    dumper = TSV::Dumper.new :key_field => "mutation_id", :fields => ["sample_id", "ref_counts", "alt_counts", "normal_cn", "minor_cn", "major_cn", "tumor_content"], :type => :list
    dumper.init(:preamble => "", :header_hash => '')

    sample_name ||= clean_name

    copynumber.add_field "Start" do |k,v|
      k.split(":")[1].to_i
    end

    copynumber.add_field "End" do |k,v|
      k.split(":")[2].to_i
    end

    vcf_info = step(:expanded_vcf).load.to_list
    af_field = vcf_info.fields.select{|f| ! f.include?("normal") && f.split(":").last == "AD" }.first
    af_pos = vcf_info.fields.index af_field
    cn_index = copynumber.range_index("Start", "End")
    TSV.traverse step(:mutation_pileup), :type => :array, :into => dumper do |line|
      chr, pos, ref, count, bases, *rest = line.split("\t")

      if vcf_info[mutation]
        ref_counts, alt_counts = vcf_info[mutation][af_pos].split(",")
      else
        base_counts = Misc.counts(bases.upcase.chars.select{|c| %w(A C T G).include? c })
        chr, pos, alt = mutation.split(":")

        alt_counts = base_counts[alt.upcase] 
        ref_counts = (base_counts.values - [alt_counts]).max
      end
      cn_key = cn_index[pos.to_i].select{|key| key.split(":").first.sub('chr', '') == chr.sub('chr', '')}.first

      if cn_key 
        values = copynumber[cn_key]
        major_cn = values[7]
        minor_cn = values[8]
      else
        major_cn = 1
        minor_cn = 1
      end

      if %w(y chry).include? chr.downcase
        normal_cn = 1
      else
        normal_cn = 2
      end

      [mutation, [sample_name, ref_counts || 1, alt_counts || 1, normal_cn, minor_cn, major_cn, tumor_content]]
    end
  end

  dep :pyclone_variant_info
  dep Sequence, :affected_genes
  task :cancer_cell_fraction => :tsv do
    affected_genes = step(:affected_genes).load

    tsv = step(:pyclone_variant_info).join.path.tsv :header_hash => "", :type => :list
    tsv.add_field "Cancer cell fraction" do |mutation,values|
      sample, ref, alt, normal, minor, major, tumor_content = values

      tumor_content = 1 if tumor_content.nil? || tumor_content.empty?
      p = tumor_content.to_f

      cn = normal.to_i
      ct = minor.to_i + major.to_i

      vaf = alt.to_f / (ref.to_i + alt.to_i)
      ccf = vaf / p * (p * ct + (1.0-p) * cn)

      ccf
    end

    organism = TSV.get_stream(self.recursive_inputs[:organism]).read
    ensg2name = Organism.identifiers(organism).index :target => "Associated Gene Name", :persist => true
    tsv.add_field "Associated Gene Name" do |mutation,values|
      ensg2name.values_at(*affected_genes[mutation]).compact
    end
  end

  # VAF / p * (p * c_t + (1-p) * c_n)
  # 
  # where VAF is the alternate allele count divided by the total
  # (alternate+reference) allele count, p is the purity and c_t (c_n) is the
  # copy number or ploidy of the tumor (normal) at the mutated site.
  # Alternatively, some of the phylogenetic reconstruction algorithms compute
  # CCF as well.



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
