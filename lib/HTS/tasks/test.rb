module HTS

  input :normal, :file, "Normal BAM", nil, :nofile => true
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38), :nofile => true
  dep :pileup, :bam => :placeholder do |jobname,options|
    deps = []

    deps << {:inputs => options.merge(:bam => options[:normal]), :jobname => jobname + '.normal'}
    deps << {:inputs => options.merge(:bam => options[:tumor]), :jobname => jobname + '.tumor'}

    deps
  end
  dep :GC_windows, :jobname => :reference
  input :normal_purity, :float, "Normal sample purity", 1 
  input :tumor_purity, :float, "Tumor sample purity", 1 
  extension 'vcf'
  task :speed_test => :text do |normal,tumor,reference,normal_purity,tumor_purity|
    Open.mkdir files_dir
    normal = dependencies.select{|dep| dep.clean_name.include? '.normal'}.first
    tumor = dependencies.select{|dep| dep.clean_name.include? '.tumor'}.first
    output = file('output')
    Misc.in_dir output do
      io_normal = CMD.cmd("zcat |sed 's/\t/#/;s/\t/#/' ", :pipe => true, :in => TSV.get_stream(normal, :noz => true))
      io_tumor = CMD.cmd("zcat |sed 's/\t/#/;s/\t/#/' ", :pipe => true, :in => TSV.get_stream(tumor, :noz => true))

      pipe = Misc.paste_streams([io_normal, io_tumor]) do |a,b|
        chr1, _sep, pos1 = a.partition("#")
        chr2, _sep, pos2 = b.partition("#")

        #chr1, pos1 = a.split("#")
        #chr2, pos2 = b.split("#")

        puts [chr1,chr2] * " "
        cmp = Misc.chr_cmp_strict(chr1, chr2)
        case cmp
        when -1
          -1
        when 1
          1
        else
          pos1.to_i <=> pos2.to_i
        end
      end

      TSV.traverse pipe, :type => :line, :bar => true do |line|
      end
    end
    nil
  end


  dep :strelka
  dep :mutect2_clean
  task :compare_variants => :tsv do

    strelka_variants = TSV.traverse step(:strelka), :into => [] do |chr, values|
      pos, id, ref, alt, *rest = values

      [chr, pos, alt] * ":"
    end

    mutect2_variants = TSV.traverse step(:mutect2_clean), :into => [] do |chr, values|
      pos, id, ref, alt, *rest = values

      [chr, pos, alt] * ":"
    end

    common = strelka_variants & mutect2_variants
    only_strelka = strelka_variants - mutect2_variants
    only_mutect2 =  mutect2_variants - strelka_variants

    
    tsv = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single, :cast => :to_i)

    tsv["Common"] = common.length
    tsv["Only Strelka"] = only_strelka.length
    tsv["Only MuTect2"] = only_mutect2.length

    tsv
  end

end
