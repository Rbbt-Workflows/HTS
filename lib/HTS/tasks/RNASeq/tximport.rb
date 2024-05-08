module HTS

  input :quantification_method, :select, "Quantification method to use", :salmon, :select_options => %w(salmon kallisto)
  dep :salmon do |jobname,options|
    method = options[:quantification_method]
    {:task => method, :inputs => options}
  end
  task :tximport => :tsv do |method|
    require 'rbbt/util/R'

    quant = dependencies.first
    organism = quant.inputs[:organism]
    tx2gene = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Gene ID"], :type => :single
    
    tx2gene_file = file('tx2g.tsv')
    data_file = file('data.tsv')
    Open.write(tx2gene_file, tx2gene.to_s(:preamble => false, :header_hash => ''))
    Open.write(data_file, quant.path.tsv(:fix => proc{|l| l.sub(/^(ENS.\d+)\.\d+/, '\1') }).to_s(:preamble => false, :header_hash => ''))
    R.run <<-EOF
rbbt.require('tximport')
rbbt.require('readr')
tx2gene = readr::read_tsv("#{tx2gene_file}")
txi <- tximport("#{data_file}", type = "#{method}", tx2gene = tx2gene)

data = NULL
for(name in names(txi)){
  if (is.array(txi[[name]])){
    if (is.null(data)){
      data = data.frame(txi[[name]])
      names(data) <- c(name)
      rownames(data) <- rownames(txi[[name]])
    } else {
      data[[name]] = txi[[name]]
    }
  }
}

rbbt.tsv.write('#{self.tmp_path}', data, key.field = "Ensembl Gene ID", extra_headers=":type=:list#:namespace=#{organism}#:cast=:to_f")
NULL
    EOF
    nil
  end

  dep_task :gene_expression, HTS, :tximport
end
