require 'rbbt-util'
require 'rbbt/workflow'
module FuSeq
  # ToDo pull files out of directory
  
  Rbbt.claim Rbbt.software.opt.FuSeq, :install, Rbbt.share.install.software.FuSeq.find

  def self.index_FASTA(file, kmer = nil, dir = nil)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.fasta_TxIndices[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".idx") || Persist.newer?(linked + '.idx', file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        args = {}
        gunzip_cmd = "gunzip -c \"#{linked}\"" +  '| sed "s/\(EN[[:alpha:]]*[[:digit:]]*\)\.[[:digit:]]/\1/g"'
        cmd = "TxIndexer -t <(#{gunzip_cmd})"
        cmd += " -k 21 " if kmer
        cmd += " -o \"#{linked}.idx\""
        cmd = "bash -c '#{cmd}'"
        CMD.cmd_log(cmd)
        #if kmer
        #  CMD.cmd_log("bash -c 'TxIndexer -t <(gunzip -c \"#{linked}\"" + '|sed \'s/\(EN[[:alpha:]]*[[:digit:]]*\)\.[[:digit:]]/\1/g\')' + " -k 21 -o \"#{linked}.idx\"'")
        #else
        #  CMD.cmd_log("bash -c 'TxIndexer -t <(gunzip -c \"#{linked}\"" + '|sed \'s/\(EN[[:alpha:]]*[[:digit:]]*\)\.[[:digit:]]/\1/g\')' + " -o \"#{linked}.idx\"'")
        #  #CMD.cmd_log("bash -c 'TxIndexer -t <(gunzip -c \"#{linked}\" | sed 's/\(EN[[:alpha:]]*[[:digit:]]*\)\.[[:digit:]]/\1/g') -o \"#{linked}.idx\"'")
        #end
      end
    end

    linked + '.idx'
  end

  def self.run(read_files, organism, kmer, cpus, output)
    Rbbt.software.opt.FuSeq.produce

    gtf_file = Organism.gene_set(organism).find
    fasta = Organism.cdna_fasta(organism).find

    index = self.index_FASTA(fasta, kmer)

    reads = read_files.collect{|f| Open.gzip?(f) ? "<(gunzip -c '#{ f }')" : "'#{ f }'" }
    gtf_file = Open.gzip?(gtf_file) ? "<(gunzip -c '#{ gtf_file }')" : "'#{ gtf_file }'" 

    cpus = Rbbt::Config.get(:cpus, :FuSeq, :fuseq, :default => 1) if cpus.nil?
    
    if reads.length == 1
      CMD.cmd_log("bash -c 'FuSeq -i \"#{index}\" -l IU -1 #{reads.first} -p \"#{cpus}\" -g #{gtf_file} -o \"#{output}\"'")
    else
      CMD.cmd_log("bash -c 'FuSeq -i \"#{index}\" -l IU -1 #{reads.first} -2 #{reads.last} -p \"#{cpus}\" -g #{gtf_file} -o \"#{output}\"'")
    end
  end

  def self.prepare_sqlite(file)
    file = file.path if Step === file
    file = file.find if Path === file
    file = File.expand_path(file)


    digest = Misc.file2md5(file)
    basename = File.basename(file)

    dir = Rbbt.var.FuSeq_sqlite[digest].find if dir.nil?
    Path.setup(dir) unless Path === dir

    linked = dir[basename].find
    if ! File.exists?(linked + ".sql") || Persist.newer?(linked + '.sql', file)

      Misc.in_dir dir do
        FileUtils.ln_s file, dir[basename] unless File.exists?(linked)
        args = {}
        CMD.cmd_log("Rscript #{Rbbt.software.opt.FuSeq.find}/R/createSqlite.R #{file} #{linked}.sql")
      end
    end

    linked + '.sql'
  end

  def self.process(results, organism, output)
    gtf_file = Organism.gene_set(organism).find
    fasta = Organism.cdna_fasta(organism).find

    sqlite = FuSeq.prepare_sqlite(gtf_file)

    grc_build = Organism.GRC_build organism

    annot = Rbbt.software.opt.FuSeq.data.glob("*" << grc_build << "*").first

    Open.mkdir output
    params = Rbbt.software.opt.FuSeq.R["params.txt"].find
    CMD.cmd_log("Rscript #{Rbbt.software.opt.FuSeq.find}/R/FuSeq.R in='#{results}' txfasta='#{fasta}' sqlite='#{sqlite}' txanno='#{annot}' out='#{output}' params='#{params}'")
  end

end

if __FILE__ == $0
  Log.severity = 0 
  iii Rbbt.software.opt.FuSeq.produce(true)
end


