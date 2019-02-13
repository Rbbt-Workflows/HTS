require 'rbbt-util'
require 'rbbt/resource'

module HISAT
  extend Resource
  self.subdir = 'share/databases/HISAT'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Rbbt.claim Rbbt.software.opt.HISAT, :install, Rbbt.share.install.software.HISAT.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HISAT.produce
end
