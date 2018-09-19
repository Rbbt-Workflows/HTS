require 'rbbt-util'
require 'rbbt/resource'

module AlleleCount
  extend Resource
  self.subdir = 'share/databases/AlleleCount'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.AlleleCount, :install, Rbbt.share.install.software.AlleleCount.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.AlleleCount.produce(true)
end
