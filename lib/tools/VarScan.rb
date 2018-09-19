require 'rbbt-util'
require 'rbbt/resource'

module VarScan
  extend Resource
  self.subdir = 'share/databases/VarScan'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Rbbt.claim Rbbt.software.opt.VarScan, :install, Rbbt.share.install.software.VarScan.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.VarScan.produce
end
