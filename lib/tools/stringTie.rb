require 'rbbt-util'
require 'rbbt/resource'

module StringTie
  extend Resource
  self.subdir = 'share/databases/StringTie'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Rbbt.claim Rbbt.software.opt.StringTie, :install, Rbbt.share.install.software.StringTie.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.StringTie.produce
end
