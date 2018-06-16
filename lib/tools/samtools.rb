require 'rbbt-util'
require 'rbbt/resource'

module Samtools
  extend Resource
  self.subdir = 'share/databases/Samtools'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.HTSLib, :install, Rbbt.share.install.software.HTSLib.find
  Rbbt.claim Rbbt.software.opt.Samtools, :install, Rbbt.share.install.software.Samtools.find

  Samtools_CMD=Rbbt.software.opt.Samtools.produce.samtools.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.HTSLib.produce 
  iif Rbbt.software.opt.Samtools.produce 
end
