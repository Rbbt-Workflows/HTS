require 'rbbt-util'
require 'rbbt/resource'

module OptiType
  
  extend Resource
  self.subdir = 'share/databases/OptiType'

  Rbbt.claim Rbbt.software.opt.OptiType, :install, Rbbt.share.install.software.OptiType.find
  Rbbt.claim Rbbt.software.opt.seqan, :install, Rbbt.share.install.software.seqan.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.OptiType.produce 
  iif Rbbt.software.opt.seqan.produce 
end
