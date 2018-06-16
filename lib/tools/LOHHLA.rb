require 'rbbt-util'
require 'rbbt/resource'

module LOHHLA
  
  extend Resource
  self.subdir = 'share/databases/LOHHLA'

  Rbbt.claim Rbbt.software.opt.Jellyfish, :install, Rbbt.share.install.software.Jellyfish.find
  Rbbt.claim Rbbt.software.opt.NovoAlign, :install, Rbbt.share.install.software.NovoAlign.find
  Rbbt.claim Rbbt.software.opt.LOHHLA, :install, Rbbt.share.install.software.LOHHLA.find
  Rbbt.claim Rbbt.software.opt.PicardTools, :install, Rbbt.share.install.software.PicardTools.find
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.NovoAlign.produce(true)
end
