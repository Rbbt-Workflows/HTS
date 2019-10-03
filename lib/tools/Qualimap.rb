require 'rbbt-util'

module Qualimap
  Rbbt.claim Rbbt.software.opt.Qualimap, :install, "https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip"

end

if __FILE__ == $0

  Rbbt.software.opt.Qualimap.produce
end
