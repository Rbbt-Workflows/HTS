require 'rbbt-util'

module Qualimap
  Rbbt.claim Rbbt.software.opt.Qualimap, :install do
    url = "https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip"
    commands =<<-EOF
install_src "$name" "$url"
rm $OPT_BIN_DIR/qualimap
ln -s $(opt_dir $name)/qualimap $OPT_BIN_DIR/qualimap
    EOF
    {:url => url, :commands => commands}
  end

  CMD.tool :qualimap, Rbbt.software.opt.Qualimap
end

if __FILE__ == $0

  Rbbt.software.opt.Qualimap.produce(true)
end
