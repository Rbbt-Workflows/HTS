require 'rbbt-util'
require 'rbbt/resource'

module HapPy
  extend Resource
  self.subdir = 'share/databases/HapPy'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  Rbbt.claim Rbbt.software.opt.HapPy, :install do
    url = "https://github.com/Illumina/hap.py.git"
    commands =<<-EOF
get_git "$name" "$url"
cd "$OPT_BUILD_DIR/$name"
unset BOOST_ROOT
python2 "install.py" "`opt_dir $name`"
make install
setup "$name"
    EOF
    {:url => url, :commands => commands}
  end

  CMD.tool "hap.py", Rbbt.software.opt.HapPy
end

Log.severity = 0 if __FILE__ == $0
iif Rbbt.software.opt.HapPy.produce if __FILE__ == $0

