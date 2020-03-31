require 'rbbt-util'

module FastQC
  Rbbt.claim Rbbt.software.opt.FastQC, :install do
    url = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
    commands =<<-EOF
install_src "$name" "$url"
chmod +x "$(opt_dir $name)/fastqc"
setup "$name"
    EOF
    {:url => url, :commands => commands}
  end

  CMD.tool :fastqc, Rbbt.software.opt.FastQC
end

if __FILE__ == $0
  Rbbt.software.opt.FastQC.produce
end
