require 'tools/control_FREEC'
module HTS

  input :config, :file, "Configuration file", nil, :nofile => true
  extension :vcf
  task :control_freeC => :text do |config|
    output = file('output')

    ControlFREEC.run(config)
    
#    Open.read(output....)
  end
end
