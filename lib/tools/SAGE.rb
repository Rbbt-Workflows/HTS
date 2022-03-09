require 'rbbt-util'
module SAGE

  Rbbt.claim Rbbt.software.opt.SAGE, :install, {:jar => "https://github.com/hartwigmedical/hmftools/releases/download/sage-v2.8/sage-2.8.jar"}

  Rbbt.claim Rbbt.share.databases.SAGE, :proc do |dirname|
    TmpFile.with_file(Open.open("https://nextcloud.hartwigmedicalfoundation.nl/s/LTiKTd8XxBqwaiC/download?path=%2FHMFTools-Resources&files=Sage&downloadStartSecret=gi9gfliagog", :mode => 'rb')) do |zip|
      TmpFile.with_file do |tmpdir|
        Misc.unzip_in_dir(zip, tmpdir)
        Open.mv File.join(tmpdir, 'Sage'), dirname
      end
    end
  end

  def self.run(tumor, reference, output, ref_version = '38', options = {})
    Rbbt.software.opt.SAGE.produce

    version = ref_version.to_s.sub('hg', '')
    version = '37' if version == '19'
    version = '37' if version == 'b37'

    case version
    when '38'
      reference_fa = GATK.prepare_FASTA(HTS.helpers[:reference_file].call 'hg38')
      assembly = 'hg38'
    when '37', '19'
      reference_fa = GATK.prepare_FASTA(HTS.helpers[:reference_file].call 'b37')
      assembly = 'hg19'
    else
      raise "Unkown reference version #{Misc.fingerprint ref_version}"
    end

    reference_fa = reference_fa.sub(/\.gz$/,'')

    type = options[:type] || 'somatic'

    Rbbt.share.databases.SAGE.produce
    hotspots = Rbbt.share.databases.SAGE[version]["KnownHotspots.#{type}.#{version}.vcf.gz"].find

    reference_name = GATK.BAM_sample_name(reference)
    tumor_name = GATK.BAM_sample_name(tumor)

    cpus =  options[:cpus] || 1

    cmd =<<-EOF.split("\n") * " "
java -Xms4G -Xmx32G -cp #{Rbbt.software.opt.jars["SAGE.jar"].find} com.hartwig.hmftools.sage.SageApplication
  -bqr_plot false
  -threads #{cpus} 
  -reference #{reference_name} -reference_bam #{reference}
  -tumor #{tumor_name} -tumor_bam #{tumor}
  -ref_genome #{reference_fa}
  -assembly #{assembly}
  -out #{output}
    EOF

    cmd << " -hotspots #{hotspots}" if hotspots
    #cmd << " -high_confidence_bed #{high_confidence_bed}" if high_confidence_bed
    #cmd << " -ensembl_data_dir #{ensembl_cache}" if ensembl_cache

    CMD.cmd_log(cmd)
  end

end

if __FILE__ == $0
  Log.with_severity 0 do
    Rbbt.software.opt.SAGE.produce
    Rbbt.share.databases.SAGE.produce(true)
  end
end
