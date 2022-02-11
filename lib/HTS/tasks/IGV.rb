require 'tools/IGV'
module HTS

  input :positions, :array, "Genomic Positions"
  input :tumor, :file, "Tumor BAM", nil, :nofile => true
  input :normal, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :depth, :integer, "Expected read depth", 50
  task :mutation_BAM_img => :array do |positions,tumor,normal,reference,depth|
    output = file('output')
    orig_reference = reference

    reference = reference_file reference

    tumor = Samtools.prepare_BAM(tumor) if tumor
    normal = Samtools.prepare_BAM(normal) if normal

    reference = Samtools.prepare_FASTA(reference)

    #reference = reference.sub('.gz', '')

    Open.mkdir files_dir

    panel_height = depth * 7 * 2
    panel_height *= 2 if normal

    margin = 20
    IGV.run <<-EOF, reference
new
genome #{reference}
snapshotDirectory #{files_dir}
load #{tumor}
load #{normal}
maxPanelHeight #{panel_height}
#{positions.collect do |position|
  chr, pos, alt = position.split(":")

  start = [pos.to_i - margin, 0].max
  eend = pos.to_i + margin + 1

  ["goto #{chr}:#{start}-#{eend}",
  "snapshot #{chr}:#{pos}:#{alt}.png"]

end * "\n"
}
exit
    EOF
    Dir.glob(File.join(files_dir, "*"))
  end

end
