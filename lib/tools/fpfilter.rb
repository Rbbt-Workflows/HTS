require 'rbbt-util'
require 'rbbt/workflow'
module FPFilter
  extend Resource
  
  Rbbt.claim Rbbt.software.opt.fpfilter, :install, Rbbt.share.install.software.fpfilter.find

  CMD.tool "bam-readcount", nil, "bash -c 'type bam-readcount'" do
    CMD.cmd('conda install bam-readcount -c bioconda')
  end

  def self.fixed_script(path=nil)
    fpfilter_cmd_original = Rbbt.software.opt.fpfilter.produce["fpfilter.pl"].find
    fpfilter_cmd_fixed = fpfilter_cmd_original.replace_extension('pl', 'fix.pl')

    if ! fpfilter_cmd_fixed.exists?
      code = fpfilter_cmd_original.read

      # New, safer implementation of read_counts_by_allele
      new_sub = <<-PERL
sub read_counts_by_allele {
    (my $line, my $allele) = @_;

    my @lineContents = split(/\\t/, $line);
    my $numContents = @lineContents;

    # Allele columns in bam-readcount begin at column index 4 (0-based):
    # chrom(0), pos(1), ref(2), depth(3), allele1(4)...
    for(my $colCounter = 4; $colCounter < $numContents; $colCounter++) {
        my $this_allele = $lineContents[$colCounter];
        next unless defined $this_allele && length $this_allele;
        # extract allele key only up to the first colon to avoid extra-colon issues
        my ($allele_key, $rest) = split(/:/, $this_allele, 2);
        if(uc($allele_key) eq uc($allele)) {
            my @alleleContents = split(/:/, $this_allele);
            # bam-readcount per-allele output should have at least 8 fields; otherwise we consider it incomplete
            return(\"\") if(@alleleContents < 8);
            # join numeric fields with tabs (fields 1..end)
            my $return_string = join(\"\\t\", @alleleContents[1..$#alleleContents]);
            return($return_string);
        }
    }

    return(\"\");
}
PERL

      # Locate the existing sub read_counts_by_allele and replace it using brace balancing
      match = code.match(/sub\\s+read_counts_by_allele\\s*\\{/)
      if match
        start_index = match.begin(0)
        # find the position of the opening brace
        open_brace_index = code.index('{', match.end(0)-1)
        if open_brace_index.nil?
          raise "Could not find opening brace for read_counts_by_allele"
        end

        # balance braces to find matching closing brace
        depth = 0
        i = open_brace_index
        last_index = nil
        while i < code.length
          ch = code.getbyte(i).chr
          if ch == '{'
            depth += 1
          elsif ch == '}'
            depth -= 1
            if depth == 0
              last_index = i
              break
            end
          end
          i += 1
        end

        if last_index.nil?
          raise "Could not find closing brace for read_counts_by_allele"
        end

        # Replace from start_index up to last_index (inclusive) with new_sub
        fixed_code = code[0...start_index] + new_sub + code[(last_index+1)..-1]
      else
        raise "read_counts_by_allele subroutine not found in original fpfilter.pl"
      end

      # Write out the fixed file
      fpfilter_cmd_fixed.write fixed_code
      FileUtils.chmod 0755, fpfilter_cmd_fixed.to_s
    end
    
    fpfilter_cmd_fixed
  end


  def self.filter(args)
    fpfilter_cmd = fixed_script
    CMD.get_tool 'bam-readcount'
    CMD.cmd_log("perl #{fpfilter_cmd}", args.to_hash.merge('add_option_dashes' => true))
  end

end

if __FILE__ == $0
  Log.severity = 0
end


