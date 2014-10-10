#!usr/bin/perl
$file1 = $ARGV[0];
open (PFAM, $file1) || die;
@pfam = <PFAM>;
close PFAM;

for ($i = 0; $i < @pfam; $i++) {
        chomp $pfam[$i];
        @linha = split (/\t/, $pfam[$i]);
        if ($i == 0) {
        $id = $linha[0];
        $score = $linha[1];
        }
        elsif ($linha[0] eq $id) {
                        if ($linha[1] >= $score){
				$score = $linha[1];
        }
	}
        else {
         @col = split(/\t/, $id);
        print "$col[0]*$col[1]\t$score\n";

        $id = $linha[0];
        $score = $linha[1];
	
        
	}
}
        print "$id\t$score\n";
exit;

