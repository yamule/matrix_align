use strict;
use warnings;

# 
open(OUT,"> nogit/large_data/id20/fas/scope20_upper.fas");
open(IN,"nogit/large_data/id20/fas/astral-scopedom-seqres-sel-gs-bib-20-2.08.fa");
while(my $ss = <IN>){
	if($ss =~ /(>[\s]*[^\s]+[\s]+[^\s]+)/){ # BLAST に与えたときに  Bad char と言われたので その対応
		$ss = $1."\n";
	}else{
		$ss = uc $ss;
		$ss =~ s/[JUZBOX]/X/g;
	}
	print OUT $ss;
}
close(IN);
close(OUT);


open(OUT,"> nogit/large_data/id20/fas/cath20_upper.fas");
open(IN,"nogit/large_data/id20/fas/cath-dataset-nonredundant-S20_assigned.fa"); #assign_cath.pl を先に RUN しておく。
while(my $ss = <IN>){
	if($ss =~ />/){
	}else{
		$ss = uc $ss;
		$ss =~ s/[JUZBOX]/X/g;
	}
	print OUT $ss;
}
close(IN);
close(OUT);
