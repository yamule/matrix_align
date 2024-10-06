use strict;
use warnings;

# makeblastdb -dbtype prot -parse_seqids -in nogit/large_data/id20/fas/scope20_upper.fas 
# blastp -query ＜cath の FASTA＞ -db ＜SCOPE の FASTA で作成した DB＞ -out nogit/large_data/id20/fas/blast/cath20_to_scope20.blastres.dat -evalue 0.001 -outfmt 6 -num_threads 12  
# で作成した BLAST の結果
my $blastres ="nogit/large_data/id20/fas/blast/cath20_to_scope20.blastres.dat";
my $infile = "nogit/large_data/id20/fas/cath20_upper.fas";
my $outfile = "nogit/large_data/id20/fas/c20rs20s.fas";

open(IN,$blastres);
my %has_hit;
while(my $ss = <IN>){
    my @ptt = split(/[\s]+/,$ss);
    my $idd = $ptt[2];
    my $eva = $ptt[$#ptt-1];
    if($eva <= 0.001){
        if($idd >= 40.0){
            $has_hit{$ptt[0]} = 100;
        }
    }
}
close(IN);

open(OUT,"> $outfile");
open(IN,"$infile");
my %processed;
my $flag = 0;
while(my $ss = <IN>){
    if($ss =~ />/){
        $flag = 1;
        if($ss =~ />([^\s]+)/){
            my $seqname = $1;
            if(defined $has_hit{$seqname}){
                $flag = 0;
                next;
            }
        }else{
            die;
        }
    }else{
        $ss = uc $ss;
    }
    
    if($flag == 1){
        print OUT $ss;
    }
}
close(IN);
close(OUT);

