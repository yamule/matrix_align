use strict;
use warnings;

# CATH の FASTA ヘッダにクラスの情報をつける
# ファイル名などは全部ハードコードしている
# レポジトリのルート直下で RUN する

my $cath_domain_list="nogit/large_data/id20/fas/cath-domain-list.txt";
my $infile = "nogit/large_data/id20/fas/cath-dataset-nonredundant-S20.fa";
my $outfile = "nogit/large_data/id20/fas/cath-dataset-nonredundant-S20_assigned.fa";

open(IN,$infile);
my %processed;
while(my $ss = <IN>){
    if($ss =~ />([^\s]+)/){
        my $pc = $1;
        my @ptt = split(/[\|\/]/,$ss);
        my $seqname = $ptt[2];

        if(defined $processed{$seqname}){
            print STDERR $seqname." already found.\n";
            die;
        }
        $processed{$seqname} = 100;
    }
}
close(IN);

my %name_to_class;
open(IN,$cath_domain_list);
# Column 1:  CATH domain name (seven characters)
# Column 2:  Class number
# Column 3:  Architecture number
# Column 4:  Topology number
# Column 5:  Homologous superfamily number
# Column 6:  S35 sequence cluster number
# Column 7:  S60 sequence cluster number
# Column 8:  S95 sequence cluster number
# Column 9:  S100 sequence cluster number
# Column 10: S100 sequence count number
# Column 11: Domain length
# Column 12: Structure resolution (Angstroms)
#            (999.000 for NMR structures and 1000.000 for obsolete PDB entries)
while(my $ss = <IN>){
    $ss = $ss;
    if($ss =~ /^#/){
        next;
    }
    my @ptt = split(/[\s]+/,$ss);

    if(!defined $processed{$ptt[0]}){
        next;
    }
    my $classs = $ptt[1].".".$ptt[2].".".$ptt[3].".".$ptt[4];
    if(defined $name_to_class{$ptt[0]}){
        if($name_to_class{$ptt[0]} ne $classs){            
            print STDERR $ss;
            die;
        }
    }
    $name_to_class{$ptt[0]} = $classs;

}
close(IN);

open(IN,$infile);
open(OUT,"> ".$outfile);
while(my $ss = <IN>){
    if($ss =~ />([^\s]+)/){
        my $pc = $1;
        my @ptt = split(/[\|\/]/,$ss);
        my $seqname = $ptt[2];
        if(!defined $name_to_class{$seqname}){
            print STDERR $seqname."\n";
            print STDERR $ss;
            die;
        }
        print OUT ">".$seqname."  ".$name_to_class{$seqname}."\n";
    }else{
        print OUT $ss;
    }
}
close(OUT);
close(IN);
