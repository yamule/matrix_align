use strict;
use warnings;

# perl scripts/calc_sp_score.pl <targetfile> <reference file>

sub load_seq{
	my $filename = $_[0];
	my @ret;
	my $len = -1;
	my $seqname = "dummy";
	open(FAS_R_IN2,$filename)or die $filename;
	while(my $tt = <FAS_R_IN2>){
		if($tt =~ /^[\s]*>/){
			if($tt =~ /^[\s]*>[\s]*([^\s]+)[\s]*/){
				push(@ret,$tt);
				$len = @ret-1;
			}else{
				push(@ret,$tt);
				$len = @ret-1;
			}
		}else{
			if($len > -1){
				$ret[$len] .= $tt;
			}
		}
	}
	close(FAS_R_IN2);
	return \@ret;
}


sub get_desc_fasta{
	my $f_dat = $_[0];
	if($f_dat =~ /^[\s]*>[\s]*[^\s]*[ \t]*([^\r^\n]*)[\r\n]*/){
		return $1;
	}
	return "";
}
sub get_name_fasta{
	my $f_dat = $_[0];
	if($f_dat =~ /^[\s]*>[\s]*([^\s]*)[\s]*/){
		return $1;
	}
	return "";
}

sub get_seq_fasta{
	my $f_dat = $_[0];
	my $f_ret = "";
	my @f_line = split(/[\r\n]+/,$f_dat);
	for(my $ii = 0;$ii < (my $len = @f_line);$ii++){
		if($f_line[$ii] !~ />/){
			$f_ret .= $f_line[$ii];
		}
	}
	return $f_ret;
}




sub get_aa_fasta{
	my $f_dat = $_[0];
	my $f_ret = "";
	my @f_line = split(/[\r\n]+/,$f_dat);
	for(my $ii = 0;$ii <= $#f_line;$ii++){
		if($f_line[$ii] !~ />/){
			$f_ret .= $f_line[$ii];
		}
	}
	$f_ret =~ s/[\s]//g;
	if($f_ret =~ /[^\-A-Za-z]/){
		die "Bad letter ".$&."\n";
	}
	return $f_ret;
}


my @seq1=@{load_seq($ARGV[0])};
my @seq2=@{load_seq($ARGV[1])};

my %aa_to_align;
foreach my $ss(@seq1){
	my $seq = uc get_aa_fasta($ss);
	my $aa = uc $seq;
	$aa =~ s/[^A-Z]//g;
	my @pseq = split(//,$seq);
	if(defined $aa_to_align{$aa}){
		die "Duplicate sequence";
	}else{
		my @tmp;
		$aa_to_align{$aa} = \@tmp;
		push(@{$aa_to_align{$aa}},\@pseq);
	}
}

my %seq_to_name;
my %seq_in_ref;
foreach my $ss(@seq2){
	my $seq = uc get_aa_fasta($ss);
	my $aa = uc $seq;
	$aa =~ s/[^A-Z]//g;
	$seq_to_name{$aa} = get_name_fasta($ss);
	my @pseq = split(//,$seq);
	if(!defined $aa_to_align{$aa}){
		die $aa." was not found in file1.\n";
	}else{
		push(@{$aa_to_align{$aa}},\@pseq);
		$seq_in_ref{$aa} = 100;
	}
}
my @ssorted = sort{$a cmp $b}keys %seq_in_ref;
my $pnum = $#ssorted+1;
my $snum = 2;# �y�A���C�Y����
for(my $ii = 0;$ii < $pnum;$ii++){
	if($snum != 2){
		die $ssorted[$ii]." has ".($#{$aa_to_align{$ssorted[$ii]}}+1)." seqs. Expected 2 .\n";
	}
}
# �������̔z���^���āA�}�b�`��Ԃł���ꍇ�A�������̃}�b�`���Ă��镶���̃C���f�N�X�[>�������̃}�b�`���Ă��镶���̃C���f�N�X
# �̃n�b�V����Ԃ�
sub posmap{
	my @s1 = @{$_[0]};
	my @s2 = @{$_[1]};
	
	if($#s1 != $#s2){
		die "The length should be the same.".join("",@s1)."\n".join("",@s2);
	}
	my %ret;
	my $pos1 = 0;
	my $pos2 = 0;
	for(my $ii = 0;$ii <= $#s1;$ii++){
		if($s1[$ii] ne "-" && $s2[$ii] ne "-"){
			$ret{$pos1} = $pos2;
		}
		if($s1[$ii] ne "-"){
			$pos1 ++;
		}
		if($s2[$ii] ne "-"){
			$pos2 ++;
		}
	}
	
	#print "$pos1 $pos2 ==\n";
	#foreach my $kk(sort{$a <=> $b}keys %ret){
	#	print $kk."\t".$ret{$kk}."\n";
	#}
	return \%ret;
}

my $norm_all = 0;
my $match_all = 0;
for(my $ii = 0;$ii < $pnum;$ii++){
	for(my $jj = $ii+1;$jj < $pnum;$jj++){
		my %targetmap = %{posmap(${$aa_to_align{$ssorted[$ii]}}[0],${$aa_to_align{$ssorted[$jj]}}[0])};
		my %refmap = %{posmap(${$aa_to_align{$ssorted[$ii]}}[1],${$aa_to_align{$ssorted[$jj]}}[1])};
		my $aname = $seq_to_name{$ssorted[$ii]};
		my $bname = $seq_to_name{$ssorted[$jj]};
		
		#my $norm = length($ssorted[$ii]);
		#if( length($ssorted[$jj]) < length($ssorted[$ii])){
		#	$norm = length($ssorted[$jj]); # min load ����̂߂�ǂ������B
		#}
		my $matcher = 0;
		my $norm = 0;
		foreach my $kk(keys %refmap){
			if(defined $targetmap{$kk}){
				if($refmap{$kk} == $targetmap{$kk}){
					$matcher ++;
				}
			}
			$norm ++;
		}
		print $aname."\t".$bname."\n";
		print $matcher/$norm;
		print "\n";

		$norm_all += $norm;
		$match_all += $matcher;
	}
}

print $match_all/$norm_all;
print "\n";
