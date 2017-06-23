#!/usr/local/bin/perl

#-------------------------------------------------------------------------------------------------------------#
# Takes in a file of positions as: chr pos base-2 base-1 base base+1 base+2    - for the sence strand 
# 
# Incorporates read data of locations to print to outfile 
#-------------------------------------------------------------------------------------------------------------#

use strict;

if ($#ARGV != 2) { die "usage: <positions_file> <parsed_reads_file> <outfile> \n";}

my @args = @ARGV;

main(@args);


sub main{
	my ($pos_file,$reads_file,$outfile)=@_;

	my %hash_pos;
	ReadInPossitions($pos_file,\%hash_pos);
	
	open(IN,"<$reads_file") or die "can't open file $reads_file \n";
	open(OUT,">$outfile") or die "can't open file $outfile \n";
	
	print OUT "chr pos base-1 base-1 base base1 base2 num_+ num_e_+ num_- num_e_- t_stat_forward t_stat_backward\n";

	my %hash_taken;
	my @arr_order;

	my $line;
	while($line=<IN>){
		chomp($line);
		my @arr = split(/ /, $line);

		# place 5: num reads +  # place 6: num errors +   # place 7: num reads -  # place 8 : num, errors -  # place9: array bscores + #place10: array scores + place11: array bscores - #places12: array scores -  
		my ($chr,$pos,$b,$dir,$qi_bef,$qi,$rest) = @arr;
		my $key = "$chr".":$pos";
		if (!(exists $hash_pos{$key})){
			die "ERROR: position taken for read is not in hash: $line \n";
		}
	
		if (!(exists $hash_taken{$key})){
			push(@arr_order,$key);
			$hash_taken{$key}=1;
		}
		my $t_arr = $hash_pos{$key};
		my $true_base = $$t_arr[2];

		if ($dir eq "+"){
			
			$$t_arr[5]++;
			if (!($true_base eq $b)){
				$$t_arr[6]++;
			}	
			if (!($qi_bef eq "NaN")){
				if (!($b eq "N")){
					if ($qi eq "NaN"){
						print "## NOTICE: quality score is NaN at line: $line \n";
					} else {
						if ($qi==0){
							print "note: quality score is 0 for line: $line \n";
						}	
						my $arr_bq_scores = $$t_arr[9];
						push(@$arr_bq_scores,$qi_bef);
						my $arr_q_scores = $$t_arr[10];
						push(@$arr_q_scores,$qi);
					}
				}
			}	 

		} elsif ($dir eq "-"){
			$$t_arr[7]++;
			if (!($true_base eq $b)){
				$$t_arr[8]++;
			}	
			if (!($qi_bef eq "NaN")){
				if (!($b eq "N")){
					if ($qi eq "NaN"){
						print "## NOTICE: quality score is NaN at line: $line \n";
					} else {
						if ($qi==0){
							print "note: quality score is 0 for line: $line \n";
						}	

						my $arr_bq_scores = $$t_arr[11];
						push(@$arr_bq_scores,$qi_bef);
						my $arr_q_scores = $$t_arr[12];
						push(@$arr_q_scores,$qi);
					}
				}
			}	 

		} else {
			die "ERROR \n";
		}

	}	

	close(IN);

	foreach my $k (@arr_order){
		my $arr_vals = $hash_pos{$k};
	
		# getting t-stat 
		my $arr_bq_scores_f = $$arr_vals[9];
		my $arr_q_scores_f = $$arr_vals[10];
		my $arr_bq_scores_b = $$arr_vals[11];
		my $arr_q_scores_b = $$arr_vals[12];
		
		my ($t_p_stat_f,$t_p_stat_b);
		if (!(0==@$arr_bq_scores_f)) {
			$t_p_stat_f = GetPairedTStat($arr_bq_scores_f,$arr_q_scores_f);
		} else {
			$t_p_stat_f="NaN";
		}
		
		if (!(0==@$arr_bq_scores_b)) {
			$t_p_stat_b = GetPairedTStat($arr_bq_scores_b,$arr_q_scores_b); 
		} else {
			#print "Having NaN as t_p_stat_b for $k \n"; 
			$t_p_stat_b="NaN";
		}
		
		
		pop(@$arr_vals);
		pop(@$arr_vals);
		pop(@$arr_vals);
		pop(@$arr_vals);

		my ($chr,$pos) = split(/:/, $k);
		my $lp = join(' ',@$arr_vals);
		print OUT "$chr $pos $lp $t_p_stat_f $t_p_stat_b\n";


	}	
	close(OUT);

}


sub GetPairedTStat{
	my ($arr_bq,$arr_q) = @_;
	
	

	my @arr_diff;
	my $len = @$arr_bq;
	if (!($len==@$arr_q)){
		die "error: arrays not of same length \n";
	} 	

	my $bin_arr_eq=1;
	my $sum_diff;
	my $sum_sq_diff;
	for(my $i=0;$i<$len;$i++){
		#my $d = $$arr_bq[$i]-$$arr_q[$i];
		my $d = $$arr_q[$i]-$$arr_bq[$i];
		if (!(0==$d)){
			$bin_arr_eq=0;
		}
		push(@arr_diff,$d);
		$sum_diff += $d;
		$sum_sq_diff += ($d*$d);
	}

	if (1==$bin_arr_eq){
		print "NOTICE: arrays are equal returning t-stat 0. \n";
		return 0;
	}

	my $sum_diff_sq = ($sum_diff*$sum_diff);

	my $enumerator = $sum_diff;
	
	if (1==$len){
		return "NaN";
	}
	my $denom  = sqrt((($len*$sum_sq_diff) - $sum_diff_sq   )/($len-1));	

	if (0==$denom){
		print "############################ NOTICE: denome is 0: $denom. \n";
		print "arr_bq:\n";
		foreach my $e (@$arr_bq){
			print "$e ";
		}
		print "arr_q: \n";
		foreach my $e (@$arr_q){
			print "$e ";
		}
		return "NaN";
	}


	my $t_stat = $enumerator/$denom;

	return $t_stat;
}




sub ReadInPossitions{
	my ($inp,$hash_p) = @_;

	open (INP,"<$inp") or die "can't open file $inp \n";

	my $line;
	while($line=<INP>){
		chomp($line);
		
		my @arr = split(/ /, $line);
		my $place = $arr[1];
		my $key = "$arr[0]".":$place";
		shift(@arr);
		shift(@arr);
		if (!(5==@arr)){
			die "ERROR: length arr not 5: $line \n";
		}
		$arr[5]=0;
		$arr[6]=0;
		$arr[7]=0;
		$arr[8]=0;
		my @new_arr1;
		$arr[9] = \@new_arr1;
		my @new_arr2;
		$arr[10] = \@new_arr2;
		my @new_arr11;
		$arr[11] = \@new_arr11;
		my @new_arr22;
		$arr[12] = \@new_arr22;

		$hash_p->{$key}=\@arr;
	}
	close(INP);

}




















