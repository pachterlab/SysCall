#!/usr/local/bin/perl
use strict;

#-------------------------------------------------------------------------------------------------------------#
# Input: fragments_sam_file  positions output
# Output:  for each fragment overlapping one of the possitions: chr pos base_read direction quality_int+1 quality_int
#-------------------------------------------------------------------------------------------------------------#

if ($#ARGV != 2) { die "usage: <infile> <possitions> <outfile>  \n";}

my @args = @ARGV;

main(@args);

sub main{
	my ($inf,$pos_file,$outf)=@_;

	my $n_not_taken=0;

	my %hash_pos;
	ReadInPossitions($pos_file,\%hash_pos);

	open(OUT,">$outf") or die "can't open file $outf \n";
	open(IN,"<$inf") or die "can't open file $inf \n";	

	my $line;
	my $c_l=0;	
	my $dir_d;
	my $chr_f;
	my ($bin_flag,$chr,$cord,$dir_1,$dir,$start_cord,$end_cord,$sequence,$len_seq,$q_s_seq);
	while($line=<IN>){
		chomp($line);
		
		if (0==($c_l % 100000)){
			print "at line $c_l \n";
		}
		$c_l++;

		my @arr = split(/\t/, $line);
		$bin_flag = $arr[1];
		$chr_f = $arr[2];
		$cord = $arr[3];
		$dir_1 = $arr[8];
		$sequence = $arr[9];
		$q_s_seq = $arr[10];

		$len_seq = length($sequence);

		if ($chr_f =~ /fa/){
			if (!($chr_f =~ /c(\d+|X|Y)\.fa/)){
				die "ERROR: SysCall can only work with chr fields such as (1) chr13 or (2) c13.fa .  \nPlease adjust your input .sam file. \n";
			} else {
				$chr = "chr"."$1";
			}
			
		} else {
			$chr=$chr_f;
		}


		if ($bin_flag & 16){
			$dir="-";
			$start_cord=$cord;
			$end_cord=$cord+ $len_seq -1;	
		} else {
			$dir="+";
			$start_cord=$cord;
			$end_cord=$cord+$len_seq-1;
		}

		
		#if ($dir_1>0){
		#	if (!($dir eq "+")){
				#print "\n\n### NOTICE: Contradicting directionality tags in read: $line \n";
		#	}
		#} elsif ($dir_1<0){
		#	if (!($dir eq "-")){
		#		print "\n\n### NOTICE: Contradicting directionality tags in read: $line \n";
		#	}
		#} 

		my $cords_arr = $hash_pos{$chr};	
		my ($take,$place_base_arr,$cord_error_arr);
		if (($end_cord < $$cords_arr[0]) || ($$cords_arr[@$cords_arr-1] < $start_cord)){ 
			$take=0;
		} else {
			($take,$place_base_arr,$cord_error_arr) = BinOverlap($start_cord,$end_cord,$cords_arr);
		}
		if (1==$take){
			my $i=0;
		
			#print "NOTICE: TAKE IS 1 \n";
			foreach my $pb (@$place_base_arr){
				my $cord_error = $$cord_error_arr[$i];
				my @arr_seq = split(//, $sequence);
				my $let = $arr_seq[$pb];
				
				my @arr_q_int = split(//, $q_s_seq); 
				my $q_asc  = $arr_q_int[$pb];
				my $q_int = ord($q_asc)-33;
				#print "q_asc: $q_asc   q_int: $q_int \n";
				my $q_asc_p_1;
				my $q_int_p_1;
				# taking the quality score of a basepair one downstream
				if ($dir eq "+"){
					if ((length($sequence)-1)==$pb){
						#print "NOTICE: error is present at the first place of read, printing NaN \n. line: $line \n"; 
						#print "NOTICE: error is present at first place in read, not taking place before.\n line: $line \n";
						$q_int_p_1="NaN";	
					} else {
						$q_asc_p_1=$arr_q_int[$pb+1];
						$q_int_p_1 = ord($q_asc_p_1)-33; 
					}
				} elsif ($dir eq "-") {
					if (0==$pb){
						#print "NOTICE: error is present at the first place of read, printing NaN \n. line: $line \n"; 
						#print "NOTICE: error is  at first place in read, not taking place before.\n line: $line \n";
						$q_int_p_1="NaN";	
					} else {
						$q_asc_p_1=$arr_q_int[$pb-1];
						$q_int_p_1 = ord($q_asc_p_1)-33; 
					}

				} else {
					die "unrecognized direction: $dir \n";
				} 
				print OUT "$chr $cord_error $let $dir $q_int_p_1 $q_int $arr[0]\n"; 
				$i++;
			}
		}



	}

	close(IN);
	close(OUT);


}

sub BinOverlap{

	my ($start,$end,$cords_arr)=@_;

	my @place_base_arr;
	my @cord_error_arr;
	my $take=0;	


	foreach my $cur_cord (@$cords_arr){
		if (($start<=$cur_cord) && ($cur_cord <=$end)){
			$take=1;
			push(@cord_error_arr,$cur_cord);
			my $place = $cur_cord-$start;
			push(@place_base_arr,$place);
		}

	}	

	return ($take,\@place_base_arr,\@cord_error_arr);

}








sub ReadInPossitions{
	my ($inp,$hash_p) = @_;

	open (INP,"<$inp") or die "can't open file $inp \n";

	my $line;
	my ($cur_chr,$cur_pos);
	while($line=<INP>){
		chomp($line);
		my @arr = split(/ /, $line);
		if (@arr!=7){
			die "ERROR: arr not 7! line: $line \n";
		}	
		my $cur_chr=$arr[0];
		my $cur_pos = $arr[1];		


		if (!(exists $hash_p->{$cur_chr})){
			my @new_arr;
			$hash_p->{$cur_chr}=\@new_arr;
		}			
		my $cur_arr = $hash_p->{$cur_chr};
		push(@$cur_arr,$cur_pos);	
	}
	
	close(INP);

	foreach my $chr (keys %$hash_p){
		my $t_arr = $hash_p->{$chr};
		my @sorted = sort {$a<=> $b} (@$t_arr);
		$hash_p->{$chr} = \@sorted;

	}

	
	#foreach my $chr (keys %$hash_p){
		#my $t_arr = $hash_p->{$chr};
		#my $l = join(' ',@$t_arr);
		#print "chr $chr : $l\n";
	#}

}





















