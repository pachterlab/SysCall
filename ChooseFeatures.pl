#!/usr/local/bin/perl

#-------------------------------------------------------------------------------------------------------------#
# read in lines. 
# Take direction in which the proportion of errors is larger. Get base-2 base-1 base0 seq  (diff_cur_dir/all_cur_dir)-(diff_other_dir/all_other_dir) #diff_cur_dir/#all_cur_dir t-test Y/N
# Takes only sites with at least one error
#-------------------------------------------------------------------------------------------------------------#



use strict;

if ($#ARGV != 1) { die "usage: <infile>  <outfile> \n";}

my @args = @ARGV;

main(@args);

sub main{
	my ($inf,$outfile)=@_;

	open (IN,"<$inf") or die "can't open file $inf \n";
	open (OUT,">$outfile") or die "can't open file $outfile \n";

	my $line;
	
	print OUT "id base-2 base-1 base diff_err_diff_dir diff_error_dir t_test\n"; 

	while($line=<IN>){
		chomp($line);
		
		my @arr = split(/ /, $line);
		if (!(13==@arr)){
			die "ERROR: not right number of elements \n";
		}
		
		my $err_f = $arr[8];
		my $all_f = $arr[7];
		my $err_b = $arr[10];
		my $all_b = $arr[9];


		if ($err_f>$all_f){
			die "ERROR: err_f > all_f : $err_f  $all_f \n";
		}
		if ($err_b>$all_b){
			die "ERROR: err_b > all_b : $err_b  $all_b \n";
		}

		
		my $t_test_f = $arr[11];
		my $t_test_b = $arr[12];
	
		my $take_f=-1;
		

		### 
		if (($err_f>0)||($err_b>0)){
			if (0==$all_f){
				$take_f=0;
			} elsif (0==$all_b){
				$take_f=1;
			} elsif (($err_f/$all_f) > ($err_b/$all_b)){
				$take_f=1;
			} elsif (($err_f/$all_f) <= ($err_b/$all_b)) {
				$take_f=0;
			}



			if (1==$take_f){
				# taking forward direction

				my $diff_dirs_errors;
				if (0==$all_b){
					$diff_dirs_errors = ($err_f/$all_f);
				} else {
					$diff_dirs_errors = ($err_f/$all_f)-($err_b/$all_b);
				}

				my $ratio_error_cdir = $err_f/$all_f;

				my $n_l = "$arr[0]:$arr[1] $arr[2] $arr[3] $arr[4] $diff_dirs_errors $ratio_error_cdir $arr[11]"; 
				print OUT "$n_l\n";


			} elsif(0==$take_f) {

				# taking reverse direction
				my $diff_dirs_errors;
				if (0==$all_f){
					$diff_dirs_errors = ($err_b/$all_b);

				} else {
					$diff_dirs_errors = ($err_b/$all_b)-($err_f/$all_f);
				}
			
				my $ratio_error_cdir = $err_b/$all_b;
				my $rev1 = Reverse($arr[6]);
				my $rev2 = Reverse($arr[5]);
				my $rev3 = Reverse($arr[4]);


				my $n_l = "$arr[0]:$arr[1] $rev1 $rev2 $rev3 $diff_dirs_errors $ratio_error_cdir $arr[12]"; 
				print OUT "$n_l\n";

			} else {
				die "ERROR: take_f not right: $take_f \n";

			}



		}
	}
}


sub Reverse{
	my ($let) = @_;
		
	if ($let eq "A"){
		return "T";
	}	
	if ($let eq "C"){
		return "G";
	}	
	if ($let eq "G"){
		return "C";
	}	
	if ($let eq "T"){
		return "A";
	}	

	die "ERROR: didn't fine letter  $let \n";

}


