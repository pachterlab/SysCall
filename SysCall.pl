#!/usr/local/bin/perl
use strict;



if ($#ARGV != 3) { die "usage:  <place_file> <read_file(sam)> <outfile> <code_dir> \n";}

my @args = @ARGV;

main(@args);

sub main{
	my ($pl_file,$read_f,$outfile,$dir)=@_;

	print "checking input files...\n";
	Check($pl_file,$read_f);
	print "done checking.\n";

	## Output a file with a line for each read that overlaps position

 	my $p_reads = $outfile.".reads_parsed";
	my $command = "perl $dir"."parse_sam_reads.pl $read_f $pl_file $p_reads"; 
	RC($command);

	## make tables from the places file 

	my $list_2 = $outfile.".table_all";
	my $command = "perl $dir"."MakeTableFile.pl $pl_file $p_reads $list_2";
	RC($command);

	# this command gets the coverrage of the sites at hand, to choose which coverage to run the classifier on
        my $cov_file = GetCoverageFile($list_2);

	my $coverage_file = "$dir"."/Coverage_Files/$cov_file";	

	## choose features for each position # takes only sites w. at least one error.  
	
	my $list_3 = $outfile.".table_chosen";
	my $command = "perl $dir"."ChooseFeatures.pl $list_2 $list_3";
	RC($command);

	## Fit model and run on test set. 
	
	my $out_fs = $outfile.".sys_errors";
	open(T,">$out_fs") or die "can't open file $out_fs for writing\n";
	close(T);
	my $out_fh = $outfile.".heterozygous";
	open(T,">$out_fh") or die "can't open file $out_fh for writing\n";
	close(T);

	my $list_3_nn = $list_3.".nn";
	my $command = "cat $list_3 | grep -v NaN >$list_3_nn";
	my $ret_nn = RC2($command); 

	print "ret_nn: $ret_nn \n";

	
	if (0==$ret_nn){	
		my $out_fs_nn = $outfile.".sys_errors.nn";
		my $out_fh_nn = $outfile.".heterozygous.nn";
		my $command = "cat $dir"."classify_nn.R | R --slave --args $list_3_nn $coverage_file $out_fs_nn $out_fh_nn";
		RC($command); 		
		my $command = "cat $out_fs_nn >>$out_fs";
		RC($command); 		
		my $command = "cat $out_fh_nn >>$out_fh";
		RC($command); 		
	}


	my $list_3_n = $list_3.".n";
	my $command = "head -1 $list_3 >$list_3_n";
	RC($command);
	my $command = "cat $list_3 | grep NaN >>$list_3_n";
	my $ret_n = RC2($command); 
	print "ret_n: $ret_n \n";

	if (0==$ret_n){
		my $out_fs_n = $outfile.".sys_errors.n";
		my $out_fh_n = $outfile.".heterozygous.n";
		my $command = "cat $dir"."classify_n.R | R --slave --args $list_3_n $coverage_file $out_fs_n $out_fh_n";
		RC($command); 		
		my $command = "cat $out_fs_n >>$out_fs";
		RC($command); 		
		my $command = "cat $out_fh_n >>$out_fh";
		RC($command); 		
	}


}

sub GetCoverageFile{
        my ($file) = @_;

	my $thresh1 = 59;
	my $thresh2 = 46;
	my $thresh3 = 33;
	my $thresh4 = 20;

	my $fc1 = "coverage_66x_train";
	my $fc2 = "coverage_53x_train";
	my $fc3 = "coverage_40x_train";
	my $fc4 = "coverage_27x_train";
	my $fc5 = "coverage_13x_train";

        open(T,"<$file") or die "can't open file $file\n";
        my $line;

        my $num_rows=0;
        my $num_reads=0;

        $line=<T>;
        while ($line=<T>){

                chomp($line);
                my @arr = split(/ /, $line);

                my $n1 = $arr[7];
                my $n2 = $arr[9];

                if (($n1 =~ /^d+$/) || ($n2 =~ /^d+$/)){
                        die "ERROR not numbers : $line\n";
                }

                $num_reads+=$n1;
                $num_reads += $n2;

                $num_rows++;
        }

        close(T);

        my $ave = $num_reads/$num_rows;
	print "## The average coverage of your candidate heterozygous sites is: $ave \n";


        my $file_cov;
        if ($ave>$thresh1){
                $file_cov = $fc1;
        } elsif ( $ave > $thresh2 ) {
                $file_cov = $fc2;
        } elsif ($ave >$thresh3) {
                $file_cov = $fc3;
        } elsif ($ave > $thresh4) {
                $file_cov = $fc4;
        } else {
                $file_cov = $fc5;
     }

        return $file_cov;

}

sub Check{
	my ($pl_file,$read_f) = @_;


	# tests that pl_file is of the right format
	open (PL,"<$pl_file") or die "can't open file $pl_file \n";
	my $line;
	while($line=<PL>){
		chomp($line);
		my @arr = split(/ /, $line);
		if (!(7==@arr)){
			die "ERROR: file $pl_file is not in right format. There should be 7 elements in line space-seperated. \n";
		}
		if (!($arr[0] =~ /chr/)){
			die "ERROR: file $pl_file is not in right format. first place should hold chr \n";
		}
		if ($arr[1] = /^(\d+)$/){
			die "ERROR: file $pl_file is not in right format. expecting an integer in second column \n";
		}
		for(my $i=2; $i<7; $i++){
			my $tmp = $arr[$i];
			if (!(1==CheckBase($tmp))){
				die "ERROR: file $pl_file is not in right format. place $i doe not have a base\n";
			}
		}
	}
	close(PL);

	# tests the read file
	print "checking read file is right format...\n";	
	open (RF,"<$read_f") or die "can't open file $pl_file \n";
	my $line;
	for (my $i=0; $i<500;$i++){
		if ($line=<RF>){
			chomp($line);
			

			my @arr = split(/\t/, $line);
			if (!(11 < @arr)){
				die "ERROR: file $read_f not of right format. There should be at least 11 elements in tab-seperated line \nline: $line\n";
			}
			if ((!(($arr[2]=~/chr/)||($arr[2]=~/c(\d+|X|Y)\.fa/)))||(!($arr[3]=~/^(\d+)$/))){
				die "ERROR: file $read_f not of right format. Either place 2 is not chr or place 3 is not coordinate for line $line \n";
			}
		}

	}
	close(RF);

}


sub CheckBase{
	my ($let) = @_;

	if (!(($let eq 'A')||($let eq 'C')||($let eq 'G')||($let eq 'T'))){
		return -1;
	} else {
		return 1;
	}

}



sub FinalOut{
	my ($place_file,$out_tmp,$out_final)=@_;

	my %hash;
	open(OF,"<$out_tmp") or die "can't open file $out_tmp \n";
	my $line;

	
	while($line=<OF>){
		chomp($line);
		
		my ($pl,$prob) = split(/ /, $line);

		$hash{$pl}=$prob;
	}
	close(OF);

	
	open(IN,"<$place_file") or die "can't open file $place_file \n";
	open(OUT,">$out_final") or die "can't open file $out_final \n";
	
	$line=<IN>;
	my $num_line=1;
	while($line=<IN>){
		chomp($line);
		
		if (exists $hash{$num_line}){
			my @arr = split(/ /, $line);
			my $pos = $arr[0];
			my $pr = $hash{$num_line};			

			print OUT "$pos $pr\n";
			#print OUT "$pos\n";

		}
		$num_line++;
	}
	close(IN);
	close(OUT);


}


sub RC{

	my ($command) = @_;
	print "running command: $command\n";
	0==system("$command") or die "command: $command failed \n";


}

sub RC2{

	my ($command) = @_;
	print "running command: $command\n";
	my $r = system("$command") ;
	return $r;

}













