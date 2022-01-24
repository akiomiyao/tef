#!/usr/bin/perl
#
# tef.pl - transposable element finder
#          A program for transposable elements detection from NGS reads.
#
# Copyright (C) 2020 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/tef
#
# Author: Akio Miyao <miyao@affrc.go.jp>
#


if ($ARGV[0] eq ""){
    print "
 tef.pl - transposable element finder
      A program for detection of active transposable elements from NGS reads.

 e.g. perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0
      perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5,th=0.7
      perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,option=clear,max_process=8

 tef.pl requires a pair of NGS reads from different conditions.
 For example, ttm2 and ttm5 are from regenerated rice individuals from callus.
 Save fastq files into ttm2/read and ttm5/read directory.

 Options shold be specified with 'name=value' separated by a comma WITHOUT space.
 Default method is the junction method.
 If tsd_size is specified, tef.pl will run with TSD method.

 If ref is specified, gziped fasta file of reference genome should be saved
 into the directory specified with ref. 
 
 At the first time of run, config file for reference will be made in the
 ref directory. The config has a list of sequence name. If you do not want to
 include unassaigned contig, change the name to NOP and save the contig file,
 and run tef.pl again.

 If ref is not specified, a or b specific transpositions without insertion 
 position on the genome will be detected.
 This function works only by TSD method.

 If tsd_size is not specified, detection will be progressed by junction method.

 If th (threthold) on the TSD method is not specified, default value 0.2 will be used.
 If a lot of noize are detected, try again with higher value e.g. th=0.7  

 tmp directory in target is not deleted and can be reused.
 If option=clear is specified, data in tmp will be cleared at 
 the begining of analysis. 

 For Linux, max process is number of CPU core. For other OS, defalt process number is 4.
 If you add max_process option, e.g. max_process=8, tef use 8 cores. 

 If data in split and/or count.tsd_size are truncated, remove the directory in
 the target and then run again.

 Author: Akio Miyao

";
    exit;
}

@nuc = ('A', 'C', 'G', 'T');
$s = {};

$start_time = time;

$wd = `pwd`;
chomp($wd);

foreach (split(',', $ARGV[0])){
    my ($name, $val) = split('=', $_);
    $$name = $val;
}

open(IN, "/proc/cpuinfo");
while(<IN>){
    $processor ++  if /processor/;
}
close(IN);

$processor = 32 if $processor > 32;
$processor = $max_process if $max_process ne "";
$processor = 4 if $processor eq "";
$th = 0.2 if $th eq "";
$tsd_size = 20 if $tsd_size eq "";

if ($sub eq ""){
    system("rm -rf $wd/$a/tmp") if -e "$wd/$a/tmp" and $option =~ /clear/;
    system("rm -rf $wd/$b/tmp") if -e "$wd/$b/tmp" and $option =~ /clear/;
    system("rm -rf $wd/$a/child") if -e "$wd/$a/child";
    system("mkdir $wd/$a/tmp") if ! -e "$wd/$a/tmp";
    system("mkdir $wd/$b/tmp") if ! -e "$wd/$b/tmp";
    system("mkdir $wd/$a/split") if ! -e "$wd/$a/split";
    system("mkdir $wd/$b/split") if ! -e "$wd/$b/split";
    system("mkdir $wd/$a/child");
    &log("job start");
    &log("Argument : $ARGV[0]");
    if ($tsd_size == 20){
	$method = "Junction";
    }else{
	$method = "TSD";
    }
    &log("Method : $method");
    &commonMethod;
    if ($tsd_size == 20){
	&junctionMethod;
    }else{
	&tsdMethod;
    }
    &join;
    system("rm -rf $wd/$a/tmp") if $option =~ /clear/;
    system("rm -rf $wd/$b/tmp") if $option =~ /clear/;
    $end_time = time;
    $elapsed_time = $end_time - $start_time;
    $hour = int($elapsed_time / 3600);
    $min = $elapsed_time % 3600;
    $sec = $min % 60;
    $min = int($min / 60);
    if ($hour >= 24){
	$day = int($hour / 24);
	$hour = $hour % 24;
    }
    
    if ($day > 1){
	$etime .= "$day days ";
    }elsif($day == 1){
	$etime .= "$day day ";
    }
    if ($hour > 1){
	$etime .= "$hour hours ";
    }elsif($hour == 1){
	$etime .= "$hour hour ";
    }
    if ($min > 1){
	$etime .= "$min minutes ";
    }elsif($min == 1){
	$etime .= "$min minute ";
    }
    if ($sec > 1){
	$etime .= "$sec seconds ";
    }elsif($sec == 1){
	$etime .= "$sec second ";
    }
    &log("Job completed: a=$a, b=$b, ref=$ref, method=$method
$etime ($elapsed_time seconds) elapsed.");
}else{
    my $file = "$wd/$a/child/child.$$";
    system("touch $file");
    &$sub;
    system("rm $file");
    exit;
}

sub commonMethod{
    if ($ref ne ""){
	&mkref;
    }
    if (! -e "$wd/$a/split/split.1"){
	&log("split to subfiles : $a");
	$cmd = "perl $0 target=$a,sub=split,a=$a &";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : commonMethod : $cmd") if $rc;
    }
    &monitorWait;
    if (! -e "$wd/$b/split/split.1"){
	&log("split to subfiles : $b");
	$cmd = "perl $0 target=$b,sub=split,a=$a &";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : commonMethod : $cmd") if $rc;
    }
    &join;
}

sub junctionMethod{
=pod
    system("mkdir $wd/$a/count.20") if ! -e "$wd/$a/count.20";
    system("mkdir $wd/$b/count.20") if ! -e "$wd/$b/count.20";
    if (&fragmentLength($a) != 40){
	&log("count subfiles : $a");
	&count($a);
    }
    if (&fragmentLength($b) != 40){
	&log("count subfiles : $b");
	&count($b);
    }
    &join;
    if (&fragmentLength($a) != 40){
	&log("merge subfiles : $a");
	&merge($a);
    }
    if (&fragmentLength($b) != 40){
	&log("merge subfiles : $b");
	&merge($b);
    }
    &join;
    &log("specific");
    &junctionSpecific;
    &log("firstmap : $a");
    &junctionFirstMap($a);
    &log("firstmap : $b");
    &junctionFirstMap($b);
    &join;
    &log("secondmap : $a");
    &junctionSecondMap($a);
    &log("secondmap : $b");
    &junctionSecondMap($b);
    &join;
    &log("sort : $a");
    &junctionSort($a);
    &log("sort : $b");
    &junctionSort($b);
    &join;
    &log("select candidate : $a");
    &junctionSelectCandidate($a);
    &log("select candidate : $b");
    &junctionSelectCandidate($b);
    &join;
    &junctionMapSelection($a);
    &junctionMapSelection($b);
    &join;
    &junctionTsdSelection($a);
    &junctionTsdSelection($b);
    &join;
    &summary($a);
    &summary($b);
    &genotype($a);
    &genotype($b);
=cut
    &verify2;
}

sub tsdMethod{
    if (&fragmentLength($a) != 20 + $tsd_size){
	&log("count subfiles : $a");
	&count($a);
    }
    if (&fragmentLength($b) != 20 + $tsd_size){
	&log("count subfiles : $b");
	&count($b);
    }
    &join;
    if (&fragmentLength($a) != 20 + $tsd_size){
	&log("merge subfiles : $a");
	&merge($a);
    }
    if (&fragmentLength($b) != 20 + $tsd_size){
	&log("merge subfiles : $b");
	&merge($b);
    }
    &join;
    &log("comm");
    &comm;
    &join;
    &log("tsdJoin");
    &tsdJoin;
    &join;
    &log("tsd");
    &tsd;
    foreach $taga (@nuc){
	foreach $tagb (@nuc){
	    foreach $tagc (@nuc){
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		&log("pair $tag");
		$cmd = "perl $0 a=$a,b=$b,sub=pair,tag=$tag,tsd_size=$tsd_size &";
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : $cmd") if $rc;
	    }
	}
    }
    &join;
    &log("sortPair");
    &sortPair;
    foreach $taga (@nuc){
	foreach $tagb (@nuc){
	    foreach $tagc (@nuc){
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		&log("selectPair $tag");
		$cmd = "perl $0 a=$a,b=$b,sub=selectPair,tag=$tag,tsd_size=$tsd_size,th=$th &";
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : $cmd") if $rc;
	    }
	}
    }
    &join;
    #    system("cat $wd/$a/tmp/pair.$a.$b.$tsd_size.* > $wd/$a/tsd_method.pair.$a.$b.$tsd_size && rm $wd/$a/tmp/pair.$a.$b.$tsd_size.*");
    $cmd = "cat $wd/$a/tmp/pair.$a.$b.$tsd_size.* > $wd/$a/tsd_method.pair.$a.$b.$tsd_size";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $cmd") if $rc;
    if($ref eq ""){
	&verify;
    }else{
	&mkQuery;
	&map;
    }
    &verify2;
}

sub genotype{
    my $target = shift;
    &monitorWait;
    chomp;
    &log("output genotype data : $target");
    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=genotypeFunc,target=$target &";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : genotype : $cmd") if $rc;
}

sub genotypeFunc{
    my (@row, $seq, $pos, %count, %pos);
    opendir(DIR, $target);
    foreach (sort readdir(DIR)){
	if (/junction_method.$target.[ACGT]/){
	    next if -z "$target/$_";
	    open(IN, "$target/$_");
	    while(<IN>){
		chomp;
		@row = split;
		open(CHR, "$ref/chr$row[1]");
		binmode(CHR);
		if ($row[4] eq "f"){
		    seek(CHR, $row[2] - 20, 0);
		}else{
		    seek(CHR, $row[2] - 21, 0);
		}
		read(CHR, $seq, 40);
		close(CHR);
		$count{$seq} = 0;
		$pos{"$row[0] $row[1] $row[2] $row[3] wt"} = $seq;
		$seq = $row[8] . $row[5];
		$count{$seq} = 0;
		$pos{"$row[0] $row[1] $row[2] $row[3] head"} = $seq;
		$seq = $row[6] . $row[9];
		$count{$seq} = 0;
		$pos{"$row[0] $row[1] $row[2] $row[3] tail"} = $seq;
	    }
	    close(IN);
	}
    }
    closedir(DIR);
    open(IN, "zcat $target/count.20/*.gz |");
    while(<IN>){
	chomp;
	@row = split;
	if (defined $count{$row[0]}){
	    if ($row[1] eq ""){
		$count{$row[0]} ++;
	    }else{
		$count{$row[0]} += $row[1];
	    }
	}
    }
    close(IN);
    opendir(DIR, $target);
    foreach (sort readdir(DIR)){
	if (/junction_method.$target.[ACGT]/){
	    next if -z "$target/$_";
	    @row = split('\.', $_);
	    print "\njunction_method.genotype.$target.$row[2].$row[3]\n";
	    open(OUT, "> $target/junction_method.genotype.$target.$row[2].$row[3]");
	    open(IN, "$target/$_");
	    while(<IN>){
		my (@row, $hpos, $tpos, $wpos, $hcount, $tcount, $wcount, $genotype);
		chomp;
		@row = split;
		$hpos = "$row[0] $row[1] $row[2] $row[3] head";
		$tpos = "$row[0] $row[1] $row[2] $row[3] tail";
		$wpos = "$row[0] $row[1] $row[2] $row[3] wt";
		$hcount = $count{$pos{$hpos}};
		$tcount = $count{$pos{$tpos}};
		$wcount = $count{$pos{$wpos}};
		if ($wcount){
		    $genotype = "H";
		}else{
		    $genotype = "M";
		}
$output = "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$hcount\t$tcount\t$wcount\t$genotype\n";
		print $output;
		print OUT $output;
	    }
	    close(IN);
	}
    }
    closedir(DIR);	
}

sub summary{
    my $target = shift;
    my ($head, $tail, $chead, $ctail, %file, %ht, %uniq, %TSD, %total, $count, $pair, %done);
    open($target, "> $target/junction_method.pair.$target");
    print "
Target: $target
Head                  Tail                  Kinds/Total\tTSDs\n";
    print $target "Target: $target
Head                  Tail                  Kinds/Total\tTSDs\n";
    opendir(DIR, $target);
    foreach (sort readdir(DIR)){
	if (/junction_method.$target.[ACGT]/){
	    my $file = $_;
	    my %tsd = ();
	    my $total = 0;
	    my $uniq = 0;
	    my $tsd = "";
	    open(IN, "$target/$file");
	    while(<IN>){
		@row = split;
		$tsd{$row[7]} ++;
		$total ++
	    }
	    close(IN);
	    foreach (sort keys %tsd){
		$uniq ++;
		$tsd .= "$_ ";
	    }
	    if ($total > 0){
		($head, $tail) = (split('\.', $file))[2,3];
		$file{"$head $tail"} = 1;
		$uniq{"$head $tail"} = $uniq;
		$TSD{"$head $tail"} = $tsd;
		$total{"$head $tail"} = $total;
	    }
	}
    }
    closedir(DIR);
    foreach (sort keys %file){
	($head, $tail) = split;
	$pair = "$head $tail";
	print "$head  $tail\t$uniq{$pair}/$total{$pair}\t$TSD{$_}\n";
	print $target "$head  $tail\t$uniq{$pair}/$total{$pair}\t$TSD{$_}\n";
    }
    
    $count = 0;
    foreach (sort keys %file){
	$count++;
    }
    
    $count = int($count / 2 + 0.5);
    
    print "Number of TE candidates: $count\n";
    print $target "Number of TE candidates: $count\n";

}

sub fragmentLength{
    my $target = shift;
    my $fragment_length = 0;
    if (-e "$target/count.$tsd_size/TTT.gz"){
	open(IN, "zcat $target/count.$tsd_size/TTT.gz |");
	while(<IN>){
	    chomp;
	    @row = split;
	    $fragment_length =length($row[0]);
	    last;
	}
	close(IN);
    }
    return $fragment_length;
}

sub junctionTsdSelection{
    my $target = shift;
    my $file;
    opendir(DIR, "$wd/$target");
    foreach $file (sort readdir(DIR)){
	if($file =~ /junction_method.map.$target/){
	    &monitorWait;
	    chomp;
	    &log("tsd selection : $target : $file");
	    $cmd = "perl $0 a=$a,b=$b,sub=junctionTsdSelectionFunc,target=$target,file=$file &";
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : junctionTsdSelection : $cmd") if $rc;
	}
    }
    closedir(DIR);
}

sub junctionTsdSelectionFunc{
    my $s = {};
    ($head, $tail) = (split('\.', $file))[3,4];
    open(IN, "$wd/$target/$file");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[3] eq "head"){
	    $fl = substr($row[5], 0, 20);
	}else{
	    $fl = substr($row[5], 20);
	}
	$s->{$row[3]}{$row[1]}{$row[2]} = $fl;
    }
    close(IN);
    open(OUT, ">$wd/$target/junction_method.$target.$head.$tail");
    foreach $chr (sort bynumber keys %{$s->{head}}){
	foreach $pos (sort bynumber keys %{$s->{head}{$chr}}){
	    $hf = $s->{head}{$chr}{$pos};
	    for ($i = 1; $i <= 38; $i++){
		$tpos = $pos - 20 + $i;
		$tf = $s->{tail}{$chr}{$tpos};
		if ($tf ne ""){
		    if ($i < 20){
			$size = 20 - $i + 1;
			
		    }else{
			$size = $i - 20 + 1;
		    }
		    $htsd = substr($hf, 20 - $size, $size);
		    $ttsd = substr($tf, 0, $size);
		    if ($htsd eq $ttsd){
			next if $htsd =~ /^A+$/;
			next if $htsd =~ /^C+$/;
			next if $htsd =~ /^G+$/;
			next if $htsd =~ /^T+$/;
			next if $htsd =~ /^(AT)+$/;
			next if $htsd =~ /^(TA)+$/;
			next if $htsd =~ /^(TA)+T$/;
			next if $htsd =~ /^A(TA)+$/;
			next if $htsd =~ /^(AC)+$/;
			next if $htsd =~ /^(CA)+$/;
			next if $htsd =~ /^(CA)+C$/;
			next if $htsd =~ /^A(CA)+$/;
			next if $htsd =~ /^(AG)+$/;
			next if $htsd =~ /^(GA)+$/;
			next if $htsd =~ /^(GA)+G$/;
			next if $htsd =~ /^A(GA)+$/;
			next if $htsd =~ /^(TG)+$/;
			next if $htsd =~ /^(GT)+$/;
			next if $htsd =~ /^(GT)+G$/;
			next if $htsd =~ /^T(GT)+$/;
			next if $htsd =~ /^(TC)+$/;
			next if $htsd =~ /^(CT)+$/;
			next if $htsd =~ /^(CT)+C$/;
			next if $htsd =~ /^T(CT)+$/;
			if ($pos > $tpos){
			    $direction = "f";
			}else{
			    $direction = "r";
			}
			$s->{res}{$chr}{$pos} = "$target\t$chr\t$pos\t$tpos\t$direction\t$head\t$tail\t$htsd\t$hf\t$tf\n";
			print $s->{res}{$chr}{$pos};
			print OUT $s->{res}{$chr}{$pos};
		    }
		}
	    }
	}
    }
    close(OUT);
#    system("rm $wd/$target/$file");
}

sub junctionMapSelectionFunc{
    open(DAT, "grep $head $wd/$target/tmp/sorted.$chr |");
    open(OUT, "> $wd/$target/tmp/map_selected.$head.$tail.$chr");
    open(HEAD, "|sort | uniq > $wd/$target/tmp/pos.head.$head.$chr");
    while(<DAT>){
	chomp;
	@row = split;
	if($row[2] =~ /$head$/){
	    if ($row[5] eq "f"){
		$pos = $row[4] + 19;
	    }else{
		$pos = $row[4] - 19;
		}
	    $row[3] = "000$row[3]";
	    $row[3] = substr($row[3], length($row[3]) - 3, 3);
	    $pos = "00000000000$pos";
	    $pos = substr($pos, length($pos) - 11, 11);
	    $row[6] = "000$row[6]";
	    $row[6] = substr($row[6], length($row[6]) - 3, 3);
	    $row[7] = "00000000000$row[7]";
	    $row[7] = substr($row[7], length($row[7]) - 11, 11);
	    print OUT "$target\t$row[3]\t$pos\thead\t$head\t$row[2]\n";
	    print HEAD "$target\t$head\t$row[6]\t$row[7]\n";
	}
    }
    close(HEAD);
    open(DAT, "grep $tail $wd/$target/tmp/sorted.$chr |");
    open(TAIL, "|sort |uniq > $wd/$target/tmp/pos.tail.$tail.$chr");
    while(<DAT>){
	chomp;
	@row = split;
	if($row[2] =~ /^$tail/){
	    $row[3] = "000$row[3]";
	    $row[3] = substr($row[3], length($row[3]) - 3, 3);
	    $row[4] = "00000000000$row[4]";
	    $row[4] = substr($row[4], length($row[4]) - 11, 11);
	    $row[6] = "000$row[6]";
	    $row[6] = substr($row[6], length($row[6]) - 3, 3);
	    $row[7] = "00000000000$row[7]";
	    $row[7] = substr($row[7], length($row[7]) - 11, 11);
	    print OUT "$target\t$row[6]\t$row[7]\ttail\t$tail\t$row[2]\n";
	    print TAIL "$target\t$tail\t$row[3]\t$row[4]\n";
	}
    }
    close(DAT);
    close(OUT);
    close(TAIL);
}

sub junctionMapSelection{
    $target = shift;
    opendir(REF, "$wd/$ref");
    foreach $chr (sort readdir(REF)){
	if ($chr =~ /^chr/){
	    open(IN, "cat $wd/$a/tmp/te.candidate $wd/$b/tmp/te.candidate |");
	    while(<IN>){
		&monitorWait;
		chomp;
		($head, $tail) = (split)[0, 1];
		&log("map selection : $target : $chr $head $tail");
		$cmd = "perl $0 a=$a,b=$b,sub=junctionMapSelectionFunc,target=$target,chr=$chr,head=$head,tail=$tail &";
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionMapSelection : $cmd") if $rc;
	    }
	    close(IN);
	}
    }
    closedir(REF);
    &join;
    open(IN, "cat $wd/$a/tmp/te.candidate $wd/$b/tmp/te.candidate |");
    while(<IN>){
	chomp;
	($head, $tail) = (split)[0, 1];
	$s = {};
	open(DAT, "cat $wd/$target/tmp/map_selected.$head.$tail.* |");
	while(<DAT>){
	    @row = split;
	    $row[1] =~ s/^0+//;
	    $row[2] =~ s/^0+//;
	    $s->{j}{$row[1]}{$row[2]} = 1;
	}
	close(DAT);
	open(DAT, "cat $wd/$target/tmp/map_selected.$head.$tail.* | sort | uniq |");
	open(OUT, "> $wd/$target/junction_method.map.$target.$head.$tail");
	while(<DAT>){
	    chomp;
	    @row = split;
	    for($i = 0; $i <= $#row; $i++){
		$row[$i] =~ s/^0+//;
	    }
	    $flag = 0;
	    for($i = -10; $i <= 10; $i++){
		next if abs($i) < 2;
		$pos = $row[2] + $i;
		if ($s->{j}{$row[1]}{$pos}){
		    $flag = 1;
		}
	    }
	    $output = join("\t", @row);
	    print OUT "$output\n" if $flag;
	}
	close(DAT);
	close(OUT);
	open(DAT, "cat $wd/$target/tmp/pos.head.$head.* |");
	open(OUT, "|sort |uniq > $wd/$target/pos.$head.$tail");
	while(<DAT>){
	    chomp;
	    @row = split;
	    print OUT "$target\t$row[2]\t$row[3]\thead\t$row[1]\n";
	}
	open(DAT, "cat $wd/$target/tmp/pos.tail.$tail* |");
	while(<DAT>){
	    chomp;
	    @row = split;
	    print OUT "$target\t$row[2]\t$row[3]\ttail\t$row[1]\n";
	}
	close(DAT);
	close(OUT);
	open(DAT, "$wd/$target/pos.$head.$tail");
	open(OUT, ">$wd/$target/junction_method.tepos.$target.$head.$tail");
	while(<DAT>){
	    chomp;
	    @row = split;
	    $row[1] =~ s/^0+//g;
	    $row[2] =~ s/^0+//g;
	    print OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\n";
	}
	close(DAT);
	close(OUT);
	system("rm $wd/$target/pos.$head.$tail");
    }
    close(IN);
    system("rm $wd/$target/tmp/pos.head.$head.* $wd/$target/tmp/pos.tail.$tail.*");
}

sub junctionSelectCommon{
    open(IN, "$wd/$a/tmp/te.candidate");
    while(<IN>){
	chomp;
	@row = split;
	$head = shift(@row);
	$tail = shift(@row);
	$tsd = join(' ', @row);
	$ht = $head . $tail;
	$candidate{$ht} = $tsd;
    }
    close(IN);
    open(IN,"$wd/$b/tmp/te.candidate");
    open(OUT, "> $wd/$a/tmp/pair.$a.$b.junction_method");
    while(<IN>){
	chomp;
	@row = split;
	$head = shift(@row);
	$tail = shift(@row);
	$tsd = join(' ', @row);
	$ht = $head . $tail;
	if($candidate{$ht}){
	    %call = ();
	    %ca = ();
	    %cb = ();
	    $specific = 0;
	    foreach (split(' ', $candidate{$ht})){
		$ca{$_} = 1;
		$call{$_} = 1;
	    }
	    foreach (@row){
		$cb{$_} = 1;
		$call{$_} = 1;
	    }
	    foreach (keys %call){
		if (($ca{$_} and ! $cb{$_}) or (! $ca{$_} and $cb{$_})){
		    $specific ++;
		}
	    }
	    if ($specific >= 2){
		print "$head\t$tail\t$candidate{$ht}\t$tsd\n";
		print OUT "$head\t$tail\t$candidate{$ht}\t$tsd\n";
	    }
	}	
    }
    close(IN);
    close(OUT);
}

sub junctionSelectCandidate{
    $target = shift;
    $s = {};
    opendir(REF, "$wd/$ref");
    foreach $file (sort readdir(REF)){
	if ($file =~ /^chr/){
	    &monitorWait;
	    &log("select candidate : $target : $file");
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSelectCandidateFunc,chr=$file &";
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : junctionSelectCandidate : $cmd") if $rc;
	}
    }
    &join;
    open(IN, "cat $wd/$target/tmp/candidate.* |");
    while(<IN>){
	chomp;
	@row = split;
	$s->{$row[0]}{$row[1]}{$row[2]} = 1;
    }
    close(IN);
    open(OUT, "> $wd/$target/tmp/te.candidate");
    foreach $head (sort keys %$s){
	foreach $tail (sort keys %{$s->{$head}}){
	    $dat =  "$head\t$tail\t";
	    $count = 0;
	    foreach $tsd (sort keys %{$s->{$head}{$tail}}){
		next if $tsd =~ /^A+$/;
		next if $tsd =~ /^C+$/;
		next if $tsd =~ /^G+$/;
		next if $tsd =~ /^T+$/;
		next if $tsd =~ /^(AT)+$/;
		next if $tsd =~ /^(TA)+$/;
		next if $tsd =~ /^(TA)+T$/;
		next if $tsd =~ /^A(TA)+$/;
		next if $tsd =~ /^(AC)+$/;
		next if $tsd =~ /^(CA)+$/;
		next if $tsd =~ /^(CA)+C$/;
		next if $tsd =~ /^A(CA)+$/;
		next if $tsd =~ /^(AG)+$/;
		next if $tsd =~ /^(GA)+$/;
		next if $tsd =~ /^(GA)+G$/;
		next if $tsd =~ /^A(GA)+$/;
		next if $tsd =~ /^(TG)+$/;
		next if $tsd =~ /^(GT)+$/;
		next if $tsd =~ /^(GT)+G$/;
		next if $tsd =~ /^T(GT)+$/;
		next if $tsd =~ /^(TC)+$/;
		next if $tsd =~ /^(CT)+$/;
		next if $tsd =~ /^(CT)+C$/;
		next if $tsd =~ /^T(CT)+$/;
		$dat .= "$tsd\t";
		$count++;
	    }
	    chomp($dat);
	    print OUT "$dat\n" if $count > 1;
	}
    }
    close(OUT);
}

sub junctionSelectCandidateFunc{
    open(IN, "$wd/$target/tmp/sorted.$chr");
    open(OUT, "> $wd/$target/tmp/candidate.$chr");
    while(<IN>){
	chomp;
	@row = split;
	$tsd_size = 20 - ($row[0] - $prev[0]);
	if ($tsd_size <= 16 and $tsd_size >= 3){
	    if ($row[1] eq "h" and $prev[1] eq "t"){
		$htsd = substr($row[2], 20 - $tsd_size, $tsd_size);
		$ttsd = substr($prev[2], 20, $tsd_size);
		if ($htsd eq $ttsd){
		    $hflanking = substr($row[2], 0, 20);
		    $head = substr($row[2], 20);
		    $tail = substr($prev[2], 0, 20);
		    $tflanking = substr($prev[2], 21);
		    $hjunction = $row[0] + 19;
		    $tjunction = $row[0] + 20;
		    print OUT "$head\t$tail\t$htsd\n";
		}
	    }elsif ($row[1] eq "t" and $prev[1] eq "h"){
		$htsd = substr($prev[2], 20 - $tsd_size, $tsd_size);
		$ttsd = substr($row[2], 20, $tsd_size);
		if ($htsd eq $ttsd){
		    $hflanking = substr($pref[2], 0, 20);
		    $head = substr($prev[2], 20);
		    $tail = substr($row[2], 0, 20);
		    $tflanking = substr($row[2], 21);
		    $hjunction = $prev[0] + 19;
		    $tjunction = $prev[0] + 20;
		    print OUT "$head\t$tail\t$htsd\n";
		}
	    }
	}
	@prev = @row;
    }
    close(IN);
    close(OUT);
}

sub junctionSort{
    my $target = shift;
    my @chr;
    &log("Sorting filtered data : $target");
    opendir(REF, "$wd/$ref");
    foreach $file (sort readdir(REF)){
	if ($file =~ /^chr/){
	    push(@chr, $file);
	    open($file, "> $wd/$target/tmp/$file");
	}
    }
    open(IN, "cat $wd/$target/tmp/map.* |");
    my ($head, $tail);
    while(<IN>){
	chomp;
	@row = split;
	$head = "chr$row[1]";
	$tail = "chr$row[4]";
	print $head "$row[2]\th\t$_\n";
	print $tail "$row[5]\tt\t$_\n";
    }
    close(IN);
    foreach $file (@chr){
	close($file);
    }
    $org_processor = $processor;
    $processor = 2 if -s "$wd/$target/tmp/$chr[0]" > 10000000000;
    foreach $file (@chr){
	&monitorWait;
	&log("Sorting $file");
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSortFunc,chr=$file &";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : junctionSort : $cmd") if $rc;
    }
    $processor = $org_processor;
}

sub junctionSortFunc{
    $cmd = "sort -k 1 -n -S 1M -T $wd/$target/tmp $wd/$target/tmp/$chr > $wd/$target/tmp/sorted.$chr && rm $wd/$target/tmp/$chr";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : junctionSortFunc : $cmd") if $rc;
}

sub junctionSecondMapFunc{
    system("rm $wd/$target/tmp/map.$tag") if -e "$wd/$target/tmp/map.$tag";
    system("rm $wd/$target/tmp/tmp.$tag") if -e "$wd/$target/tmp/tmp.$tag";
    $cmd = "cat $wd/$target/tmp/first.$tag.* | sort -S 1M -T $wd/$target/tmp > $wd/$target/tmp/second.$tag";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : junctionSecondMapFunc : $cmd") if $rc;
    open(IN, "zcat $wd/$ref/ref20.$tag.gz |join $wd/$target/tmp/second.$tag - |");
    while(<IN>){
	chomp;
        @row = split;
	if ($prev ne $row[1]){
	    close($tag);
	    if (! $flag and -e "$wd/$target/tmp/tmp.$tag"){
		$cmd = "cat $wd/$target/tmp/tmp.$tag >> $wd/$target/tmp/map.$tag";
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSecondMapFunc : $cmd") if $rc;
	    }
	    $flag = 0;
	    open($tag, "> $wd/$target/tmp/tmp.$tag");
	}
	$prev = $row[1];
	print $tag "$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\n";
	if (abs($row[3] - $row[6]) <= 25){
	    $flag = 1;
	}
    }
    close(IN);
    if (! $flag and -e "$wd/$target/tmp/tmp.$tag"){
	$cmd = "cat $wd/$target/tmp/tmp.$tag >> $wd/$target/tmp/map.$tag";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : junctionSecondMapFunc : $cmd") if $rc;
    }
    system("rm $wd/$target/tmp/tmp.$tag");
    system("rm $wd/$target/tmp/first.$tag");
    system("rm $wd/$target/tmp/first.$tag.*");
    system("rm $wd/$target/tmp/second.$tag");
}

sub junctionFirstMapFunc{
    open(IN, "$wd/$target/tmp/specific.$tag");
    open(OUT, "| sort -S 1M -T $wd/$target/tmp > $wd/$target/tmp/target.$tag");
    while(<IN>){
	$seq = substr($_, 0, 20);
	next if $seq =~ /ATATATATATATATATAT/;
	next if $seq =~ /ACACACACACACACACAC/;
	next if $seq =~ /TGTGTGTGTGTGTGTGTG/;
	next if $seq =~ /AGAGAGAGAGAGAGAGAG/;
	next if $seq =~ /TCTCTCTCTCTCTCTCTC/;
	print OUT "$seq $_";
    }
    close(IN);
    close(OUT);
    open(IN, "zcat $wd/$ref/ref20.$tag.gz |join $wd/$target/tmp/target.$tag - |");
    open(OUT, "> $wd/$target/tmp/first.$tag");
    while(<IN>){
	@row = split;
	$seq = substr($row[1], 20, 20);
	next if $seq =~ /ATATATATATATATATAT/;
	next if $seq =~ /ACACACACACACACACAC/;
	next if $seq =~ /TGTGTGTGTGTGTGTGTG/;
	next if $seq =~ /AGAGAGAGAGAGAGAGAG/;
	next if $seq =~ /TCTCTCTCTCTCTCTCTC/;
	print OUT "$seq\t$row[1]\t$row[2]\t$row[3]\t$row[4]\n";
    }
    close(IN);
    close(OUT);
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $subtag = join('', @tag);
		open($subtag, "> $wd/$target/tmp/first.$subtag.$tag");
	    }
	}
    }
    open(IN, "$wd/$target/tmp/first.$tag");
    while(<IN>){
	$subtag = substr($_, 0, 3);
	print $subtag $_;
    }
    close(IN);
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $subtag = join('', @tag);
		close($subtag);
	    }
	}
    }
}

sub junctionSecondMap{
    my $target = shift;
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		&log("Mapping $target : $tag");
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSecondMapFunc,tag=$tag &";
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSecondMap : $cmd") if $rc;
            }
        }
    }
}

sub junctionFirstMap{
    my $target = shift;
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag &";
		&log($cmd);
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionFirstMap : $cmd") if $rc;
            }
        }
    }
}

sub junctionSpecificFunc{
    $cmd = "bash -c 'join -v 1 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz) > $wd/$a/tmp/specific.$tag.tmp'";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub : $cmd") if $rc;
    open(IN, "$wd/$a/tmp/specific.$tag.tmp");
    open(OUT, "> $wd/$a/tmp/specific.$tag");
    while(<IN>){
	chomp;
	@row = split;
	next if $row[0] =~/AAAAAAAAAA|CCCCCCCCCC|GGGGGGGGGG|TTTTTTTTTT/;
	if ($row[1] > 1){
	    print OUT "$row[0]\n";
	}
    }
    close(IN);
    close(OUT);
    $cmd = "bash -c 'join -v 2 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz)  > $wd/$b/tmp/specific.$tag.tmp'";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub : $cmd") if $rc;
    open(IN, "$wd/$b/tmp/specific.$tag.tmp");
    open(OUT, "> $wd/$b/tmp/specific.$tag");
    while(<IN>){
	chomp;
	@row = split;
	next if $row[0] =~/AAAAAAAAAA|CCCCCCCCCC|GGGGGGGGGG|TTTTTTTTTT/;
	if ($row[1] > 1){
	    print OUT "$row[0]\n";
	}
    }
    close(IN);
    close(OUT);
    system("rm wd/$a/tmp/specific.$tag.tmp wd/$b/tmp/specific.$tag.tmp");
}

sub junctionSpecific{
    my $target = shift;
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
                &monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tsd_size=$tsd_size,tag=$tag &";
                &log($cmd);
                $rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSpecific : $cmd") if $rc;
            }
        }
    }
}

sub mapFunc{
    open(IN, "zcat $ref/ref20.$tag | join $a/tmp/mapquery.$tag -|");
    open(OUT, "> $wd/$a/tmp/tsd_method.mapped.$a.$b.$tsd_size.$tag");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[2] eq "head"){
	    if ($row[6] eq "f"){
		$row[5] += 19;
	    }else{
		$row[5] -= 19;
	    }
	}
	print OUT "$row[1]\t$row[6]\t$row[4]\t$row[5]\t$row[2]\t$row[3]\t$row[0]\n";
    }
    close(IN);
    close(OUT);

}

sub map{
    &log("map");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		open($tag, "> $wd/$a/tmp/mapquery.$tag")
	    }
	}
    }

    open(IN, "sort $wd/$a/tmp/mapquery.tmp | uniq |");
    while(<IN>){
	$tag = substr($_, 0, 3);
	print $tag $_;
    }
    close(IN);

    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		close($tag)
	    }
	}
    }
    
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
                &monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mapFunc,tsd_size=$tsd_size,tag=$tag &";
                &log($cmd);
                $rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : mapFunc : $cmd") if $rc;
            }
        }
    }
    &join;
    open(IN, "cat $wd/$a/tmp/tsd_method.mapped.$a.$b.$tsd_size.* |sort |");
    while(<IN>){
	chomp;
	@row = split;
	if ($row[3] - $prev == $tsd_size - 1){
	    @prev = split('\t', $pline);
	    if ($row[1] eq "f"){
		$tsd = substr($row[6], 20 - $tsd_size, $tsd_size);
		$ttsd = substr($prev[6], 0, $tsd_size);
		if ($tsd eq $ttsd){
		    $s->{$row[0]}{$row[2]}{$row[3]}{$row[1]} = "$row[0]\t$row[2]\t$row[3]\t$prev[3]\t$row[1]\t$row[5]\t$prev[5]\t$tsd\t$row[6]\t$prev[6]";
		    $flanking{"$row[6] $prev[6]"} ++;
		    $ht{"$row[5] $prev[5]"} .= " $tsd";
		}
	    }else{
		$tsd = substr($prev[6], 20 - $tsd_size, $tsd_size);
		$ttsd = substr($row[6], 0, $tsd_size);
		if ($tsd eq $ttsd){
		    $s->{$row[0]}{$row[2]}{$row[3]}{$row[1]} = "$row[0]\t$row[2]\t$prev[3]\t$row[3]\t$row[1]\t$prev[5]\t$row[5]\t$tsd\t$prev[6]\t$row[6]";
		    $flanking{"$prev[6] $row[6]"} ++;
		    $ht{"$prev[5] $row[5]"} .= " $tsd";
		}
	    }
	}
	$pline = $_;
	$prev = $row[3];
    }
    close(IN);
    my $ht;
    foreach $ht (sort keys %ht){
	my ($tsd, %tsd, $total_tsd, $uniq_tsd);
	foreach $tsd (sort split('\ ', $ht{$ht})){
	    if ($tsd ne ""){
		$total_tsd ++;
		$tsd{$tsd} ++;
	    }
	}
	foreach (sort keys %tsd){
	    $uniq_tsd ++;
	}
	if ($uniq_tsd >= 3 or ($uniq_tsd / $total_tsd > 0.5 and $total_tsd > 1)){
	    $passed{$ht} = 1;
	    &log("## $ht $uniq_tsd / $total_tsd (unique/total TSD)");
	}
    }
    open(OUT, "> $wd/$a/tsd_method.result.$a.$b.$tsd_size");
    foreach $line (sort keys %$s){
        foreach $chr (sort bynumber keys %{$s->{$line}}){
            foreach $pos (sort bynumber keys %{$s->{$line}{$chr}}){
		foreach $direction (sort bynumber keys %{$s->{$line}{$chr}{$pos}}){
	   	    @row = split('\t', $s->{$line}{$chr}{$pos}{$direction});
		    if ($passed{"$row[5] $row[6]"}){
			if ($flanking{"$row[8] $row[9]"} == 1){
			    print $s->{$line}{$chr}{$pos}{$direction} . "\tuniq\n";
			    print OUT $s->{$line}{$chr}{$pos}{$direction} . "\tunique\n";
			}else{
			    print $s->{$line}{$chr}{$pos}{$direction} . "\n";
			    print OUT $s->{$line}{$chr}{$pos}{$direction} . "\n";
			}
		    }
		}
            }
	}
    }
    close(OUT);
}

sub searchQuery{
    if ($b eq ""){
	$cmd = "grep $seq $wd/$a/split/split.$number > $wd/$a/tmp/selected.$type.$seq.$number";
    }else{
	$cmd = "grep $seq $wd/$b/split/split.$number > $wd/$b/tmp/selected.$type.$seq.$number";
    }
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : searchQuery : $cmd") if $rc;
}

sub mkQuery{
    open(OUT, "> $wd/$a/tmp/mapquery.tmp");
    open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    while(<IN>){
	($head, $tail) = split;
	$head{$head} = 1;
	$tail{$tail} = 1;
    }
    close(IN);

    opendir(DIR, "$wd/$a/split");
    foreach (readdir(DIR)){
	if (/^split.[0-9]+$/){
	    $num = (split('\.', $_))[1];
	    if ($last < $num){
		$last = $num;
	    }
	}
    }
    foreach ($i = 1; $i <= $last; $i++){
	foreach $head (sort keys %head){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,sub=searchQuery,number=$i,type=head,seq=$head &";
	    &log($cmd);
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : $sub : $cmd") if $rc;
	}
        foreach $tail (sort keys %tail){
            &monitorWait;
            $cmd = "perl $0 a=$a,sub=searchQuery,number=$i,type=tail,seq=$tail &";
            &log($cmd);
            $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : $sub : $cmd") if $rc;
        }
	&join;
	opendir(DIR, "$wd/$a/tmp");
	foreach $file (readdir(DIR)){
	    if ($file =~ /^sele/){
		($type, $seq) = (split('\.', $file))[1,2];
		if (-s "$wd/$a/tmp/$file" > 0){
		    open(IN, "$wd/$a/tmp/$file");
		    while(<IN>){
			chomp;
			$pos = index($_, $seq);
			if ($type eq "head"){
			    $flanking = substr($_, $pos - 20, 20);
			}else{
			    $flanking = substr($_, $pos + 20, 20);
			}
			if (length($flanking) == 20 and $flanking !~ /N/){
			    print OUT "$flanking $a $type $seq\n";
			}
		    }
		    close(IN);
		}
		system("rm $wd/$a/tmp/$file");
	    }
	}
	closedir(DIR);
    }
    $last = 0;
    opendir(DIR, "$wd/$b/split");
    foreach (readdir(DIR)){
	if (/^split.[0-9]+$/){
	    $num = (split('\.', $_))[1];
	    if ($last < $num){
		$last = $num;
	    }
	}
    }
    foreach ($i = 1; $i <= $last; $i++){
	foreach $head (sort keys %head){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,sub=searchQuery,number=$i,type=head,seq=$head &";
	    &log($cmd);
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : $sub : $cmd") if $rc;
	}
	foreach $tail (sort keys %tail){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,sub=searchQuery,number=$i,type=tail,seq=$tail &";
	    &log($cmd);
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : $sub : $cmd") if $rc;
	}
	&join;
	opendir(DIR, "$wd/$b/tmp");
	foreach $file (readdir(DIR)){
	    if ($file =~ /^sele/){
		($type, $seq) = (split('\.', $file))[1,2];
		if (-s "$wd/$b/tmp/$file" > 0){
		    open(IN, "$wd/$b/tmp/$file");
		    while(<IN>){
			chomp;
			$pos = index($_, $seq);
			if ($type eq "head"){
			    $flanking = substr($_, $pos - 20, 20);
			}else{
			    $flanking = substr($_, $pos + 20, 20);
			}
			if (length($flanking) == 20 and $flanking !~ /N/){
			    print OUT "$flanking $b $type $seq\n";
			}
		    }
		    close(IN);
		}
		system("rm $wd/$b/tmp/$file");
	    }
	}
	closedir(DIR);
    }
    close(OUT);
}

sub mkref{
    if (! -d "$wd/$ref"){
	system("mkdir $ref");
	print "
$ref directory has been created.
Save gzipped fasta file of reference genome into the $ref directory.
The fasta files should be compressed and .gz extention is required.

After saving the gz file, run $0 again with same options.

";
	exit;
    }
    chdir "$wd/$ref";
    if (! -e "config"){
	&log("making config file of in $ref directory.");
	$cmd = "zcat *.gz|grep '>' | sed  -e 's/^>//' > config.tmp";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : mkref : $cmd") if $rc;
	open(IN, "config.tmp");
	open(OUT, "> config");
	while(<IN>){
	    print;
	    chomp;
	    s/^chr//i;
	    $_ += 0 if /^[0-9]$/;
	    @row = split;
	    print OUT "$row[0]\n";
	}
	close(IN);
	close(OUT);
	system("rm config.tmp");
	print "
 Edit config file in $ref directory.

 Replace lines to 'NOP' for the chromosome entry which is not required.

 And run again $0 script.

";
	exit;
    }
    if (! -e "$wd/$ref/ref20.TTT.gz"){
	open(IN, "config");
	while(<IN>){
	    chomp;
	    push(@chr, $_);
	}
	close(IN);
	open(IN, "zcat *.gz|");
	while(<IN>){
	    if (/^>/){
		$chr++;
		$chr_name = $chr[$chr -1];
		close(CHR);
		if ($chr_name ne "NOP" and $chr_name ne ""){
		    &log("Making chr$chr_name file");
		    $chr_name =~ s/^chr//i;
		    open(CHR, "> chr$chr_name");
		}
	    }elsif($chr_name ne "NOP" and $chr_name ne ""){
		chomp;
		y/a-z/A-Z/;
		y/BDEFHIJKLMOPQRSUVWXYZ/N/;
		print CHR;
	    }
	}
	close(IN);
	close(CHR);
	chdir $wd;
	&mk20;
    }
    chdir $wd;
}

sub sort20mer{
    $cmd = "sort -S 1M -T $wd/$ref/tmp $wd/$ref/tmp/ref20.$tag.* |gzip -c > $wd/$ref/ref20.$tag.gz";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : sort20mer : $cmd") if $rc;
    system("rm $wd/$ref/tmp/ref20.$tag.*");
}

sub mk20mer{
    my ($fin, $fout, $i, $nuc, @tag, $tag, $forward, $length, $comp, $fpos, $rpos, $notstd);
    return if $chr eq "N";
    
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$fout = "$tag.$chr";
		open($fout, "> $wd/$ref/tmp/ref20.$tag.$chr");
	    }
	}
    }

    my $file = "$wd/$ref/chr$chr";
    $fin = "fin.$chr";
    open($fin, $file);
    binmode($fin);
    $length = -s $file;
    for ($i = 0; $i <= $length - 20; $i++){
	seek($fin, $i, 0);
	read($fin, $forward, 20);
	if ($forward !~ /N/){
	    $comp = complement($forward);
	    $fpos = $i + 1;
	    $rpos = $fpos + 20 -1;
	    $tag = substr($forward, 0, 3) . ".$chr";
	    print $tag "$forward\t$chr\t$fpos\tf\n";
	    $tag = substr($comp, 0, 3) . ".$chr";
	    print $tag "$comp\t$chr\t$rpos\tr\n";
	}
    }
    close($fin);
    
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$fout = "$tag.$chr";
		close($fout);
	    }
	}
    }
}

sub mk20{
    system("mkdir $wd/$ref/tmp") if ! -d "$wd/$ref/tmp";
    if (! -e "$wd/$ref/ref20.TTT.gz"){
	&log("Making 20mer position file.");
	foreach $i (@chr){
	    next if $i eq "NOP";
	    &monitorWait;
	    &log("Processing chr$i");
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mk20mer,chr=$i &";
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : mk20 : $cmd") if $rc;
	}
	&join;
    }
    my (@row, %tag);
    opendir(DIR, "$wd/$ref/tmp/");
    foreach (readdir(DIR)){
	if (/^ref20/){
	    @row = split('\.', $_);
	    $tag{$row[1]} = 1;
	}
    }
    closedir(DIR);
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		if ($tag{$tag}){
		    &monitorWait;
		    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=sort20mer,tag=$tag &";
		    &log($cmd);
		    $rc = system($cmd);
		    $rc = $rc >> 8;
		    &log("ERROR : verify : $cmd") if $rc;
		}
	    }
	}
    }    
    &join;
    system("rm -r $wd/$ref/tmp")
}

sub verify2{
    opendir(DIR, "$a/split");
    foreach (sort readdir(DIR)){
	if (/^split/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$a,sub=verifyFunc,tsd_size=$tsd_size,file=$_ &";
	    &log("$cmd");
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : verify :$cmd") if $rc;
	}
    }
    closedir(DIR);
    opendir(DIR, "$b/split");
    foreach (sort readdir(DIR)){
	if (/^split/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$b,sub=verifyFunc,tsd_size=$tsd_size,file=$_ &";
	    &log($cmd);
	    $rc = system($cmd);
	    $rc = $rc >> 8;
	    &log("ERROR : verify : $cmd") if $rc;
	}
    }
    closedir(DIR);
    &join;
    opendir(DIR, "$a/tmp");
    foreach (sort readdir(DIR)){
	if (/verify/){
	    open(IN, "$a/tmp/$_");
	    while(<IN>){
		chomp;
		@row = split;
		$s->{$row[0]}{$row[1]}{$row[2]} = $row[3];
	    }
	    close(IN);
	}
    }
    foreach (sort readdir(DIR)){
	if (/verify/){
	    open(IN, "$b/tmp/$_");
	    while(<IN>){
		chomp;
		@row = split;
		$s->{$row[0]}{$row[1]}{$row[2]} = $row[3];
	    }
	close(IN);
	}
    }
    
    open(OUT, "|sort -T $wd/$a/tmp |uniq > $wd/$a/tmp/query");
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/tmp/te.candidate");
    }
    while(<IN>){
	print;
	@row = split;
	foreach $upstream (sort keys %{$s->{$row[0]}{H}}){
	    $htsd = $s->{$row[0]}{H}{$upstream};
	    foreach $downstream (sort keys %{$s->{$row[1]}{T}}){
		$ttsd = $s->{$row[1]}{T}{$downstream};
		if ($htsd eq $ttsd){
		    $tsdsize = length($htsd);
		    $wildtype = $upstream . substr($downstream, $tsdsize, 20);
		    $head = $upstream . $row[0];
		    $tail = $row[1] . substr($downstream, 0, 20);
		    print OUT "$head\n";
		    print OUT "$tail\n";
		    print OUT "$wildtype\n";
		}
	    }
	}
    }
    close(IN);
    close(OUT);
    $cmd = "zcat $wd/$a/count.20/*.gz | join - $wd/$a/tmp/query > $wd/$a/tmp/verify.count";
    &log($cmd);
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify :$cmd") if $rc;
    
    $cmd = "zcat $wd/$b/count.20/*.gz | join - $wd/$a/tmp/query > $wd/$b/tmp/verify.count";
    &log($cmd);
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify : $cmd") if $rc;
    
    open(IN, "$wd/$a/tmp/verify.count");
    while(<IN>){
	chomp;
	@row = split;
	$a{$row[0]} = $row[1];
	
    }
    open(IN, "$wd/$b/tmp/verify.count");
    while(<IN>){
	chomp;
	@row = split;
	$b{$row[0]} = $row[1];
    }
    
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
	open(OUT, "> $wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/tmp/te.candidate");
	open(OUT, "> $wd/$a/junction_method.verify.$a.$b");
    }
    while(<IN>){
	@row = split;
	foreach $upstream (sort keys %{$s->{$row[0]}{H}}){
	    $htsd = $s->{$row[0]}{H}{$upstream};
	    foreach $downstream (sort keys %{$s->{$row[1]}{T}}){
		$ttsd = $s->{$row[1]}{T}{$downstream};
		if ($htsd eq $ttsd){
		    $tsdsize = length($htsd);
		    $wildtype = $upstream . substr($downstream, $tsdsize, 20);
		    $head = $upstream . $row[0];
		    $downfl = substr($downstream, 0, 20);
		    $tail = $row[1] . $downfl;
		    if (($a{$head} == 0 and $a{$tail} == 0 and $a{$wildtype} > 0 and $b{$head} > 0 and $b{$tail} > 0) or ($b{$head} == 0 and $b{$tail} == 0 and $a{$wildtype} > 0 and $a{$head} > 0 and $a{$tail} > 0)){
			$result = "$row[0]\t$row[1]\t$upstream\t$downfl\t$a{$head}\t$a{$tail}\t$a{$wildtype}\t$b{$head}\t$b{$tail}\t$b{$wildtype}\n";
			print $result;
			print OUT $result;
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);
}

sub verifyFunc{
    my (@row, $tsdsize, %seq, %tsdsize, $length, $seq, $tsd, $upstream, $downstream);
    print "$target $file\n";
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/tmp/te.candidate");
    }
    while(<IN>){
	@row = split;
	$tsdsize = length($row[2]);
	$seq{$row[0]} = "H";
	$seq{$row[1]} = "T";
	$tsdsize{$row[0]} = $tsdsize;
	$tsdsize{$row[1]} = $tsdsize;
    }
    close(IN);
    open(OUT, "> $wd/$target/tmp/verify.$file");
    open(IN, "$wd/$target/split/$file");
    while(<IN>){
	chomp;
	$length = length($_) if $length eq "";
	for($i = 0; $i <= $length - 20; $i++){
	    $seq = substr($_, $i, 20);
	    if ($seq{$seq}){
		$tsdsize = $tsdsize{$seq};
		if ($seq{$seq} eq "H"){
		    $upstream = substr($_, $i - 20, 20);
		    if(length($upstream) == 20){
			$tsd = substr($_, $i - $tsdsize, $tsdsize);
			print OUT "$seq\tH\t$upstream\t$tsd\n";
		    }
		}elsif($seq{$seq} eq "T"){
		    $downstream = substr($_, $i + 20, 20 + $tsdsize);
		    if(length($downstream) == 20 + $tsdsize){
			$tsd = substr($_, $i +20, $tsdsize);
			print OUT "$seq\tT\t$downstream\t$tsd\n";
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);
}

sub verify{
    &log("Verify");
    open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    while(<IN>){
	@row = split;
	$seq{$row[0]} = "H";
	$seq{$row[1]} = "T";
    }
    close(IN);
    
    $command = "cat";
    opendir(DIR, "$wd/$a/read");
    foreach (sort readdir(DIR)){
	if (/gz$/){
	    $command = "zcat";
	    last;
	}
	if (/bz2$/){
	    $command = "bzcat";
	    last;
	}
	if (/xz$/){
	    $command =  "xzcat";
	    last;
	}
    }
    close(DIR);
    
    $count = 0;
    open(IN, "$command $wd/$a/read/* |");
    while(<IN>){
	$count++;
	if ($count % 4 == 2){
	    next if /N/;
	    chomp;
	    $length = length($_);
	    for($i = $tsd_size; $i <= $length -20; $i++){
		$kmer = substr($_, $i, 20);
		$utsd = substr($_, $i - $tsd_size, $tsd_size);
		$dtsd = substr($_, $i + 20, $tsd_size);
		if ($seq{$kmer} eq "H"){
		    $upstream = substr($_, $i - 20, 20);
		    next if length($upstream) != 20;
		    $s->{$utsd}{H}{$upstream}{$kmer} ++;
		}elsif($seq{$kmer} eq "T"){
		    $downstream = substr($_, $i + 20, 20);
		    next if length($downstream) != 20;
		    $s->{$dtsd}{T}{$downstream}{$kmer} ++;
		}
	    }
	}
    }
    close(IN);
    &log("Verify. output tsd.$a.$b.$tsd_size");
    open(TSD, "> $wd/$a/tsd.$a.$b.$tsd_size");
    foreach $tsd (sort keys %{$s}){
	foreach $type (sort keys %{$s->{$tsd}}){
	    foreach $flanking (sort keys %{$s->{$tsd}{$type}}){
		foreach $kmer (sort keys %{$s->{$tsd}{$type}{$flanking}}){
		    if ($s->{$tsd}{$type}{$flanking}{$kmer} > 3){
			if ($tsd eq $ptsd and $type eq "T"){
			    $flanking = substr($flanking, $tsd_size);
			    $target = $pflanking . $flanking;
			    print TSD "$tsd $pkmer $kmer $target\n";
			    $target{$target} = "$pkmer\t$kmer";
			}
			$ptsd = $tsd;
			$ptype = $type;
			$pflanking = $flanking;
			$pkmer = $kmer;
		    }
		}
	    }
	}
    }
    close(TSD);
    $command = "cat";
    opendir(DIR, "$wd/$b/read");
    foreach (sort readdir(DIR)){
	if (/gz$/){
	    $command = "zcat";
	    last;
	}
	if (/bz2$/){
	    $command = "bzcat";
	    last;
	}
	if (/xz$/){
	    $command =  "xzcat";
	    last;
	}
    }
    close(DIR);
    $count = 0;
    &log("Verify. output transposon.$a.$b.$tsd_size");
    open(OUT, "> $wd/$a/tsd_method.transposon.$a.$b.$tsd_size");
    open(IN, "$command $wd/$b/read/* |");
    while(<IN>){
	$count++;
	if ($count % 4 == 2){
	    next if /N/;
	    for($i = 0; $i < length($_) - 40 + $tsd_size; $i++){
		my $seq = substr($_, $i, 40 - $tsd_size);
		if (length($seq) == 40 -$tsd_size and defined $target{$seq}){
		    print "$target{$seq}\t$seq\t$_";
		    print OUT "$target{$seq}\t$seq\t$_";
		}
	    }
	}
    }
}

sub tsd{
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$aout = "A-$tag";
		$bout = "B-$tag";
		open($aout, "| sort -S 1M -T $wd/$a/tmp | uniq > $wd/$a/tmp/tsd.head.$tag");
		open($bout, "| sort -S 1M -T $wd/$b/tmp | uniq > $wd/$b/tmp/tsd.head.$tag");
	    }
	}
    }
    open(IN, "join $wd/$a/tmp/head.select $wd/$b/tmp/head.select|");
    while(<IN>){
	chomp;
	@row = split;
	$tag = substr($row[1], 0, 3);
	$aout = "A-$tag";
	print $aout "$row[1]\t$row[0]\t$row[2]\n";
	$tag = substr($row[3], 0, 3);
	$bout = "B-$tag";
	print $bout "$row[3]\t$row[0]\t$row[4]\n";
    }
    close(IN);
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$aout = "A-$tag";
		$bout = "B-$tag";
		close($aout);
		close($bout);
	    }
	}
    }
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$aout = "A-$tag";
		$bout = "B-$tag";
		open($aout, "| sort -S 1M -T $wd/$a/tmp | uniq > $wd/$a/tmp/tsd.tail.$tag");
		open($bout, "| sort -S 1M -T $wd/$b/tmp | uniq > $wd/$b/tmp/tsd.tail.$tag");
	    }
	}
    }
    open(IN, "join $wd/$a/tmp/tail.select $wd/$b/tmp/tail.select|");
    while(<IN>){
	chomp;
	@row = split;
	$tag = substr($row[1], 0, 3);
	$aout = "A-$tag";
	print $aout "$row[1]\t$row[0]\t$row[2]\n";
	$tag = substr($row[3], 0, 3);
	$bout = "B-$tag";
	print $bout "$row[3]\t$row[0]\t$row[4]\n";
    }
    close(IN);
    foreach $nuc (@nuc){
	$tag[0] = $nuc;
	foreach $nuc (@nuc){
	    $tag[1] = $nuc;
	    foreach $nuc (@nuc){
		$tag[2] = $nuc;
		$tag = join('', @tag);
		$aout = "A-$tag";
		$bout = "B-$tag";
		close($aout);
		close($bout);
	    }
	}
    }
}

sub pair{
    open(OUT, "> $wd/$a/tmp/p.$tag");
    open(IN, "join $wd/$a/tmp/tsd.head.$tag $wd/$a/tmp/tsd.tail.$tag|");
    while(<IN>){
	chomp;
	@row = split;
	print OUT "$row[1]" . "_" . "$row[3]\t$row[0]\t$row[2]\t$row[4]\n";
    }
    close(IN);
    close(OUT);
    system("rm $wd/$a/tmp/tsd.*.$tag");

    open(OUT, "> $wd/$b/tmp/p.$tag");
    open(IN, "join $wd/$b/tmp/tsd.head.$tag $wd/$b/tmp/tsd.tail.$tag|");
    while(<IN>){
	chomp;
	@row = split;
	print OUT "$row[1]" . "_" . "$row[3]\t$row[0]\t$row[2]\t$row[4]\n";
    }
    close(IN);
    close(OUT);
    system("rm $wd/$b/tmp/tsd.*.$tag");
}

sub sortPair{
   foreach $taga (@nuc){
       foreach $tagb (@nuc){
	   foreach $tagc (@nuc){
	       $tag = $taga . $tagb . $tagc;
	       open($tag, "> $wd/$a/tmp/pair.$tag.tmp");
	   }
       }
   }
   open(IN, "cat $wd/$a/tmp/p.*|");
   while(<IN>){
       $tag = substr($_, 0, 3);
       print $tag $_;
   }
   close(IN);
   foreach $taga (@nuc){
       foreach $tagb (@nuc){
	   foreach $tagc (@nuc){
	       $tag = $taga . $tagb . $tagc;
	       close($tag);
	   }
       }
   }

   foreach $taga (@nuc){
       foreach $tagb (@nuc){
           foreach $tagc (@nuc){
               $tag = $taga . $tagb . $tagc;
	       &monitorWait;
	       &log("making ref20.$tag.gz");
	       $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=sortPairFunc,tag=$tag,target=$a &";
	       $rc = system($cmd);
	       $rc = $rc >> 8;
	       &log("ERROR : $sub : $cmd") if $rc;
	   }
       }
   }

   foreach $taga (@nuc){
       foreach $tagb (@nuc){
           foreach $tagc (@nuc){
               $tag = $taga . $tagb . $tagc;
               open($tag, "> $wd/$b/tmp/pair.$tag.tmp");
           }
       }
   }
   open(IN, "cat $wd/$b/tmp/p.*|");
   while(<IN>){
       $tag = substr($_, 0, 3);
       print $tag $_;
   }
   close(IN);
   foreach $taga (@nuc){
       foreach $tagb (@nuc){
	   foreach $tagc (@nuc){
	       $tag = $taga . $tagb . $tagc;
	       close($tag);
	       &monitorWait;
	       &log("sort pair.$tag in $b");
	       $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=sortPairFunc,tag=$tag,target=$b &";
               $rc = system($cmd);
	       $rc = $rc >> 8;
	       &log("ERROR : $sub : $cmd") if $rc;
	   }
       }
   }
   system("rm $wd/$a/tmp/p.* $wd/$b/tmp/p.*");
}

sub sortPairFunc{
    $cmd = "sort -S 1M -T $wd/$target/tmp $wd/$target/tmp/pair.$tag.tmp > $wd/$target/tmp/pair.$tag && rm $wd/$target/tmp/pair.$tag.tmp";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : sortPairFunc : $cmd") if $rc;
}

sub selectPair{
    open(IN, "join $wd/$a/tmp/pair.$tag $wd/$b/tmp/pair.$tag |");
    open(OUT, "> $wd/$a/tmp/pair.$a.$b.$tsd_size.$tag");
    while(<IN>){
	chomp;
	@row = split;
	if ($prev ne $row[0]){
	    foreach (sort keys %c){
		if ($a{$_}){
		    $outa .=  "$_ ";
		    $counta ++;
		    $total ++;
		}
		if ($b{$_}){
		    $outb .=  "$_ ";
		    $countb ++;
		    $total ++;
		}
		if ($a{$_} and $b{$_}){
		    $countc ++;
		    $total ++;
		}
	    }

	    if (($counta > 2 and $countb > 2) and $countc / $total < 0.3 and $counta / $countah > $th and $counta / $countat > $th and $countb / $countbh > $th and $countb / $countbt > $th){
		($head, $tail) = split('_', $prev);
		($ha = $head) =~ y/A//cd;
		($hc = $head) =~ y/C//cd;
		($hg = $head) =~ y/G//cd;
		($ht = $head) =~ y/T//cd;
		($ta = $tail) =~ y/A//cd;
		($tc = $tail) =~ y/C//cd;
		($tg = $tail) =~ y/G//cd;
		($tt = $tail) =~ y/T//cd;
		$hecount = 0;
		$tecount = 0;
		$hecount ++ if length($ha) < 2;
		$hecount ++ if length($hc) < 2;
		$hecount ++ if length($hg) < 2;
		$hecount ++ if length($ht) < 2;
		$tecount ++ if length($ta) < 2;
		$tecount ++ if length($tc) < 2;
		$tecount ++ if length($tg) < 2;
		$tecount ++ if length($tt) < 2;

		if ($hecount <= 1 and $tecount <= 1){ 
		    $output = "$head\t$tail\t" . $outa . "\t" . $outb . "\n";
		    print $output;
		    print OUT $output;
		}
	    }
	    %a = ();
	    %b = ();
	    %c = ();
	    $outa = "";
	    $outb = "";
	    $counta = 0;
	    $countb = 0;
	    $countc = 0;
	    $total = 0;
	}
	$prev = $row[0];
	$a{$row[1]} = 1;
	$b{$row[4]} = 1;
	$c{$row[1]} = 1;
	$c{$row[4]} = 1;
	$countah = $row[2];
	$countat = $row[3];
	$countbh = $row[5];
	$countbt = $row[6];
    }
    close(IN);
    close(OUT);
    system("rm $wd/$a/tmp/pair.$tag $wd/$b/tmp/pair.$tag");
}

sub tsdJoin{
    &monitorWait;
    $cmd = "perl $0 target=$a,sub=countHeadTail,type=head,a=$a &";
    &log($cmd);
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub :$cmd") if $rc;
    &monitorWait;
    $cmd = "perl $0 target=$a,sub=countHeadTail,type=tail,a=$a &";
    &log($cmd);
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub : $cmd") if $rc;
    &monitorWait;
    $cmd = "perl $0 target=$b,sub=countHeadTail,type=head,a=$a &";
    &log($cmd);
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub : $cmd") if $rc;
    &monitorWait;
    $cmd = "perl $0 target=$b,sub=countHeadTail,type=tail,a=$a &";
    &log($cmd);
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : $sub : $cmd") if $rc;
}

sub countHeadTail{
    my (@row, $prev, $count);
    $cmd = "cat $wd/$target/tmp/*.$type | sort -S 1M -T $wd/$target/tmp > $wd/$target/tmp/$type";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : countHeadTail : $cmd") if $rc;
    open(OUT, "> $wd/$target/tmp/$type.count");
    open(IN, "$wd/$target/tmp/$type");
    while(<IN>){
	@row = split;
	if ($prev ne $row[0] and $prev ne ""){
	    print OUT "$prev\t$count\n";
	    $count = 0;
	}
	$prev = $row[0];
	$count++;
    }
    close(IN);
    print OUT "$prev\t$count\n";
    close(OUT); 
    open(OUT, "> $wd/$target/tmp/$type.select");
    open(IN, "join $wd/$target/tmp/$type $wd/$target/tmp/$type.count|");
    while(<IN>){
	@row = split;
	if ($row[2] > 1){
	    print OUT "$_";
	}
    }
    close(IN);
    close(OUT); 
}

sub commFunc{
    open(AH, "> $wd/$a/tmp/$tag.head");
    open(BH, "> $wd/$b/tmp/$tag.head");
    open(AT, "> $wd/$a/tmp/$tag.tail");
    open(BT, "> $wd/$b/tmp/$tag.tail");
    open(IN, "bash -c 'comm -3 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz) '|");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
	next if ($row[0] ne "" and  $row[1] <= 4) or ($row[0] eq "" and $row[2] <= 4);
	next if ($row[0] =~ /A{8,}/ or $row[1] =~ /A{8,}/);
	next if ($row[0] =~ /C{8,}/ or $row[1] =~ /C{8,}/);
	next if ($row[0] =~ /G{8,}/ or $row[1] =~ /G{8,}/);
	next if ($row[0] =~ /T{8,}/ or $row[1] =~ /T{8,}/);
	if ($row[0] ne ""){
	    $head = substr($row[0], $tsd_size);
	    $tsd  = substr($row[0], 0, $tsd_size);
	    print AH "$head\t$tsd\n";
	    $tail = substr($row[0], 0, 20);
	    $tsd = substr($row[0], 20, $tsd_size);
	    print AT "$tail\t$tsd\n";
	}
	if ($row[1] ne ""){
	    $head = substr($row[1], $tsd_size);
	    $tsd  = substr($row[1], 0, $tsd_size);
	    print BH "$head\t$tsd\n";
	    $tail = substr($row[1], 0, 20);
	    $tsd = substr($row[1], 20, $tsd_size);
	    print BT "$tail\t$tsd\n";
	}
    }
    close(AH);
    close(BH);
    close(AT);
    close(BT);
}

sub comm{
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,sub=commFunc,tag=$tag,tsd_size=$tsd_size &";
		&log("common : $tag : $tsd_size");
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : comm : $cmd") if $rc;
	    }
	}
    }
}

sub mergeFunc{
    my ($num, @num, @sorted, $a, $b, $c, $filea, $fileb, $filec, @row);

    opendir(DIR, "$wd/$target/count.$tsd_size/");
    foreach (readdir(DIR)){
	if (/^$tag.[0-9]+/){
	    ($num = $_) =~ y/0-9//cd;
	    push(@num, $num);
	}
    }
    @sorted = sort bynumber @num;
    
    
    $a = $sorted[0];
    $b = $a + 1;
    $c = $sorted[$#sorted];
    $c = $c + 1;
    while(1){
	last if $a == $c or $b == $c;
	$filea = "$tag.$a.gz";
	$fileb = "$tag.$b.gz";
	$filec = "$tag.$c.gz";
	$cmd = "bash -c 'join -a 1 -a 2 <(zcat $wd/$target/count.$tsd_size/$filea) <(zcat $wd/$target/count.$tsd_size/$fileb)' | awk '{print \$1 \"\t\" \$2 + \$3}' |gzip -c > $wd/$target/count.$tsd_size/$filec && rm $wd/$target/count.$tsd_size/$filea $wd/$target/count.$tsd_size/$fileb";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : mergeFunc : $cmd") if $rc;
	$a += 2;
	$b += 2;
	$c ++;
    }

    open(IN, "zcat $wd/$target/count.$tsd_size/$filec |");
    open(OUT, "|gzip -c > $wd/$target/count.$tsd_size/$tag.gz");
    while(<IN>){
	@row = split;
	print OUT "$row[0]\t$row[1]\n";
    }
    close(IN);
    close(OUT);
    system("rm $wd/$target/count.$tsd_size/$tag.*.gz");
}

sub merge{
    my $target = shift;
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		$cmd = "perl $0 target=$target,sub=mergeFunc,tsd_size=$tsd_size,tag=$tag,a=$a &";
		&log($cmd);
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : merge :$cmd") if $rc;
	    }
	}
    }
}

sub countFunc{
    my ($taga, $tagb, $tagc, $tag, $count, $nuc);
    open(IN, "$wd/$target/split/split.$number");
    open(SEQ, "|sort -S 1M -T $wd/$target/tmp |uniq -c |awk '{print \$2 \"\t\" \$1}' > $wd/$target/count.$tsd_size/count.$number");
    while(<IN>){
	    chomp;
	    &mkKmer($_);
	    &mkKmer(&complement($_));
    }
    close(IN);
    close(SEQ);
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		$fout = $tag . ".out";
		open($fout, "|gzip -c > $wd/$target/count.$tsd_size/$tag.$number.gz");
	    }
	}
    }
    
    open(IN, "$wd/$target/count.$tsd_size/count.$number");
    while(<IN>){
	$tag = substr($_, 0, 3);
	$fout = $tag . ".out";
	print $fout "$_";
    }
    close(IN);
    sleep 10;
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		$fout = $tag . ".out";
		close($fout);
	    }
	}
    }
    system("rm $wd/$target/count.$tsd_size/count.$number");
}

sub mkKmer{
    my $read = shift;
    my $length = length($read);
    my ($i, $seq);
    for($i = 0; $i < $length; $i++){
        $seq = substr($read, $i, 20 + $tsd_size);
        last if length($seq) != 20 + $tsd_size;
        next if $seq =~ /N/;
        print SEQ "$seq\n";
    }
}

sub count{
    my $target = shift;
    my ($last, $i, $cmd);
    opendir(DIR, "$wd/$target/split");
    foreach (readdir(DIR)){
	if (/^split/){
	    $number = (split('\.', $_))[1];
	    if ($last < $number){
		$last = $number;
	    }
	}
    }
    system("mkdir $wd/$target/count.$tsd_size") if ! -d "$wd/$target/count.$tsd_size";
    foreach ($i = 1; $i <= $last; $i++){
	&monitorWait;
	$cmd = "perl $0 target=$target,sub=countFunc,number=$i,a=$a,tsd_size=$tsd_size &";
	&log($cmd);
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : count : $cmd") if $rc;
    }
}

sub split{
    my $number = 1;
    my $command = "cat";
    opendir(DIR, "$wd/$target/read");
    foreach (sort readdir(DIR)){
	if (/gz$/){
	    $command = "zcat";
	    last;
	}
	if (/bz2$/){
	    $command = "bzcat";
	    last;
	}
	if (/xz$/){
	    $command =  "xzcat";
	    last;
	}
    }
    close(DIR);
    open(IN, "$command $wd/$target/read/* |");
    open(OUT, "> $wd/$target/split/split.$number");
    while(<IN>){
	$count++;
	if($count == 2){
		$lines++;
		print OUT;
	}elsif($count == 4){
	    $count = 0;
	}
	if ($lines == 1000000){
	    close(OUT);
	    $number++;
	    open(OUT, "> $target/split/split.$number");
	    $lines = 0;
	}
    }
    close(IN);
    close(OUT);
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = length($seq);
    my $out = "";
    while($i > 0){
        $i--;
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }else{
            $out .= "N";
        }
    }
    return $out;
}

sub monitorWait{
    while(1){
        my $count = 0;
        $wait_time = 0.1 if $wait_time eq "";
        select undef, undef, undef, $wait_time;
        opendir(CDIR, "$wd/$a/child");
        foreach(readdir(CDIR)){
            if (/child/){
                $count++;
            }
        }
        closedir(CDIR);
        if ($count > $processor){
	    $wait_time = 0.5;
	}else{
	    $wait_time = 0.1;
	}	    
        if ($processor > $count){
            return 1;
        }
    }
}

sub join{
    while(1){
        my $count = 0;
        sleep 1;
        opendir(CDIR, "$wd/$a/child");
        foreach(readdir(CDIR)){
            if (/child/){
                $count++;
            }
        }
        closedir(CDIR);
        if ($count == 0){
            return 1;
        }
    }
}

sub log{
    my $message = shift;
    my $now = `date "+%Y-%m-%d %H:%M:%S"`;
    chomp($now);
    $message =~ s/\&$//;
    $message = "$now : $message\n";
    $| = 1;
    if ($tsd_size == 20){
	open(LOG, ">> $wd/$a/junction_method.log");
    }else{
	open(LOG, ">> $wd/$a/tsd_method.log.$tsd_size");
    }
    flock(LOG,2);
    print $message;
    print LOG $message;
    close(LOG);
}

sub bynumber{
    $a <=> $b;
}
