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
      perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,option=clear,max_process=8,sort_tmp=/mnt/ssd/tmp

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
 sort_tmp=directory_for_sort is the option for temporary directory for sort command.
 If sort_tmp is specified to fast disk, e.g. SSD, sorting will be accelerated.

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
$sort_tmp = "$wd/$a/tmp" if $sort_tmp eq "";
$th = 0.2 if $th eq "";
$tsd_size = 20 if $tsd_size eq "";
if ($tsd_size == 20){
    $method = "Junction";
}else{
    $method = "TSD";
}

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
    &log("Method : $method");
    &commonMethod;
    if ($tsd_size == 20){
	&junctionMethod;
    }else{
	&tsdMethod;
    }
    &join;
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
    &count($a) if &fragmentLength($a) != 20 + $tsd_size;
    &count($b) if &fragmentLength($b) != 20 + $tsd_size;
    &join;
    &merge($a) if &fragmentLength($a) != 20 + $tsd_size;
    &merge($b) if &fragmentLength($b) != 20 + $tsd_size;
    &join;
}

sub junctionMethod{
    &junctionSpecific;
    &junctionFirstMap($a);
    &junctionFirstMap($b);
    &join;
    &junctionSecondMap($a);
    &junctionSecondMap($b);
    &join;
    &junctionSort($a);
    &junctionSort($b);
    &join;
    &junctionSelectCandidate($a);
    &junctionSelectCandidate($b);
    &join;
    &verify;
    &map;
}

sub tsdMethod{
    &tsdHeadTail;
    &verify;
    if ($ref ne ""){
	&map;
    }
}

sub junctionSelectCandidate{
    $target = shift;
    opendir(REF, "$wd/$ref");
    foreach $file (sort readdir(REF)){
	if ($file =~ /^chr/){
	    &monitorWait;
	    &log("select candidate : $target : $file");
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSelectCandidateFunc,chr=$file,sort_tmp=$sort_tmp &";
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
    &log("junctionSort : Sorting filtered data : $target");
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
    system("rm $wd/$target/tmp/map.*");
    $org_processor = $processor;
    $processor = 2 if -s "$wd/$target/tmp/$chr[0]" > 10000000000;
    foreach $file (@chr){
	&monitorWait;
	&log("junctionSort: Sorting $file");
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSortFunc,chr=$file,sort_tmp=$sort_tmp &";
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : junctionSort : $cmd") if $rc;
    }
    $processor = $org_processor;
}

sub junctionSortFunc{
    $cmd = "sort -k 1 -n -S 1M -T $sort_tmp $wd/$target/tmp/$chr > $wd/$target/tmp/sorted.$chr && rm $wd/$target/tmp/$chr";
    $rc = system($cmd);
    $rc = $rc >> 8;
    &log("ERROR : junctionSortFunc : $cmd") if $rc;
}

sub junctionSecondMap{
    my $target = shift;
    &log("junctionSecondMap : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSecondMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
		&log("$cmd");
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSecondMap : $cmd") if $rc;
            }
        }
    }
    system("rm $wd/$target/tmp/specific.* $wd/$target/tmp/target.*");
}

sub junctionSecondMapFunc{
    system("rm $wd/$target/tmp/map.$tag") if -e "$wd/$target/tmp/map.$tag";
    system("rm $wd/$target/tmp/tmp.$tag") if -e "$wd/$target/tmp/tmp.$tag";
    $cmd = "cat $wd/$target/tmp/first.$tag.* | sort -S 1M -T $sort_tmp > $wd/$target/tmp/second.$tag";
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

sub junctionFirstMap{
    my $target = shift;
    &log("junctionFirstMap : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
		&log($cmd);
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionFirstMap : $cmd") if $rc;
            }
        }
    }
}

sub junctionFirstMapFunc{
    open(IN, "$wd/$target/tmp/specific.$tag");
    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$target/tmp/target.$tag");
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

sub junctionSpecific{
    my $target = shift;
    &log("junctionSpecific : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
                &monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tsd_size=$tsd_size,tag=$tag,sort_tmp=$sort_tmp &";
                &log($cmd);
                $rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSpecific : $cmd") if $rc;
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
	if ($row[1] > 4){
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
    system("rm $wd/$a/tmp/specific.$tag.tmp $wd/$b/tmp/specific.$tag.tmp");
}

sub map{
    system("rm $wd/$a/tmp/map.*") if -e "$wd/$a/tmp/map.AAA";
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		open($tag, "|sort -S 1M -T $sort_tmp > $wd/$a/tmp/mapquery.$tag")
	    }
	}
    }
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/junction_method.verify.$a.$b");
    }
    while(<IN>){
	@row = split;
	$tag = substr($row[0], 0, 3);
	print $tag "$row[0]\n";
	$te{$row[0]} = "head";
	$tag = substr($row[1], 0, 3);
	print $tag "$row[1]\n";
	$te{$row[1]} = "tail";
	$tag = substr($row[2], 0, 3);
	print $tag "$row[2]\n";
	$tag = substr($row[4], 0, 3);
	print $tag "$row[4]\n";

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
		if (-s "$wd/$a/tmp/mapquery.$tag" > 0){
		    &monitorWait;
		    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mapQuery,tag=$ta,sort_tmp=$sort_tmpg &";
		    &log($cmd);
		    $rc = system($cmd);
		    $rc = $rc >> 8;
		    &log("ERROR : verify : $cmd") if $rc;
		}
	    }
	}
    }
    &join;

    open(IN, "cat $wd/$a/tmp/map.* |");
    while(<IN>){
	chomp;
	@row = split;
	$s->{map}{$row[0]}{$row[1]}{$row[2]} = $row[3];
    }
    close(IN);
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/junction_method.verify.$a.$b");
    }
    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$a/tmp/mapped.tmp");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
	for ($i = 5; $i <= 10; $i++){
	    $row[$i] += 0;
	}
	$tmp = join("\t", @row);
	foreach $uchr (sort keys %{$s->{map}{$row[2]}}){
	    foreach $upos (sort keys %{$s->{map}{$row[2]}{$uchr}}){
		foreach $dchr (sort keys %{$s->{map}{$row[4]}}){
		    foreach $dpos (sort keys %{$s->{map}{$row[4]}{$dchr}}){
			if ($uchr eq $dchr and abs($upos - $dpos) < 37 and $s->{map}{$row[2]}{$uchr}{$upos} eq $s->{map}{$row[4]}{$dchr}{$dpos}){
			    $direction = $s->{map}{$row[2]}{$uchr}{$upos};
			    if ($direction eq "f"){
				$upos += 19;
			    }else{
				$upos -= 19;
			    }
			    if ($row[5] > 3 and $row[6] > 3){
				if ($row[7] <= 2){
				    $genotype = "M";
				}else{
				    $genotype = "H";
				}
				$uchr = "000" . $uchr if $uchr !~ /^0/;
				$uchr = substr($uchr, length($uchr) - 3, 3);
				$upos = "00000000000" . $upos if $upos !~ /^0/;
				$upos = substr($upos, length($upos) - 11, 11);
				print OUT "$a\t$uchr\t$upos\t$dpos\t$direction\t$tmp\t$genotype\n";
			    }
			    if ($row[8] > 3 and $row[9] > 3){
				if ($row[10] <= 2){
				    $genotype = "M";
				}else{
				    $genotype = "H";
				}
				$uchr = "000" . $uchr if $uchr !~ /^0/;
				$uchr = substr($uchr, length($uchr) - 3, 3);
				$upos = "00000000000" . $upos if $upos !~ /^0/;
				$upos = substr($upos, length($upos) - 11, 11);
				print OUT "$b\t$uchr\t$upos\t$dpos\t$direction\t$tmp\t$genotype\n";
			    }
			}
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);
    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.genotype.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.genotype.$a.$b");
    }
    open(IN, "$wd/$a/tmp/mapped.tmp");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
	$row[1] =~ s/^0+//;
	$row[2] =~ s/^0+//;
	$result = join("\t", @row) . "\n";
	print $result;
	print OUT $result;
    }
    close(IN);
    close(OUT);


    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$a/tmp/temap");
    foreach $seq (sort keys %{$s->{map}}){
	if ($te{$seq}){
	    foreach $chr (sort keys %{$s->{map}{$seq}}){
		foreach $pos (sort keys %{$s->{map}{$seq}{$chr}}){
		    $direction = $s->{map}{$seq}{$chr}{$pos};
		    $ochr = "000" . $chr;
		    $ochr = substr($ochr, length($ochr) - 3, 3);
		    $opos = "00000000000" . $pos;
		    $opos = substr($opos, length($opos) - 11, 11);
		    print OUT "$ochr\t$opos\t$te{$seq}\t$seq\t$direction\n";
		}
	    }
	}
    }
    close(OUT);

    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.tepos.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.tepos.$a.$b");
    }
    open(IN, "$wd/$a/tmp/temap");
    while(<IN>){
	@row = split;
	$row[0] =~ s/^0+//g;
	$row[1] =~ s/^0+//g;
	$output = join("\t", @row) . "\n";
	print $output;
	print OUT$output;
    }
    close(IN);
    close(OUT);
}

sub mapQuery{
    $cmd = "zcat $wd/$ref/ref20.$tag.gz | join - $wd/$a/tmp/mapquery.$tag > $wd/$a/tmp/map.$tag";
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify :$cmd") if $rc;
}

sub junctionSecondMap{
    my $target = shift;
    &log("junctionSecondMap : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSecondMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
		&log("$cmd");
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSecondMap : $cmd") if $rc;
            }
        }
    }
    system("rm $wd/$target/tmp/specific.* $wd/$target/tmp/target.*");
}

sub junctionSecondMapFunc{
    system("rm $wd/$target/tmp/map.$tag") if -e "$wd/$target/tmp/map.$tag";
    system("rm $wd/$target/tmp/tmp.$tag") if -e "$wd/$target/tmp/tmp.$tag";
    $cmd = "cat $wd/$target/tmp/first.$tag.* | sort -S 1M -T $sort_tmp > $wd/$target/tmp/second.$tag";
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

sub junctionFirstMap{
    my $target = shift;
    &log("junctionFirstMap : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
		&log($cmd);
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionFirstMap : $cmd") if $rc;
            }
        }
    }
}

sub junctionFirstMapFunc{
    open(IN, "$wd/$target/tmp/specific.$tag");
    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$target/tmp/target.$tag");
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

sub junctionSpecific{
    my $target = shift;
    &log("junctionSpecific : $target");
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
                &monitorWait;
		$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tsd_size=$tsd_size,tag=$tag,sort_tmp=$sort_tmp &";
                &log($cmd);
                $rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : junctionSpecific : $cmd") if $rc;
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
    system("rm $wd/$a/tmp/specific.$tag.tmp $wd/$b/tmp/specific.$tag.tmp");
}

sub map{
    system("rm $wd/$a/tmp/map.*") if -e "$wd/$a/tmp/map.AAA";
    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		open($tag, "|sort -S 1M -T $sort_tmp > $wd/$a/tmp/mapquery.$tag")
	    }
	}
    }
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/junction_method.verify.$a.$b");
    }
    while(<IN>){
	@row = split;
	$tag = substr($row[0], 0, 3);
	print $tag "$row[0]\n";
	$te{$row[0]} = "head";
	$tag = substr($row[1], 0, 3);
	print $tag "$row[1]\n";
	$te{$row[1]} = "tail";
	$tag = substr($row[2], 0, 3);
	print $tag "$row[2]\n";
	$tag = substr($row[4], 0, 3);
	print $tag "$row[4]\n";

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
		if (-s "$wd/$a/tmp/mapquery.$tag" > 0){
		    &monitorWait;
		    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mapQuery,tag=$tag,sort_tmp=$sort_tmp &";
		    &log($cmd);
		    $rc = system($cmd);
		    $rc = $rc >> 8;
		    &log("ERROR : verify : $cmd") if $rc;
		}
	    }
	}
    }
    &join;

    open(IN, "cat $wd/$a/tmp/map.* |");
    while(<IN>){
	chomp;
	@row = split;
	$s->{map}{$row[0]}{$row[1]}{$row[2]} = $row[3];
    }
    close(IN);
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/junction_method.verify.$a.$b");
    }
    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$a/tmp/mapped.tmp");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
	for ($i = 5; $i <= 10; $i++){
	    $row[$i] += 0;
	}
	$tmp = join("\t", @row);
	foreach $uchr (sort keys %{$s->{map}{$row[2]}}){
	    foreach $upos (sort keys %{$s->{map}{$row[2]}{$uchr}}){
		foreach $dchr (sort keys %{$s->{map}{$row[4]}}){
		    foreach $dpos (sort keys %{$s->{map}{$row[4]}{$dchr}}){
			if ($uchr eq $dchr and abs($upos - $dpos) < 37 and $s->{map}{$row[2]}{$uchr}{$upos} eq $s->{map}{$row[4]}{$dchr}{$dpos}){
			    $direction = $s->{map}{$row[2]}{$uchr}{$upos};
			    if ($direction eq "f"){
				$upos += 20;
				$dpos -= 1;
			    }else{
				$upos -= 20;
				$dpos += 1;
			    }
			    if ($row[5] > 3 and $row[6] > 3){
				if ($row[7] <= 2){
				    $genotype = "M";
				}else{
				    $genotype = "H";
				}
				$uchr = "000" . $uchr if $uchr !~ /^0/;
				$uchr = substr($uchr, length($uchr) - 3, 3);
				$upos = "00000000000" . $upos if $upos !~ /^0/;
				$upos = substr($upos, length($upos) - 11, 11);
				print OUT "$a\t$uchr\t$upos\t$dpos\t$direction\t$tmp\t$genotype\n";
			    }
			    if ($row[8] > 3 and $row[9] > 3){
				if ($row[10] <= 2){
				    $genotype = "M";
				}else{
				    $genotype = "H";
				}
				$uchr = "000" . $uchr if $uchr !~ /^0/;
				$uchr = substr($uchr, length($uchr) - 3, 3);
				$upos = "00000000000" . $upos if $upos !~ /^0/;
				$upos = substr($upos, length($upos) - 11, 11);
				print OUT "$b\t$uchr\t$upos\t$dpos\t$direction\t$tmp\t$genotype\n";
			    }
			}
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);
    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.genotype.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.genotype.$a.$b");
    }
    open(IN, "$wd/$a/tmp/mapped.tmp");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
	$row[1] =~ s/^0+//;
	$row[2] =~ s/^0+//;
	$result = join("\t", @row) . "\n";
	print $result;
	print OUT $result;
    }
    close(IN);
    close(OUT);


    open(OUT, "| sort -S 1M -T $sort_tmp > $wd/$a/tmp/temap");
    foreach $seq (sort keys %{$s->{map}}){
	if ($te{$seq}){
	    foreach $chr (sort keys %{$s->{map}{$seq}}){
		foreach $pos (sort keys %{$s->{map}{$seq}{$chr}}){
		    $direction = $s->{map}{$seq}{$chr}{$pos};
		    $ochr = "000" . $chr;
		    $ochr = substr($ochr, length($ochr) - 3, 3);
		    $opos = "00000000000" . $pos;
		    $opos = substr($opos, length($opos) - 11, 11);
		    print OUT "$ochr\t$opos\t$te{$seq}\t$seq\t$direction\n";
		}
	    }
	}
    }
    close(OUT);

    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.tepos.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.tepos.$a.$b");
    }
    open(IN, "$wd/$a/tmp/temap");
    while(<IN>){
	@row = split;
	$row[0] =~ s/^0+//g;
	$row[1] =~ s/^0+//g;
	$output = join("\t", @row) . "\n";
	print $output;
	print OUT$output;
    }
    close(IN);
    close(OUT);
}

sub mapQuery{
    $cmd = "zcat $wd/$ref/ref20.$tag.gz | join - $wd/$a/tmp/mapquery.$tag > $wd/$a/tmp/map.$tag";
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify :$cmd") if $rc;
}

sub verify{
    opendir(DIR, "$a/split");
    foreach (sort readdir(DIR)){
	if (/^split/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$a,sub=verifyFunc,tsd_size=$tsd_size,file=$_,sort_tmp=$sort_tmp &";
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
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$b,sub=verifyFunc,tsd_size=$tsd_size,file=$_,sort_tmp=$sort_tmp &";
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
	if (/verify.split/){
	    open(IN, "$a/tmp/$_");
	    while(<IN>){
		chomp;
		@row = split;
		$s->{$row[0]}{$row[1]}{$row[2]}{$row[3]} = 1;
	    }
	    close(IN);
	}
    }
    opendir(DIR, "$b/tmp");
    foreach (sort readdir(DIR)){
	if (/verify.split/){
	    open(IN, "$b/tmp/$_");
	    while(<IN>){
		chomp;
		@row = split;
		$s->{$row[0]}{$row[1]}{$row[2]}{$row[3]} = 1;
	    }
	close(IN);
	}
    }

    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/query");
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    }else{
	open(IN, "$wd/$a/tmp/te.candidate");
    }
    while(<IN>){
	@row = split;
	foreach $upstream (sort keys %{$s->{$row[0]}{H}}){
	    foreach $htsd (sort keys %{$s->{$row[0]}{H}{$upstream}}){
		foreach $downstream (sort keys %{$s->{$row[1]}{T}}){
		    foreach $ttsd (sort keys %{$s->{$row[1]}{T}{$downstream}}){
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
	}
    }
    close(IN);
    close(OUT);

    foreach $nuc (@nuc){
        $tag[0] = $nuc;
        foreach $nuc (@nuc){
            $tag[1] = $nuc;
            foreach $nuc (@nuc){
                $tag[2] = $nuc;
                $tag = join('', @tag);
		open($tag, "> $wd/$a/tmp/query.$tag")
	    }
	}
    }
    open(IN, "$wd/$a/tmp/query");
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
		if (-s "$wd/$a/tmp/query.$tag" > 0){
		    &monitorWait;
		    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=countQuery,tag=$tag,sort_tmp=$sort_tmp &";
		    &log($cmd);
		    $rc = system($cmd);
		    $rc = $rc >> 8;
		    &log("ERROR : verify : $cmd") if $rc;
		}
	    }
	}
    }
    &join;
    open(IN, "cat $wd/$a/tmp/verify.count.* |");
    while(<IN>){
	chomp;
	@row = split;
	$a{$row[0]} = $row[1];
	
    }
    open(IN, "cat $wd/$b/tmp/verify.count.* |");
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
	    foreach $htsd (sort keys %{$s->{$row[0]}{H}{$upstream}}){
		foreach $downstream (sort keys %{$s->{$row[1]}{T}}){
		    foreach $ttsd (sort keys %{$s->{$row[1]}{T}{$downstream}}){
			if ($htsd eq $ttsd){
			    $tsdsize = length($htsd);
			    $wildtype = $upstream . substr($downstream, $tsdsize, 20);
			    $head = $upstream . $row[0];
			    $downfl = substr($downstream, 0, 20);
			    $tail = $row[1] . $downfl;
			    if (($a{$head} > 2 and $a{$tail} > 2 and $b{$wildtype} > 2) or ($b{$head} > 2 and $b{$tail} > 2 and $a{$wildtype} > 2)){
				$result = "$row[0]\t$row[1]\t$upstream\t$htsd\t$downfl\t$a{$head}\t$a{$tail}\t$a{$wildtype}\t$b{$head}\t$b{$tail}\t$b{$wildtype}\n";
				print $result;
				print OUT $result;
				if ($a{$head} > 2){
				    $s->{tecount}{$a}{"$row[0]\t$row[1]"} ++;
				    $s->{tetsd}{$a}{"$row[0]\t$row[1]"}{$htsd} ++;
				}
				if ($b{$head} > 2){
				    $s->{tecount}{$b}{"$row[0]\t$row[1]"} ++;
				    $s->{tetsd}{$b}{"$row[0]\t$row[1]"}{$htsd} ++;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);
    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.summary.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.summary.$a.$b");
    }
    foreach $line (sort keys %{$s->{tecount}}){
	foreach $te (sort keys %{$s->{tecount}{$line}}){
	    $total = $s->{tecount}{$line}{$te};
	    $count = 0;
	    foreach (sort keys %{$s->{tetsd}{$line}{$te}}){
		$count ++;
	    }
	    print OUT "$line\t$te\t$count / $total\n";
	}
    }
    close(OUT);
}

sub countQuery{
    $cmd = "zcat $wd/$a/count.20/$tag.gz | join - $wd/$a/tmp/query.$tag > $wd/$a/tmp/verify.count.$tag";
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify :$cmd") if $rc;
    
    $cmd = "zcat $wd/$b/count.20/$tag.gz | join - $wd/$a/tmp/query.$tag > $wd/$b/tmp/verify.count.$tag";
    $rc = system("$cmd");
    $rc = $rc >> 8;
    &log("ERROR : verify : $cmd") if $rc;
}

sub verifyFunc{
    my (@row, $tsdsize, %seq, %tsdsize, $length, $seq, $tsd, $upstream, $downstream);
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
	$s->{tsdsize}{$row[0]}{$tsdsize} = 1;
	$s->{tsdsize}{$row[1]}{$tsdsize} = 1;
	$seq{complement($row[0])} = "T";
	$seq{complement($row[1])} = "H";
	$s->{tsdsize}{complement($row[0])}{$tsdsize} = 1;
	$s->{tsdsize}{complement($row[1])}{$tsdsize} = 1;
	
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
		foreach $tsdsize (sort keys %{$s->{tsdsize}{$seq}}){		
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
    }
    close(IN);
    close(OUT);
}

sub tsdHeadTail{
    &log("detect head and tail sequences");
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		$cmd = "perl $0 a=$a,b=$b,sub=tsdHeadTailFunc,tag=$tag,tsd_size=$tsd_size,sort_tmp=$sort_tmp &";
		&log("tsdHeadTail : $tag : $tsd_size");
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : tsdHeadTail : $cmd") if $rc;
	    }
	}
    }
    &join;
    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/head.select");
    open(IN, "cat $wd/$a/tmp/*.head | sort -S 1M -T $sort_tmp |");
    while(<IN>){
	chomp;
	@row  = split;
	if ($row[0] eq $prev[0]){
	    print OUT "$prev[1]\t$prev[0]\t$prev[2]\n
$row[1]\t$row[0]\t$row[2]\n";
	}
 	@prev = @row;      
    }
    close(IN);
    close(OUT);
    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/tail.select");
    open(IN, "cat $wd/$a/tmp/*.tail | sort -S 1M -T $sort_tmp |");
    while(<IN>){
	chomp;
	@row  = split;
	if ($row[0] eq $prev[0]){
	    print OUT "$prev[1]\t$prev[0]\t$prev[2]\n
$row[1]\t$row[0]\t$row[2]\n";
	}
	@prev = @row;
    }
    close(IN);
    close(OUT);
    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$b/tmp/head.select");
    open(IN, "cat $wd/$b/tmp/*.head | sort -S 1M -T $sort_tmp |");
    while(<IN>){
	chomp;
	@row  = split;
	if ($row[0] eq $prev[0]){
	    print OUT "$prev[1]\t$prev[0]\t$prev[2]\n
$row[1]\t$row[0]\t$row[2]\n";
	}
	@prev = @row;
    }
    close(IN);
    close(OUT);
    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$b/tmp/tail.select");
    open(IN, "cat $wd/$b/tmp/*.tail | sort -S 1M -T $sort_tmp |");
    while(<IN>){
	chomp;
	@row  = split;
	if ($row[0] eq $prev[0]){
	    print OUT "$prev[1]\t$prev[0]\t$prev[2]\n
$row[1]\t$row[0]\t$row[2]\n";
	}
	@prev = @row;
    }
    close(IN);
    close(OUT);

    &log("join head.select tail.select : $a");
    open(OUT, "|sort -S 1M -T $sort_tmp > $wd/$a/tmp/pair");
    open(IN, "join $wd/$a/tmp/head.select $wd/$a/tmp/tail.select|");
    while(<IN>){
	chomp;
	@row = split;
	print OUT "$row[1]_$row[3]\t$row[0]\t$row[2]\t$row[4]\n";
    }
    close(IN);
    close(OUT);
    &log("join head.select tail.select : $b");
    open(OUT, "|sort -S 1M -T $sort_tmp > $wd/$b/tmp/pair");
    open(IN, "join $wd/$b/tmp/head.select $wd/$b/tmp/tail.select|");
    while(<IN>){
	chomp;
	@row = split;
	print OUT "$row[1]_$row[3]\t$row[0]\t$row[2]\t$row[4]\n";
    }
    close(IN);
    close(OUT);
    
    &log("tsdHeadTail : tsd list : $a");
    open(OUT, "> $wd/$a/tmp/tsd");
    open(IN, "$wd/$a/tmp/pair");
    while(<IN>){
	chomp;
	@row = split;
	if ($prev eq $row[0]){
	    $tsd .= $ptsd . ":";
	}elsif ($tsd ne ""){
	    print OUT "$prev\t$tsd" . "$ptsd\n";
	    $tsd = "";
	}
	$prev = $row[0];
	$ptsd = $row[1];
    }
    print OUT "$prev\t$tsd" . "ptsd" if $tsd ne "";
    close(IN);
    close(OUT);
    
    &log("tsdHeadTail : tsd list : $b");
    open(OUT, "> $wd/$b/tmp/tsd");
    open(IN, "$wd/$b/tmp/pair");
    while(<IN>){
	chomp;
	@row = split;
	if ($prev eq $row[0]){
	    $tsd .= $ptsd . ":";
	}elsif ($tsd ne ""){
	    print OUT "$prev\t$tsd" . "$ptsd\n";
	    $tsd = "";
	}
	$prev = $row[0];
	$ptsd = $row[1];
    }
    print OUT "$prev\t$tsd" . "ptsd" if $tsd ne "";
    close(IN);
    close(OUT);
    open(OUT, "> $wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    open(IN, "join $wd/$a/tmp/tsd $wd/$b/tmp/tsd|");
    while(<IN>){
	chomp;
	@row = split;
	($head, $tail) = split('_', $row[0]);
	$row[1] =~ s/:/ /g;
	$row[2] =~ s/:/ /g;
	print OUT "$head\t$tail\t$row[1]\t$row[2]\n";
    }
}

sub tsdHeadTailFunc{
    open(AH, "> $wd/$a/tmp/$tag.head");
    open(AT, "> $wd/$a/tmp/$tag.tail");
    open(IN, "bash -c 'join -a 1 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz) '|");
    while(<IN>){
	chomp;
	@row = split;
	next if ($row[2] ne "" or $row[1] <= 4);
	next if ($row[0] =~ /A{8,}/);
	next if ($row[0] =~ /C{8,}/);
	next if ($row[0] =~ /G{8,}/);
	next if ($row[0] =~ /T{8,}/);
	$head = substr($row[0], $tsd_size);
	$tsd  = substr($row[0], 0, $tsd_size);
	print AH "$head\t$tsd\t$row[1]\n";
	$tail = substr($row[0], 0, 20);
	$tsd = substr($row[0], 20, $tsd_size);
	print AT "$tail\t$tsd\t$row[1]\n";
    }
    close(IN);
    close(AH);
    close(AT);
    open(BH, "> $wd/$b/tmp/$tag.head");
    open(BT, "> $wd/$b/tmp/$tag.tail");
    open(IN, "bash -c 'join -a 2 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz) '|");
    while(<IN>){
	chomp;
	@row = split;
	next if ($row[2] ne "" or $row[1] <= 4);
	next if ($row[0] =~ /A{8,}/);
	next if ($row[0] =~ /C{8,}/);
	next if ($row[0] =~ /G{8,}/);
	next if ($row[0] =~ /T{8,}/);
	$head = substr($row[0], $tsd_size);
	$tsd  = substr($row[0], 0, $tsd_size);
	print BH "$head\t$tsd\t$row[1]\n";
	$tail = substr($row[0], 0, 20);
	$tsd = substr($row[0], 20, $tsd_size);
	print BT "$tail\t$tsd\t$row[1]\n";
    }
    close(IN);
    close(BH);
    close(BT);
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

sub merge{
    my $target = shift;
    &log("merge subfiles : $target");
    foreach $nuc (@nuc){
	$taga = $nuc;
	foreach $nuc (@nuc){
	    $tagb = $nuc;
	    foreach $nuc (@nuc){
		$tagc = $nuc;
		$tag = $taga . $tagb . $tagc;
		&monitorWait;
		$cmd = "perl $0 target=$target,sub=mergeFunc,tsd_size=$tsd_size,tag=$tag,a=$a,sort_tmp=$sort_tmp &";
		&log($cmd);
		$rc = system($cmd);
		$rc = $rc >> 8;
		&log("ERROR : merge :$cmd") if $rc;
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

sub count{
    my $target = shift;
    my ($last, $i, $cmd);
    &log("count subfiles : $target");
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
	$cmd = "perl $0 target=$target,sub=countFunc,number=$i,a=$a,tsd_size=$tsd_size,sort_tmp=$sort_tmp &";
	&log($cmd);
	$rc = system($cmd);
	$rc = $rc >> 8;
	&log("ERROR : count : $cmd") if $rc;
    }
}

sub countFunc{
    my ($taga, $tagb, $tagc, $tag, $count, $nuc);
    open(IN, "$wd/$target/split/split.$number");
    open(SEQ, "|sort -S 1M -T $sort_tmp |uniq -c |awk '{print \$2 \"\t\" \$1}' > $wd/$target/count.$tsd_size/count.$number");
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

 Replace line to 'NOP' for the chromosome entry which is not required.

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
    $cmd = "sort -S 1M -T $sort_tmp $wd/$ref/tmp/ref20.$tag.* |gzip -c > $wd/$ref/ref20.$tag.gz";
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
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mk20mer,chr=$i,sort_tmp=$sort_tmp &";
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
		    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=sort20mer,tag=$tag,sort_tmp=$sort_tmp &";
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
