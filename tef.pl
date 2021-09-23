#!/usr/bin/perl
#
# tef.pl - transposable element finder
#          A program for transposable elements detection from NGS reads.
#
# Copyright (C) 2020 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/tef
#
# Author: MIYAO Akio <miyao@affrc.go.jp>
#


if ($ARGV[0] eq ""){
    print "
 tef.pl - transposable element finder
          A program for detection of active transposable elements from NGS reads.

 e.g. perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0
      perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,method=tsd
      perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5,th=0.7

 tef.pl requires a set of NGS reads from different conditions.
 For example, ttm2 and ttm5 are from regenerated rice individuals from callus.
 Save fastq files into ttm2/read and ttm5/read directory.

 Options shold be specified with 'name=value' separated by a comma without space.
 Default method is the junction method.
 For the tsd method, option 'method=tsd' is required.

 If ref is specified, fasta file of reference genome sequence compressed
 by gzip should be saved into the directory specified with ref. 
 
 If ref is not specified, a or b specific transpositions will be detected.
 This function works only by tsd method.

 If tsd_size is not specified, detection will be progressed with
 default tsd_size 5.

 If th is not specified, default value 0.7 will be used.
 If no TE is detected, try again with lowere value e.g. th=0.3  

 Author: MIYAO Akio

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
$tsd_size = 5 if $tsd_size eq "";
$th = 0.7 if $th eq "";

if ($sub eq ""){
    system("rm -rf $wd/$a/tmp") if -e "$wd/$a/tmp";
    system("rm -rf $wd/$b/tmp") if -e "$wd/$b/tmp";
    system("rm -rf $wd/$a/child") if -e "$wd/$a/child";
    system("mkdir $wd/$a/tmp");
    system("mkdir $wd/$b/tmp");
    system("mkdir $wd/$a/child");
    system("rm $wd/$a/log") if -e "$wd/$a/log";
    open(LOG, "> $wd/$a/log");
    foreach (split(',', $ARGV[0])){
	my ($name, $val) = split('=', $_);
	print LOG "$name : $val\n";
    }
    close(LOG);
    &report("job start");
    &commonMethod;
    if ($method eq "" or $method =~ /jun/){
	print LOG "method : junction\n" if $method eq "";;
	$method = "junctionMethod";
	&junctionMethod;
    }else{
	print LOG "method : tsd\n";
	$method = "tsdMethod";
	&tsdMethod;
    }
    &waitChild;
    system("rm -rf $wd/$a/tmp");
    system("rm -rf $wd/$b/tmp");
    system("rm -rf $wd/$a/child");
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
    report("Job completed: Analysis has been done. $a $b $ref $method
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
    report("split to subfiles : $a");
    system("perl $0 target=$a,sub=split,a=$a,wd=$wd &");
    &canFork;
    report("split to subfiles : $b");
    system("perl $0 target=$b,sub=split,a=$a,wd=$wd &");
    &waitChild;
}

sub junctionMethod{
    report("count subfiles : $a");
    &count($a);
    report("count subfiles : $b");
    &count($b);
    report("merge subfiles : $a");
    &waitChild;
    &merge($a);
    report("merge subfiles : $b");
    &merge($b);
    &waitChild;
    report("specific");
    &junctionSpecific;
    report("firstmap : $a");
    &junctionFirstMap($a);
    report("firstmap : $b");
    &junctionFirstMap($b);
    &waitChild;
    report("secondmap : $a");
    &junctionSecondMap($a);
    report("secondmap : $b");
    &junctionSecondMap($b);
    &waitChild;
    report("sort : $a");
    &junctionSort($a);
    report("sort : $b");
    &junctionSort($b);
    &waitChild;
    report("select candidate : $a");
    &junctionSelectCandidate($a);
    report("select candidate : $b");
    &junctionSelectCandidate($b);
    &waitChild;
    &junctionMapSelection($a);
    &junctionMapSelection($b);
    &waitChild;
    &junctionTsdSelection($a);
    &junctionTsdSelection($b);
}

sub tsdMethod{
    report("count subfiles : $a");
    &count($a);
    report("count subfiles : $b");
    &count($b);
    &waitChild;
    report("merge subfiles : $a");
    &merge($a);
    report("merge subfiles : $b");
    &merge($b);
    &waitChild;
    report("comm");
    &comm;
    &waitChild;
    report("join");
    &join;
    &waitChild;
    report("tsd");
    &tsd;
    foreach $taga (@nuc){
	foreach $tagb (@nuc){
	    foreach $tagc (@nuc){
		$tag = $taga . $tagb . $tagc;
	       &canFork;
		report("pair $tag");
		system("perl $0 a=$a,b=$b,sub=pair,tag=$tag,tsd_size=$tsd_size,wd=$wd &");
	    }
	}
    }
    &waitChild;
    report("sortPair");
    &sortPair;
    foreach $taga (@nuc){
	foreach $tagb (@nuc){
	    foreach $tagc (@nuc){
		$tag = $taga . $tagb . $tagc;
		&canFork;
		report("selectPair $tag");
		system("perl $0 a=$a,b=$b,sub=selectPair,tag=$tag,tsd_size=$tsd_size,th=$th,wd=$wd &");
	    }
	}
    }
    &waitChild;
    system("cat $wd/$a/tmp/pair.$a.$b.$tsd_size.* > $wd/$a/tsd_method.pair.$a.$b.$tsd_size && rm $wd/$a/tmp/pair.$a.$b.$tsd_size.*");
    if($ref eq ""){
	&verify;
	report("Job end. $a $b");
	exit;
    }
    &mkQuery;
    &map;
}

sub junctionTsdSelection{
    my $target = shift;
    my $file;
    opendir(DIR, "$wd/$target");
    foreach $file (sort readdir(DIR)){
	if($file =~ /junction_method.map.$target/){
	    &canFork;
	    chomp;
	    report("tsd selection : $target : $file");
	    system("perl $0 a=$a,b=$b,sub=junctionTsdSelectionFunc,target=$target,file=$file,wd=$wd &");
	    
	}
    }
    closedir(DIR);
}

sub junctionTsdSelectionFunc{
    my $s = {};
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
    foreach $chr (sort bynumber keys %{$s->{head}}){
	foreach $pos (sort bynumber keys %{$s->{head}{$chr}}){
	    $hf = $s->{head}{$chr}{$pos};
	    for ($i = 1; $i <= 38; $i++){
		$tf = $s->{tail}{$chr}{$pos - 20 + $i};
		if ($tf ne ""){
		    if ($i < 20){
			$size = 20 - $i + 1;
			
		    }else{
			$size = $i - 20 + 1;
		    }
		    $htsd = substr($hf, 20 - $size, $size);
		    $ttsd = substr($tf, 0, $size);
		    if ($htsd eq $ttsd){
			$s->{tsd}{$chr}{$pos} = $htsd;
			$s->{tsd}{$chr}{$pos - 20 + $i} = $ttsd;
		    }
		}
	    }
	}
    }
    open(IN, "$wd/$target/$file");
    ($head, $tail) = (split('\.', $file))[3,4];
    open(OUT, ">$wd/$target/junction_method.$target.$head.$tail");
    while(<IN>){
	@row = split;
	if ($s->{tsd}{$row[1]}{$row[2]}){
	    $output =  "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$s->{tsd}{$row[1]}{$row[2]}\t$row[5]\n";
	    print $output;
	    print OUT $output;
	}
    }
    close(IN);
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
	    open(IN, "$wd/$target/tmp/te.candidate");
	    while(<IN>){
		&canFork;
		chomp;
		($head, $tail) = (split)[0, 1];
		report("map selection : $target : $chr $head $tail");
		system("perl $0 a=$a,b=$b,sub=junctionMapSelectionFunc,target=$target,chr=$chr,head=$head,tail=$tail,wd=$wd &");
	    }
	    close(IN);
	}
    }
    closedir(REF);
    &waitChild;
    open(IN, "$wd/$target/tmp/te.candidate");
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
#	system("rm $wd/$target/pos.$head.$tail");
    }
    close(IN);
#    system("rm $wd/$target/tmp/pos.head.$head.* $wd/$target/tmp/pos.tail.$tail.*");
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
	    &canFork;
	    &report("select candidate : $target : $file");
	    system("perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSelectCandidateFunc,chr=$file,wd=$wd &");
	}
    }
    &waitChild;
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
	$tsd_size = $row[0] - $prev[0] - 10;
	if ($tsd_size <= 10 and $tsd_size > 3){
	    if ($row[1] eq "h" and $prev[1] eq "t"){
		$htsd = substr($row[2], 15, $tsd_size);
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
		$htsd = substr($prev[2], 15, $tsd_size);
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
    &report("Sorting filtered data : $target");
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
	&canFork;
	&report("Sorting $file");
	system("perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSortFunc,chr=$file,wd=$wd &");
    }
    $processor = $org_processor;
}

sub junctionSortFunc{
    system("sort -k 1 -n -S 1M -T $wd/$target/tmp $wd/$target/tmp/$chr > $wd/$target/tmp/sorted.$chr && rm $wd/$target/tmp/$chr");
}

sub junctionSecondMapFunc{
    system("rm $wd/$target/tmp/map.$tag") if -e "$wd/$target/tmp/map.$tag";
    system("rm $wd/$target/tmp/tmp.$tag") if -e "$wd/$target/tmp/tmp.$tag";
    system("cat $wd/$target/tmp/first.$tag.* | sort -S 1M -T $wd/$target/tmp > $wd/$target/tmp/second.$tag");
    open(IN, "zcat $wd/$ref/ref20.$tag.gz |join $wd/$target/tmp/second.$tag - |");
    while(<IN>){
	chomp;
        @row = split;
	if ($prev ne $row[1]){
	    close($tag);
	    if (! $flag and -e "$wd/$target/tmp/tmp.$tag"){
		system("cat $wd/$target/tmp/tmp.$tag >> $wd/$target/tmp/map.$tag");
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
	system("cat $wd/$target/tmp/tmp.$tag >> $wd/$target/tmp/map.$tag");
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
		&canFork;
		&report("Mapping $target : $tag");
		system("perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSecondMapFunc,tag=$tag,wd=$wd &");
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
		&canFork;
		&report("perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag,wd=$wd");
		system("perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag,wd=$wd &");
            }
        }
    }
}

sub junctionSpecificFunc{
    system("bash -c 'join -v 1 <(zcat $wd/$a/tmp/$tag.gz) <(zcat $wd/$b/tmp/$tag.gz) > $wd/$a/tmp/specific.$tag'");
    system("bash -c 'join -v 2 <(zcat $wd/$a/tmp/$tag.gz) <(zcat $wd/$b/tmp/$tag.gz) > $wd/$b/tmp/specific.$tag'");
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
                &canFork;
                &report("perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tag=$tag,wd=$wd");
                system("perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tag=$tag,wd=$wd &");
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
    report("map");
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
                &canFork;
                &report("perl $0 a=$a,b=$b,ref=$ref,sub=mapFunc,tag=$tag,wd=$wd");
                system("perl $0 a=$a,b=$b,ref=$ref,sub=mapFunc,tag=$tag,wd=$wd &");
            }
        }
    }
    &waitChild;
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
		}
	    }else{
		$tsd = substr($prev[6], 20 - $tsd_size, $tsd_size);
		$ttsd = substr($row[6], 0, $tsd_size);
		if ($tsd eq $ttsd){
		    $s->{$row[0]}{$row[2]}{$row[3]}{$row[1]} = "$row[0]\t$row[2]\t$prev[3]\t$row[3]\t$row[1]\t$prev[5]\t$row[5]\t$tsd\t$prev[6]\t$row[6]";
		    $flanking{"$prev[6] $row[6]"} ++;
		}
	    }
	}
	$pline = $_;
	$prev = $row[3];
    }
    close(IN);
    
    open(OUT, "> $wd/$a/tsd_method.result.$a.$b.$tsd_size");
    foreach $line (sort keys %$s){
        foreach $chr (sort bynumber keys %{$s->{$line}}){
            foreach $pos (sort bynumber keys %{$s->{$line}{$chr}}){
		foreach $direction (sort bynumber keys %{$s->{$line}{$chr}{$pos}}){
	   	    @row = split('\t', $s->{$line}{$chr}{$pos}{$direction});
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
    close(OUT);
}

sub searchQuery{
    if ($b eq ""){
	$cmd = "grep $seq $wd/$a/tmp/split.$number > $wd/$a/tmp/selected.$type.$seq.$number";
    }else{
	$cmd = "grep $seq $wd/$b/tmp/split.$number > $wd/$b/tmp/selected.$type.$seq.$number";
    }
    system($cmd);
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

    opendir(DIR, "$wd/$a/tmp");
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
	    &canFork;
	    $cmd = "perl $0 a=$a,sub=searchQuery,number=$i,type=head,seq=$head,wd=$wd &";
	    &report($cmd);
	    system($cmd);
	}
        foreach $tail (sort keys %tail){
            &canFork;
            $cmd = "perl $0 a=$a,sub=searchQuery,number=$i,type=tail,seq=$tail,wd=$wd &";
            &report($cmd);
            system($cmd);
        }
	&waitChild;
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
    opendir(DIR, "$wd/$b/tmp");
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
	    &canFork;
	    $cmd = "perl $0 a=$a,b=$b,sub=searchQuery,number=$i,type=head,seq=$head,wd=$wd &";
	    &report($cmd);
	    system($cmd);
	}
	foreach $tail (sort keys %tail){
	    &canFork;
	    $cmd = "perl $0 a=$a,b=$b,sub=searchQuery,number=$i,type=tail,seq=$tail,wd=$wd &";
	    &report($cmd);
	    system($cmd);
	}
	&waitChild;
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
	report("making config file of in $ref directory.");
	system("zcat *.gz|grep '>' | sed  -e 's/^>//' > config.tmp");
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
		    print "mkref. processing $chr_name\n";
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
    system("sort -S 1M -T $wd/$ref/tmp $wd/$ref/tmp/ref20.$tag.* |gzip -c > $wd/$ref/ref20.$tag.gz");
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
	&report("Making 20mer position file.");
	foreach $i (@chr){
	    next if $i eq "NOP";
	    &canFork;
	    &report("Processing chr$i");
	    system("perl $0 a=$a,b=$b,ref=$ref,sub=mk20mer,chr=$i,wd=$wd &");
	}
	&waitChild;
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
		    &canFork;
		    &report("making ref20.$tag.gz");
		    system("perl $0 a=$a,b=$b,ref=$ref,sub=sort20mer,tag=$tag,wd=$wd &");
		}
	    }
	}
    }    
    &waitChild;
    system("rm -r $wd/$ref/tmp")
}

sub verify{
    report("Verify");
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
    report("Verify. output tsd.$a.$b.$tsd_size");
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
    report("Verify. output transposon.$a.$b.$tsd_size");
    open(OUT, "> $wd/$a/transposon.$a.$b.$tsd_size");
    open(IN, "$command $wd/$b/read/* |");
    while(<IN>){
	$count++;
	if ($count % 4 == 2){
	    next if /N/;
	    foreach $target (keys %target){
		if (/$target/){
		    print "$target{$target}\t$target\t$_";
		    print OUT "$target{$target}\t$target\t$_";
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
	       &canFork;
	       &report("making ref20.$tag.gz");
	       system("perl $0 a=$a,b=$b,ref=$ref,sub=sortPairFunc,tag=$tag,target=$a,wd=$wd &");
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
	       &report("sort pair.$tag in $b");
               system("perl $0 a=$a,b=$b,ref=$ref,sub=sortPairFunc,tag=$tag,target=$b,wd=$wd &");
	   }
       }
   }
   system("rm $wd/$a/tmp/p.* $wd/$b/tmp/p.*");
}

sub sortPairFunc{
    system("sort -S 1M -T $wd/$target/tmp $wd/$target/tmp/pair.$tag.tmp > $wd/$target/tmp/pair.$tag && rm $wd/$target/tmp/pair.$tag.tmp");
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

sub join{
    &canFork;
    $cmd = "perl $0 target=$a,sub=countHeadTail,type=head,a=$a,wd=$wd &";
    &report($cmd);
    system($cmd);
    &canFork;
    $cmd = "perl $0 target=$a,sub=countHeadTail,type=tail,a=$a,wd=$wd &";
    &report($cmd);
    system($cmd);
    &canFork;
    $cmd = "perl $0 target=$b,sub=countHeadTail,type=head,a=$a,wd=$wd &";
    &report($cmd);
    system($cmd);
    &canFork;
    $cmd = "perl $0 target=$b,sub=countHeadTail,type=tail,a=$a,wd=$wd &";
    &report($cmd);
    system($cmd);
}

sub countHeadTail{
    my (@row, $prev, $count);
    system("cat $wd/$target/tmp/*.$type | sort -S 1M -T $wd/$target/tmp > $wd/$target/tmp/$type");
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
    open(IN, "bash -c 'comm -3 <(zcat $wd/$a/tmp/$tag.gz) <(zcat $wd/$b/tmp/$tag.gz) '|");
    while(<IN>){
	chomp;
	@row = split('\t', $_);
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
		&canFork;
		$cmd = "perl $0 a=$a,b=$b,sub=commFunc,tag=$tag,tsd_size=$tsd_size,wd=$wd &";
		&report("common : $tag : $tsd_size");
		system($cmd);
	    }
	}
    }
}

sub mergeFunc{
    my ($num, @num, @sorted, $a, $b, $c, $filea, $fileb, $filec, @row);

    chdir "$wd/$target/tmp";
    
    opendir(DIR, "$wd/$target/tmp/");
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
	$now = `date`;
	chomp($now);
	system("bash -c 'join -a 1 -a 2 <(zcat $filea) <(zcat $fileb)' | awk '{print \$1 \"\t\" \$2 + \$3}' |gzip -c > $filec && rm $filea $fileb");
	$a += 2;
	$b += 2;
	$c ++;
    }

    open(IN, "zcat $filec |");
    open(OUT, "|gzip -c > $tag.gz");
    while(<IN>){
	@row = split;
	next if $row[0] =~/AAAAAAAAAA|CCCCCCCCCC|GGGGGGGGGG|TTTTTTTTTT/;
	if ($row[1] >= 5){
	    print OUT "$row[0]\n";
	}
    }
    close(IN);
    close(OUT);
    system("rm $tag.*.gz");
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
		&canFork;
		$cmd = "perl $0 target=$target,sub=mergeFunc,tag=$tag,a=$a,wd=$wd &";
		&report("perl $0 target=$target,sub=mergeFunc,tag=$tag,a=$a,wd=$wd");
		system($cmd);
	    }
	}
    }
}

sub countFunc{
    my ($taga, $tagb, $tagc, $tag, $count, $nuc);
    chdir "$wd/$target/tmp";
    open(IN, "split.$number");
    open(SEQ, "|sort -S 1M -T $wd/$target/tmp |uniq -c |awk '{print \$2 \"\t\" \$1}' > count.$number");
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
		open($fout, "|gzip -c > $tag.$number.gz");
	    }
	}
    }
    
    open(IN, "count.$number");
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
    system("rm count.$number");
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
	print "!!! hit CTACATGTCATTGATATACCTGTTAAATATATATACAAGC\n\n" if $seq eq "CTACATGTCATTGATATACCTGTTAAATATATATACAAGC"; 
    }
}

sub count{
    my $target = shift;
    my ($last, $i, $cmd);
    $tsd_size = 20 if $method eq "junctionMethod";
    opendir(DIR, "$wd/$target/tmp");
    foreach (readdir(DIR)){
	if (/^split/){
	    $number = (split('\.', $_))[1];
	    if ($last < $number){
		$last = $number;
	    }
	}
    }
    foreach ($i = 1; $i <= $last; $i++){
	&canFork;
	$cmd = "perl $0 target=$target,sub=countFunc,number=$i,a=$a,tsd_size=$tsd_size,wd=$wd &";
	&report($cmd);
	system($cmd);
    }
}

sub split{
    system("mkdir $wd/$target/tmp") if ! -e "$wd/$target/tmp";
    $number = 1;
    
    $command = "cat";
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
    open(OUT, "> $wd/$target/tmp/split.$number");
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
	    open(OUT, "> $target/tmp/split.$number");
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

sub canFork{
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

sub waitChild{
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

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    $message = "$now : $message\n";
    open(REPORT, ">> $wd/$a/log");
    flock(REPORT,2);
    print $message;
    print REPORT $message;
    close(REPORT);
}

sub bynumber{
    $a <=> $b;
}
