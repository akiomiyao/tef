#!/usr/bin/perl
#
# tef.pl - transposable element finder
#          A program for transposable elements detection from NGS reads.
#
# Copyright (C) 2020 National Agriculture and Food Research Organization. All Rights Reserved.
#
# License: refer to https://github.com/akiomiyao/tef
#
# Author: Akio Miyao
#


if ($ARGV[0] eq ""){
    print "
tef.pl - Transposable Element Finder

A program for detecting active transposable elements from NGS reads.

Usage Examples:
 perl tef.pl a=ttm2,ref=IRGSP1.0
 perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0
 perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,tsd_size=5,th=0.7
 ulimit -n 4096 && perl tef.pl a=ttm2,b=ttm5,ref=IRGSP1.0,debug=yes,max_process=8,sort_tmp=/mnt/ssd/tmp

Update: With the new feature, you can now omit the b parameter. When b is omitted, TEF
 will use the reference genome specified by ref as the second set of data. This means
 TEF will compare the NGS reads specified by a directly with the reference genome
 to detect transpositions.

Requirements: tef.pl requires a pair of NGS reads from different conditions. 
 For example, ttm2 and ttm5 are from regenerated rice individuals from callus.
 Save FASTQ files into the ttm2/read and ttm5/read directories.
 If NGS reads specified by 'b' are omitted, tef.pl compares the reads specified
 by 'a' with the reference genome sequence.

Options: Options should be specified with 'name=value' separated by a comma WITHOUT spaces.
 The default method is the junction method.
 If tsd_size is specified, tef.pl will run with the TSD method.

Reference Genome: If ref is specified, a gzipped FASTA file of the reference genome
 should be saved into the directory specified by ref. On the first run, a config
 file for the reference will be created in the ref directory.
 The config file contains a list of sequence names. If you do not want to include
 unassigned contigs, change the name to NOP, save the config file,
 and run tef.pl again.

Detection Methods: If ref is not specified, specific transpositions without
 insertion positions on the genome will be detected.
 This function works only with the TSD method.
 If tsd_size is not specified, detection will proceed using the junction method.
 If th (threshold) for the TSD method is not specified, the default value of 0.2 will
 be used. If a lot of noise is detected, try again with a higher value, e.g., th=0.7.

Temporary Directories: Because temporary directories in targets grow to a huge size,
 they will be deleted at the end of the analysis. If debug=yes is added to the
 options, temporary directories will not be deleted. Additionally, to save disk
 space, the compress=yes option can be used to compress temporary files.
 sort_tmp=directory_for_sort is the option for specifying a temporary directory
 for the sort command. If sort_tmp is specified to a fast disk, e.g., SSD,
 sorting will be accelerated. Setting all directories on an SSD disk is recommended.

Processing: For Linux, max_process is the number of CPU cores. For other OS,
 the default process number is 4. If you add the max_process option,
 e.g., max_process=8, tef.pl will use 8 cores.

Troubleshooting: If data in split and/or count.tsd_size are truncated, remove the
 directory in the target and then run again. TEF may open more than 1000 files simultaneously.
 If a file open error is displayed, please execute ulimit -n 4096 before running tef.pl.

Web Page: https://akiomiyao.github.io/tef/
Code: https://github.com/akiomiyao/tef
Paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05011-3
Question: https://github.com/akiomiyao/tef/issues

Author: Akio Miyao

";
    exit;
}

$start_time = time;

$s = {};
@nuc = ('A', 'C', 'G', 'T');
foreach $nuc (@nuc){
    my $nuca = $nuc;
    foreach $nuc (@nuc){
	my $nucb = $nuc;
	foreach $nuc (@nuc){
	    my $nucc = $nuc;
	    my $tag = $nuca . $nucb . $nucc;
	    push(@tag, $tag);
	}
    }
}

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

$processor = 20 if $processor > 20;
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
    system("rm -rf $wd/$a/child") if -e "$wd/$a/child";
    system("mkdir $wd/$a/tmp") if ! -e "$wd/$a/tmp";
    $b = $ref if $b eq "";
    system("mkdir $wd/$b/tmp") if ! -e "$wd/$b/tmp";
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
    system("rm -r $wd/$a/child $wd/$a/tmp $wd/$b/tmp") if ! $debug;
    system("rm $wd/$a/junction_method.all.$a.$b") if ! $debug;
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

    my $fragment_length = 20 + $tsd_size;
    if (&fragmentLength($a) != $fragment_length){
	&count($a);
    }
    if (&fragmentLength($b) != $fragment_length){
	&count($b);
    }
}

sub junctionMethod{
    &junctionSpecific;
    &join;
    &junctionFirstMap($a);
    &junctionFirstMap($b);
    &join;
    &junctionSecondMap($a);
    &junctionSecondMap($b);
    &join;
    &junctionSort($a);
    &junctionSort($b);
    &join;
    &junctionCandidate;
    &junctionCountCandidate;
    &toVcf;
}

sub tsdMethod{
    &tsdHeadTail;
    &verify;
    if ($ref ne ""){
	&map;
    }
}

sub toVcf{
    my ($tsd_size, $rnuc, $direction, $gt, $hcount, $tcount, $wcount, $total, $alt, $wt, $ad, $af);
    &log("toVcf : output vcf format file : junction_method.$a.$b.vcf");
    $timestamp = `date '+%Y-%m-%d %H:%M:%S %z'`;
    chomp($timestamp);
    $filedate = (split('\ ', $timestamp))[0];
    $filedate =~ y/-//d;
    open(IN, "$wd/$a/junction_method.genotype.$a.$b");
    open(OUT, "> $wd/$a/junction_method.$a.$b.vcf");
    print OUT "##fileformat=VCFv4.3
##fileDate=$filedate
##source=<PROGRAM=tif.pl,target=$target,reference=$ref>
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=MEINFO,Number=9,Type=String,Description=\"Movile element info of the form ME_HEAD_SEQ,ME_TAIL_SEQ,JUNCTION_POS_OF_HEAD,JUNCTION_POS_OF_TAIL,TSD_SIZE,TSD_SEQUENCE,DIRECTION,COUNT_OF_READS_WITH_JUNCTION_OF_HEAD,COUNT_OF_READS_WITH_JUNCTION_OF_TAIL,COUNT_OF_NOINSERTION_READS\">
##ALT=<ID=INS,Description=\"Insertion of a mobile element\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=.,Type=String,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Read Depth\">
##created=<TIMESTAMP=\"$timestamp\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$a\t$b\n";
    while(<IN>){
	chomp;
	@row = split;
	$rnuc = substr($row[7], length($row[7]) - 1, 1);
	$tsd_size = length($row[8]);
	if ($row[4] eq "f"){
	    $direction = "forward";
	}elsif($row[4] eq "r"){
	    $diredtion= "reverse";
	}
	if ($row[16] eq "H"){
	    $gt = "1/0";
	}elsif($row[16] eq "M"){
	    $gt = "1/1";
	}
	if ($row[0] eq $a){
	    $hcount = $row[10];
	    $tcount = $row[11];
	    $wcount = $row[12];
	    $total = $row[10] + $row[11] + $row[12];
	    $alt = $row[10] + $row[11];
	    $wt = $row[12];
	    $ad = "$wt,$alt";
	    $af = int($alt * 1000 / $total) / 1000;
	    print OUT "$row[1]\t$row[2]\t.\t$rnuc\t$rnuc<INS>\t.\t.\tGT=$gt;AF=$af;DP=$total;MEINFO=$row[5],$row[6],$row[2],$row[3],$tsd_size,$row[8],$direction,$hcount,$tcount,$wcount;SVTYPE=INS\tGT:AD:DP\t$gt:$ad:$total\t.\n";
	}else{
	    $hcount = $row[13];
	    $tcount = $row[14];
	    $wcount = $row[15];
	    $total = $row[13] + $row[14] + $row[15];
	    $alt = $row[13] + $row[14];
	    $wt = $row[15];
	    $ad = "$wt,$alt";
	    $af = int($alt * 1000 / $total) / 1000;
	    print OUT "$row[1]\t$row[2]\t.\t$rnuc\t$rnuc<INS>\t.\t.\tGT=$gt;AF=$af;DP=$total;MEINFO=$row[5],$row[6],$row[2],$row[3],$tsd_size,$row[8],$direction,$hcount,$tcount,$wcount;SVTYPE=INS\tGT:AD:DP\t.\t$gt:$ad:$total\n";
	}
    }
    close(IN);
    close(OUT);
}

sub junctionCandidate{
    opendir(REF, "$wd/$ref");
    foreach $file (sort readdir(REF)){
	if ($file =~ /^chr/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionCandidateFunc,chr=$file,sort_tmp=$sort_tmp &";
	    &log($cmd);
	    system($cmd);
	}
    }
    &join;
}

sub junctionCandidateFunc{
    open(CHR, "$wd/$ref/$chr");
    binmode(CHR);
    ($chrnum = $chr)=~ s/^chr//;
    open(OUT, "| gzip > $wd/$a/tmp/candidate.$chr.gz");
    open(IN, "cat $wd/$a/tmp/sorted.$chr.f $wd/$a/tmp/sorted.$chr.r |");
    while(<IN>){
	chomp;
	push(@dat, $_);
	next if $i++ < 40;
	$dat = shift(@dat);
	@row = split('\t', $dat);
	foreach (@dat){
	    @pos = split;
	    $size = $pos[0] - $row[0];
	    next if $size > 16;
	    next if $size == 0;
	    $size ++;
	    if ($row[1] eq "h"){
		if ($pos[1] eq "t"){
		    $htsd = substr($row[2], 20 - $size, $size);
		    $ttsd = substr($pos[2], 20, $size);
		    if ($htsd eq $ttsd){
			$head = substr($row[2], 20, 20);
			$hflanking = substr($row[2], 0, 20);
			$tail = substr($pos[2], 0, 20);
			$tflanking = substr($pos[2], 20, 20);
			seek(CHR, $row[0] - 21, 0);
			read(CHR, $wt, 40);
			next if $wt =~ /N/;
			print OUT "$row[2]\t$chrnum\t$row[0]\t$pos[0]\tah\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]
$pos[2]\t$chrnum\t$row[0]\t$pos[0]\tat\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]
$wt\t$chrnum\t$row[0]\t$pos[0]\taw\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]\n";
		    }
		}
	    }else{
		if ($pos[1] eq "h"){
		    $htsd = substr($pos[2], 20 - $size, $size);
		    $ttsd = substr($row[2], 20, $size);
		    if ($htsd eq $ttsd){
			$head = substr($pos[2], 20, 20);
			$hflanking = substr($pos[2], 0, 20);
			$tail = substr($row[2], 0, 20);
			$tflanking = substr($row[2], 20, 20);
			seek(CHR, $pos[0] - 21, 0);
			read(CHR, $wt, 40);
			next if $wt =~ /N/;
			print OUT "$pos[2]\t$chrnum\t$pos[0]\t$row[0]\tah\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]
$row[2]\t$chrnum\t$pos[0]\t$row[0]\tat\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]
$wt\t$chrnum\t$pos[0]\t$row[0]\taw\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]\n";
		    }
		}
	    }
	}
    }
    close(IN);
    system("rm $wd/$a/tmp/sorted.$chr.f $wd/$a/tmp/sorted.$chr.r");
    open(IN, "cat $wd/$b/tmp/sorted.$chr.f $wd/$b/tmp/sorted.$chr.r |");
    while(<IN>){
	chomp;
	push(@dat, $_);
	next if $i++ < 40;
	$dat = shift(@dat);
	@row = split('\t', $dat);
	foreach (@dat){
	    @pos = split;
	    $size = $pos[0] - $row[0];
	    next if $size > 16;
	    next if $size == 0;
	    $size ++;
	    if ($row[1] eq "h"){
		if ($pos[1] eq "t"){
		    $htsd = substr($row[2], 20 - $size, $size);
		    $ttsd = substr($pos[2], 20, $size);
		    if ($htsd eq $ttsd){
			$head = substr($row[2], 20, 20);
			$hflanking = substr($row[2], 0, 20);
			$tail = substr($pos[2], 0, 20);
			$tflanking = substr($pos[2], 20, 20);
			seek(CHR, $row[0] - 21, 0);
			read(CHR, $wt, 40);
			next if $wt =~ /N/;
			print OUT "$row[2]\t$chrnum\t$row[0]\t$pos[0]\tbh\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]
$pos[2]\t$chrnum\t$row[0]\t$pos[0]\tbt\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]
$wt\t$chrnum\t$row[0]\t$pos[0]\tbw\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[5]\n";
		    }
		}
	    }else{
		if ($pos[1] eq "h"){
		    $htsd = substr($pos[2], 20 - $size, $size);
		    $ttsd = substr($row[2], 20, $size);
		    if ($htsd eq $ttsd){
			$head = substr($pos[2], 20, 20);
			$hflanking = substr($pos[2], 0, 20);
			$tail = substr($row[2], 0, 20);
			$tflanking = substr($row[2], 20, 20);
			seek(CHR, $pos[0] - 21, 0);
			read(CHR, $wt, 40);
			next if $wt =~ /N/;
			print OUT "$pos[2]\t$chrnum\t$pos[0]\t$row[0]\tbh\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]
$row[2]\t$chrnum\t$pos[0]\t$row[0]\tbt\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]
$wt\t$chrnum\t$pos[0]\t$row[0]\tbw\t$head\t$tail\t$hflanking\t$htsd\t$tflanking\t$row[8]\n";
		    }
		}
	    }
	}
    }
    close(IN);
    open(IN, "$wd/$b/tmp/sorted.$chr.f $wd/$b/tmp/sorted.$chr.r");
    close(OUT);
}

sub junctionCountCandidate{
    my (@chr, $chr);
    opendir(DIR, "$wd/$ref");
    foreach (sort readdir(DIR)){
	if (/^chr/){
	    ($chr = $_) =~ s/^chr//;
	    push(@chr, $chr);
	}
    }
    &log("junctionCountCandidate: $a : sorting candidates");
    foreach $tag (@tag){
	open($tag, "|gzip > $wd/$a/tmp/tmpa.$tag.gz ");
    }
    open(IN, "zcat $wd/$a/tmp/candidate.*.gz |");
    while(<IN>){
	$tag = substr($_, 0, 3);
	print $tag $_;
    }
    close(IN);
    $org_processor = $processor;
    $processor = 4;
    foreach $tag (@tag){
	close($tag);
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionSortCandidateFunc,tag=$tag,sort_tmp=$sort_tmp &";
	&log("$cmd");
	system($cmd);
    }
    &join;
    $processor = $org_processor;
    
    foreach $tag (@tag){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionCountCandidateFunc,tag=$tag,sort_tmp=$sort_tmp &";
	&log("$cmd");
	system($cmd);
    }
    &join;

    foreach $tag (@tag){
	open(IN, "$wd/$a/tmp/count.$tag");
	foreach $chr (@chr){
	    open($chr, "> $wd/$a/tmp/count.$chr.$tag");
	}
	while(<IN>){
	    $chr = (split)[0];
	    print $chr $_;
	}
	foreach $chr (@chr){
	    close($chr);
	}
	close($tag);
	system("rm $wd/$a/tmp/count.$tag") if ! $debug;
    }
    system("rm $wd/$a/tmp/candidate.*") if ! $debug;

    &log("junctionCountCandidate : making junction_method.summary.$a.$b and junction_method.all.$a.$b");
    foreach $chr (@chr){
	&log("junctionCountCandidate : making junction_method.genotype.$a.$b.$chr");
	open(OUT, "| sort  -S 1M -T $sort_tmp > $wd/$a/tmp/junction_method.genotype.tmp.$chr");
	open(IN, "cat $wd/$a/tmp/count.$chr.* | sort  -S 1M -T $sort_tmp | uniq |");
	while(<IN>){
	    chomp;
	    @row = split;
	    if ($row[3] =~ /h$/){
		$ah = $row[10];
		$bh = $row[11];
	    }elsif($row[3] =~ /t$/){
		$at = $row[10];
		$bt = $row[11];
	    }elsif($row[3] =~ /w$/){
		$aw = $row[10];
		$bw = $row[11];
		$genotype = "";
		if ($ah > 0 and $at > 0 and $bh == 0 and $bt == 0 and $bw > 0){
		    if ($aw == 0){
			$genotype = "M";
		    }else{
			$genotype = "H";
		    }
		}elsif($ah == 0 and $at == 0 and $aw > 0 and $bh > 0 and $bt > 0){
		    if ($bw == 0){
			$genotype = "M";
		    }else{
			$genotype = "H";
		    }
		}
		if ($genotype){
		    if ($row[3] =~ /^a/){
			$target = $a;
		    }else{
			$target = $b;
		    }
		    $row[0] = "0000$row[0]";
		    $row[0] = substr($row[0], length($row[0]) - 4, 4);
		    $row[1] = "0000000000$row[1]";
		    $row[1] = substr($row[1], length($row[1]) - 10, 10);
		    $row[2] = "0000000000$row[2]";
		    $row[2] = substr($row[2], length($row[2]) - 10, 10);
		    print OUT "$target\t$row[0]\t$row[1]\t$row[2]\t$row[9]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$ah\t$at\t$aw\t$bh\t$bt\t$bw\t$genotype\n";
		}
		undef $ah;
		undef $at;
		undef $at;
		undef $aw;
		undef $bh;
		undef $bt;
		undef $bw;
	    }
	}
	close(IN);
	close(OUT);
    }
    system("rm $wd/$a/tmp/count.*") if ! $debug;

    &log("junctionCountCandidate : making junction_method.all.$a.$b");
    open(OUT, "> $wd/$a/junction_method.all.$a.$b");
    open(IN, "cat $wd/$a/tmp/junction_method.genotype.tmp.* |sort -S 1M -T $sort_tmp |");
    while(<IN>){
	chomp;
	@row = split;
	$row[1] =~ s/^0*// ;
	$row[2] += 0;
	$row[3] += 0;
	print OUT join("\t", @row) . "\n";
    }
    close(IN);
    close(OUT);

    &log("junctionCountCandidate : making telist");
    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/telist");
    open(IN, "$wd/$a/junction_method.all.$a.$b");
    while(<IN>){
	@row = split;
	next if length($row[8]) > 10;
	print OUT "$row[5]\t$row[6]\t$row[0]\t$row[8]\n";
    }
    close(IN);
    close(OUT);

    &log("junctionCountCandidate : making telist.count");
    open(OUT, "|uniq -c > $wd/$a/tmp/telist.count");
    open(IN, "$wd/$a/tmp/telist");
    while(<IN>){
	chomp;
	@row = split;
	print OUT "$row[0] $row[1]\n";
    }
    close(IN);
    close(OUT);

    open(HEAD, "| sort -k 1 -k 2 -n -S 1M -T $sort_tmp > $wd/$a/tmp/telist.head");
    open(TAIL, "| sort -k 1 -k 2 -n -S 1M -T $sort_tmp > $wd/$a/tmp/telist.tail");
    open(IN, "$wd/$a/tmp/telist.count");
    while(<IN>){
	chomp;
	@row = split;
	next if $row[0] < 2;
	$ht{"$row[1] $row[2]"} = " ";
    }
    close(IN);

    open(IN, "$wd/$a/tmp/telist");
    while(<IN>){
	chomp;
	@row = split;
	if ($ht{"$row[0] $row[1]"} ne ""){
	    $ht{"$row[0] $row[1]"} .= "$row[3] ";
	}
    }
    close(IN);

    foreach (sort keys %ht){
	$tsd = $ht{$_};
	$tsd =~ s/^\ //;
	$tsd =~ s/\ $//;
	$flag = 1;
	@tsd = split(' ', $tsd);
	$flag = 0 if $#tsd < 1;
	$size = length($tsd[0]);
	for ($i = 1; $i <= $#tsd; $i++){
	    if (length($tsd[$i]) != $size){
		$flag = 0;
	    }
	}
	if ($flag){
	    my ($head, $tail) = split;
	    $tsd_count = $#tsd + 1;
	    print HEAD "$head $tsd_count $tail\n";
	    print TAIL "$tail $tsd_count $head\n";
	}
    }
    close(IN);
    close(HEAD);
    close(TAIL);

    %select = ();
    open(TAIL, "$wd/$a/tmp/telist.tail");
    while(<TAIL>){
	chomp;
	@row = split;
	$select{$row[0]} = $row[2];
    }
    close(TAIL);

    %exclude = ();
    open(TAIL, "$wd/$a/tmp/telist.tail");
    while(<TAIL>){
	chomp;
	@row = split;
	if ($select{$row[0]} ne $row[2]){
	    $exclude{$row[2]} = $row[0];
	}
    }
    close(TAIL);
    %select = ();
    open(HEAD, "$wd/$a/tmp/telist.head");
    while(<HEAD>){
	chomp;
	@row = split;
	$select{$row[0]} = $row[2];
    }
    close(HEAD);

    foreach $head (sort keys %select){
	$tail = $select{$head};
	if ($exclude{$head} ne $tail){
	    $tsd = $ht{"$head $tail"};
	    @tsd = split("\ ", $tsd);
	    @std = split('', $tsd[0]);
	    $similar = 0;
	    for($i = 1; $i <= $#tsd; $i++){
		$hit = 0;
		@tg = split('', $tsd[$i]);
		for($j = 0; $j <= $#std; $j++){
		    if ($std[$j] eq $tg[$j]){
			$hit++;
		    }
		}
		if ($hit / $#std > 0.6){
		    $similar++;
		}
	    }
	    if ($similar / ($#tsd + 1) < 0.5){
		$final{$head} = $tail;
	    }
	}
    }

    &log("junctionCountCandidate : output junction_method.summary.$a.$b");
    $acount = 0;
    $bcount = 0;
    $tea = "";
    $teb = "";
    open(OUT, "> $wd/$a/junction_method.summary.$a.$b");
    open(IN, "$wd/$a/tmp/telist");
    while(<IN>){
	chomp;
	@row = split;
	if ($prev[0] ne $row[0] or $prev[1] ne $row[1]){
	    if ($prev[0] ne "" and $final{$prev[0]} eq $prev[1]){
		$tea =~ s/\ $//;
		$teb =~ s/\ $//;
		print OUT "$prev[0]\t$prev[1]\t$acount\t$bcount\t$tea\t$teb\n";
		print "$prev[0]\t$prev[1]\t$acount\t$bcount\t$tea\t$teb\n";
	    }
	    $acount = 0;
	    $bcount = 0;
	    $tea = "";
	    $teb = "";
	    if ($row[2] eq $a){
		$acount ++;
		$tea .= "$row[3] "
	    }else{
		$bcount ++;
		$teb .= "$row[3] "
	    }
	}else{
	    if ($row[2] eq $a){
		$acount ++;
		$tea .= "$row[3] "
	    }else{
		$bcount ++;
		$teb .= "$row[3] "
	    }
	}
	@prev = @row;
    }
    close(IN);
    if ($row[0] ne $prev[0] and $final{$row[0]} eq $row[1]){
	if ($row[2] eq $a){
	    $acount ++;
	    $tea .= "$row[3] "
	}else{
	    $bcount ++;
	    $teb .= "$row[3] "
	}
	$tea =~ s/\ $//;
	$teb =~ s/\ $//;
	print OUT "$row[0]\t$row[1]\t$acount\t$bcount\t$tea\t$teb\n";
	print "$row[0]\t$row[1]\t$acount\t$bcount\t$tea\t$teb\n";
    }
    close(OUT);
    
    &log("junctionCountCandidate: making junction_method.genotype.$a.$b");
    open(OUT, "> $wd/$a/junction_method.genotype.$a.$b");
    open(IN, "$wd/$a/junction_method.all.$a.$b");
    while(<IN>){
	@row = split;
	if ($final{$row[5]} eq $row[6]){
	    print OUT;
	}
    }
    close(IN);
    close(OUT);

    open(OUT, "|sort  -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/te.list");
    open(IN, "$wd/$a/junction_method.summary.$a.$b");
    while(<IN>){
	@row = split;
	next if $row[0] eq "Head";
	print OUT "$row[0]\thead
$row[1]\ttail\n";
    }
    close(IN);
    close(OUT);

    open(IN, "$wd/$a/tmp/te.list");
    while(<IN>){
	$tag = substr($_, 0, 3);
	if ($tag ne $prev){
	    open(OUT, "> $wd/$a/tmp/telist.$tag");
	    push(@tetag, $tag);
	}
	print OUT;
	$prev = $tag;
    }
    close(IN);
    close(OUT);

    foreach $tag (@tetag){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionTePosFunc,tag=$tag,sort_tmp=$sort_tmp &";
	&log("$cmd");
	system($cmd);
    }
    &join;

    open(OUT, "|sort -S 1M -T $sort_tmp >  $wd/$a/tmp/tmp.tepos");
    open(IN, "cat $wd/$a/tmp/tepos.* |");
    while(<IN>){
	chomp;
	@row = split;
	$row[2] = "0000$row[2]";
	$row[2] = substr($row[2], length($row[2]) - 4, 4);
	$row[3] = "0000000000$row[3]";
	$row[3] = substr($row[3], length($row[3]) - 10, 10);
	print OUT "$row[2]\t$row[3]\t$row[0]\t$row[1]\t$row[4]\n";
    }
    close(IN);
    close(OUT);

    open(OUT, "> $wd/$a/junction_method.tepos.$a.$b");
    open(IN, "$wd/$a/tmp/tmp.tepos");
    while(<IN>){
	chomp;
	@row = split;
	$row[0] =~ s/^0*//;
	$row[1] =~ s/^0*//;
	print OUT join("\t", @row) . "\n";
    }
    close(IN);
    close(OUT);
}

sub junctionSortCandidateFunc{
    system ("zcat $wd/$a/tmp/tmpa.$tag |sort  -S 1M -T $sort_tmp |gzip > $wd/$a/tmp/tmp.$tag.gz");
    system ("rm $wd/$a/tmp/tmpa.$tag.gz");
}

sub junctionCountCandidateFunc{
    $cmd = "bash -c 'join -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.2 -e 0 <(zcat $wd/$a/tmp/tmp.$tag.gz) <(zcat $wd/$a/count.20/$tag.gz) > $wd/$a/tmp/tmpa.$tag'";
    system($cmd);
#    system("zcat $wd/$a/count.20/$tag.gz |join -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 2.2 -e 0 $wd/$a/tmp/tmp.$tag - > $wd/$a/tmp/tmpa.$tag");

#    $cmd = "bash -c 'join -a 1 -o 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.2 -e 0 <(zcat $wd/$a/tmp/tmpa.$tag.gz) <(zcat $wd/$b/count.20/$tag.gz) > $wd/$a/tmp/count.$tag'";
    system($cmd);

    system("zcat $wd/$b/count.20/$tag.gz |join -a 1 -o 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.2 -e 0 $wd/$a/tmp/tmpa.$tag - > $wd/$a/tmp/count.$tag");
    system("rm $wd/$a/tmp/tmp.$tag.gz $wd/$a/tmp/tmpa.$tag") if ! $debug;
}

sub junctionTePosFunc{
    return if ! -e "$wd/$a/tmp/telist.$tag";
    system("zcat $wd/$ref/ref20.$tag.gz | join $wd/$a/tmp/telist.$tag - > $wd/$a/tmp/tepos.$tag");
}

sub junctionSort{
    my $target = shift;
    my @chr;
    &log("junctionSort : Sorting filtered data : $target");
    opendir(REF, "$wd/$ref");
    foreach $file (sort readdir(REF)){
	if ($file =~ /^chr/){
	    push(@chr, $file);
	    $filename = $file . ".f";
	    open($filename, "> $wd/$target/tmp/$filename");
	    $filename = $file . ".r";
	    open($filename, "> $wd/$target/tmp/$filename");
	}
    }
    open(IN, "cat $wd/$target/tmp/map.* |");
    my ($head, $tail);
    while(<IN>){
	chomp;
	@row = split;
	$head = "chr$row[1]";
	$tail = "chr$row[4]";
	if($row[3] eq "f"){
	    $pos = $row[2] + 19;
	}else{
	    $pos = $row[2] - 19;
	}
	$head = $head . "." . $row[3];
	print $head "$pos\th\t$row[0]\t$row[1]\t$pos\t$row[3]\t$row[4]\t$row[5]\t$row[6]\n";
	$tail = $tail . "." . $row[6];
	print $tail "$row[5]\tt\t$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\n";
    }
    close(IN);
    system("rm $wd/$target/tmp/map.*") if ! $debug;
    foreach $file (@chr){
	close("$file.f");
	close("$file.r");
    }
    $org_processor = $processor;
    $processor = 2 if -s "$wd/$target/tmp/$chr[0]" > 10000000000;
    foreach $file (@chr){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSortFunc,chr=$file,sort_tmp=$sort_tmp &";
	&log($cmd);
	system($cmd);
    }
    $processor = $org_processor;
}

sub junctionSortFunc{
    $cmd = "sort -k 1 -n -S 1M -T $sort_tmp $wd/$target/tmp/$chr.f > $wd/$target/tmp/sorted.$chr.f";
    system($cmd);
    system("rm $wd/$target/tmp/$chr.f") if ! $debug;
    $cmd = "sort -k 1 -n -S 1M -T $sort_tmp $wd/$target/tmp/$chr.r > $wd/$target/tmp/sorted.$chr.r";
    system($cmd);
    system("rm $wd/$target/tmp/$chr.r") if ! $debug;
}

sub junctionSecondMap{
    my $target = shift;
    &log("junctionSecondMap : $target");
    foreach $tag (@tag){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionSecondMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
	&log("$cmd");
	system($cmd);
    }
}

sub junctionSecondMapFunc{
    system("rm $wd/$target/tmp/map.$tag") if -e "$wd/$target/tmp/map.$tag";
    $cmd = "cat $wd/$target/tmp/first.$tag.* | sort -S 1M -T $sort_tmp > $wd/$target/tmp/second.$tag";
    system($cmd);
    system("rm $wd/$target/tmp/first.$tag.*") if ! $debug;
    open(OUT, "> $wd/$target/tmp/map.$tag");
    open(IN, "zcat $wd/$ref/ref20.$tag.gz |join $wd/$target/tmp/second.$tag - |");
    while(<IN>){
	chomp;
        @row = split;
	if ($prev ne $row[1] and abs($row[3] - $row[6]) > 25){
	    print OUT "$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\n";
	    $prev = $row[1];
	}
    }
    close(IN);
    system("rm $wd/$target/tmp/second.$tag") if ! $debug;
}

sub junctionFirstMap{
    my $target = shift;
    &log("junctionFirstMap : $target");
    foreach $tag (@tag){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$target,sub=junctionFirstMapFunc,tag=$tag,sort_tmp=$sort_tmp &";
	&log($cmd);
	system($cmd);
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
    system("rm $wd/$target/tmp/specific.$tag") if ! $debug;
    
    open(IN, "zcat $wd/$ref/ref20.$tag.gz |join $wd/$target/tmp/target.$tag - |");
    open(OUT, "> $wd/$target/tmp/first.$tag");
    while(<IN>){
	@row = split;
	$seq = substr($row[1], 20, 20);
	print OUT "$seq\t$row[1]\t$row[2]\t$row[3]\t$row[4]\n" if $prev ne $seq;
	$prev = $seq;
    }
    close(IN);
    close(OUT);
    system("rm $wd/$target/tmp/target.$tag") if ! $debug;
    
    foreach $subtag (@tag){
	open($subtag, "> $wd/$target/tmp/first.$subtag.$tag");
    }
    open(IN, "$wd/$target/tmp/first.$tag");
    while(<IN>){
	$subtag = substr($_, 0, 3);
	print $subtag $_;
    }
    close(IN);
    foreach $subtag (@tag){
	close($subtag);
    }
    system("rm $wd/$target/tmp/first.$tag") if ! $debug;
}

sub junctionSpecific{
    foreach $tag (@tag){
	&monitorWait;
	$cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=junctionSpecificFunc,tsd_size=$tsd_size,tag=$tag,sort_tmp=$sort_tmp &";
	&log($cmd);
	system($cmd);
    }
}

sub junctionSpecificFunc{
    $cmd = "bash -c 'join -v 1 <(zcat $wd/$a/count.$tsd_size/$tag.gz) <(zcat $wd/$b/count.$tsd_size/$tag.gz) > $wd/$a/tmp/specific.$tag.tmp'";
    system($cmd);
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
    system($cmd);
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
    foreach $tag (@tag){
	open($tag, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/mapquery.$tag")
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
    foreach $tag (@tag){
	close($tag)
    }
    foreach $tag (@tag){
	if (-s "$wd/$a/tmp/mapquery.$tag" > 0){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=mapQuery,tag=$tag,sort_tmp=$sort_tmp &";
	    &log($cmd);
	    system($cmd);
	}
    }
    &join;

    open(IN, "cat $wd/$a/tmp/map.* |");
    while(<IN>){
	next if /TATATATATATATATATA/;
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
	print OUT join("\t", @row) . "\n";
    }
    close(IN);
    close(OUT);
}

sub mapQuery{
    $cmd = "zcat $wd/$ref/ref20.$tag.gz | join - $wd/$a/tmp/mapquery.$tag > $wd/$a/tmp/map.$tag";
    system("$cmd");
}

sub verify{
    opendir(DIR, "$a/split");
    foreach (sort readdir(DIR)){
	if (/^split/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$a,sub=verifyFunc,tsd_size=$tsd_size,file=$_,sort_tmp=$sort_tmp &";
	    &log("$cmd");
	    system($cmd);
	}
    }
    closedir(DIR);
    opendir(DIR, "$b/split");
    foreach (sort readdir(DIR)){
	if (/^split/){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,target=$b,sub=verifyFunc,tsd_size=$tsd_size,file=$_,sort_tmp=$sort_tmp &";
	    &log($cmd);
	    system($cmd);
	}
    }
    closedir(DIR);
    &join;

    &log("verify : making query sequences");
    opendir(DIR, "$a/tmp");
    foreach (sort readdir(DIR)){
	if (/verify.split/){
	    open(IN, "$a/tmp/$_");
	    while(<IN>){
		chomp;
		@row = split;
		$s->{$row[0]}{$row[1]}{$row[3]}{$row[2]} = 1;
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
		$s->{$row[0]}{$row[1]}{$row[3]}{$row[2]} = 1;
	    }
	    close(IN);
	}
    }

    open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tmp/query");
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    }else{
	open(IN, "cat $wd/$a/tmp/te.candidate $wd/$b/tmp/te.candidate |");
    }
    while(<IN>){
	next if /ACACACACACAC/;
	next if /AGAGAGAGAGAG/;
	next if /ATATATATATAT/;
	next if /CGCGCGCGCGCG/;
	next if /CTCTCTCTCTCT/;
	next if /GTGTGTGTGTGT/;
	@row = split;
	foreach $htsd (sort keys %{$s->{$row[0]}{H}}){
	    foreach $upstream (sort keys %{$s->{$row[0]}{H}{$htsd}}){
		foreach $downstream (sort keys %{$s->{$row[1]}{T}{$htsd}}){
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
    
    foreach $tag (@tag){
	open($tag, "> $wd/$a/tmp/query.$tag")
    }
    open(IN, "$wd/$a/tmp/query");
    while(<IN>){
	$tag = substr($_, 0, 3);
	print $tag $_;

    }
    close(IN);
    foreach $tag (@tag){
	close($tag)
    }

    foreach $tag (@tag){
	if (-s "$wd/$a/tmp/query.$tag" > 0){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=countQuery,tag=$tag,sort_tmp=$sort_tmp &";
	    &log($cmd);
	    system($cmd);
	}
    }
    &join;

    &log("verify : making verify file");
    open(IN, "cat $wd/$a/tmp/verify.count.* |");
    while(<IN>){
	chomp;
	@row = split;
	$a{$row[0]} += $row[1];
	
    }
    open(IN, "cat $wd/$b/tmp/verify.count.* |");
    while(<IN>){
	chomp;
	@row = split;
	$b{$row[0]} += $row[1];
    }

    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
	open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/tsd_method.verify.$a.$b.$tsd_size");
    }else{
	open(IN, "cat $wd/$a/tmp/te.candidate $wd/$b/tmp/te.candidate |");
	open(OUT, "|sort -S 1M -T $sort_tmp |uniq > $wd/$a/junction_method.verify.$a.$b");
    }
    while(<IN>){
	next if /ACACACACACAC/;
	next if /AGAGAGAGAGAG/;
	next if /ATATATATATAT/;
	next if /CGCGCGCGCGCG/;
	next if /CTCTCTCTCTCT/;
	next if /GTGTGTGTGTGT/;
	@row = split;
	foreach $htsd (sort keys %{$s->{$row[0]}{H}}){
	    foreach $upstream (sort keys %{$s->{$row[0]}{H}{$htsd}}){
		foreach $downstream (sort keys %{$s->{$row[1]}{T}{$htsd}}){
		    $tsdsize = length($htsd);
		    $wildtype = $upstream . substr($downstream, $tsdsize, 20);
		    $head = $upstream . $row[0];
		    $downfl = substr($downstream, 0, 20);
		    $tail = $row[1] . $downfl;
		    if (($a{$head} > 0 and $a{$tail} > 0 and $b{$head} == 0 and $b{$tail} == 0 and $b{$wildtype} > 2) or ($a{$head} == 0 and $a{$tail} == 0 and $b{$head} > 0 and $b{$tail} > 0 and $a{$wildtype} > 2)){
			$result = "$row[0]\t$row[1]\t$upstream\t$htsd\t$downfl\t$a{$head}\t$a{$tail}\t$a{$wildtype}\t$b{$head}\t$b{$tail}\t$b{$wildtype}\n";
			print OUT $result;
			if ($a{$head} > 0){
			    $s->{tecount}{"$row[0]\t$row[1]"}{$a} ++;
			    $s->{tetsd}{"$row[0]\t$row[1]"}{$a}{$htsd} ++;
			}
			if ($b{$head} > 0){
			    $s->{tecount}{"$row[0]\t$row[1]"}{$b} ++;
			    $s->{tetsd}{"$row[0]\t$row[1]"}{$b}{$htsd} ++;
			}
		    }
		}
	    }
	}
    }
    close(IN);
    close(OUT);

    &log("verify : making summay file");
    if ($method eq "TSD"){
	open(OUT, "> $wd/$a/tsd_method.summary.$a.$b.$tsd_size");
    }else{
	open(OUT, "> $wd/$a/junction_method.summary.$a.$b");
    }
    print OUT "Head\tTail\t$a\t$b\n";
    foreach $te (sort keys %{$s->{tecount}}){
	$totala = $s->{tecount}{$te}{$a};
	$totalb = $s->{tecount}{$te}{$b};
	$counta = 0;
	$countb = 0;
	foreach (sort keys %{$s->{tetsd}{$te}{$a}}){
	    $counta ++;
	}
	foreach (sort keys %{$s->{tetsd}{$te}{$b}}){
	    $countb ++;
	}
	print OUT "$te\t$counta\t$countb\n";
    }
    close(OUT);
}

sub countQuery{
    $cmd = "zcat $wd/$a/count.20/$tag.gz | join - $wd/$a/tmp/query.$tag > $wd/$a/tmp/verify.count.$tag";
    system("$cmd");
    
    $cmd = "zcat $wd/$b/count.20/$tag.gz | join - $wd/$a/tmp/query.$tag > $wd/$b/tmp/verify.count.$tag";
    system("$cmd");
}

sub verifyFunc{
    my (@row, $tsdsize, %seq, %tsdsize, $length, $seq, $tsd, $upstream, $downstream);
    if ($method eq "TSD"){
	open(IN, "$wd/$a/tsd_method.pair.$a.$b.$tsd_size");
    }else{
	open(IN, "cat $wd/$a/tmp/te.candidate $wd/$b/tmp/te.candidate |");
    }
    while(<IN>){
	next if /ACACACACACAC/;
	next if /AGAGAGAGAGAG/;
	next if /ATATATATATAT/;
	next if /CGCGCGCGCGCG/;
	next if /CTCTCTCTCTCT/;
	next if /GTGTGTGTGTGT/;
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
	$length = length($_);
	for($i = 0; $i <= $length - 20; $i++){
	    $seq = substr($_, $i, 20);
	    if ($seq{$seq}){
		foreach $tsdsize (sort keys %{$s->{tsdsize}{$seq}}){		
		    if ($seq{$seq} eq "H"){
			$upstream = substr($_, $i - 20, 20);
			next if $upstream =~ /N/;
			if(length($upstream) == 20){
			    $tsd = substr($_, $i - $tsdsize, $tsdsize);
			    print OUT "$seq\tH\t$upstream\t$tsd\n";
			}	    
		    }elsif($seq{$seq} eq "T"){
			$downstream = substr($_, $i + 20, 20 + $tsdsize);
			next if $downstream =~ /N/;
			if(length($downstream) == 20 + $tsdsize){
			    $tsd = substr($_, $i +20, $tsdsize);
			    print OUT "$seq\tT\t$downstream\t$tsd\n";
			}
		    }
		}
	    }
	}
	$complement = &complement($_);
	for($i = 0; $i <= $length - 20; $i++){
	    $seq = substr($complement, $i, 20);
	    if ($seq{$seq}){
		foreach $tsdsize (sort keys %{$s->{tsdsize}{$seq}}){		
		    if ($seq{$seq} eq "H"){
			$upstream = substr($_, $i - 20, 20);
			next if $upstream =~ /N/;
			if(length($upstream) == 20){
			    $tsd = substr($_, $i - $tsdsize, $tsdsize);
			    print OUT "$seq\tH\t$upstream\t$tsd\n";
			}	    
		    }elsif($seq{$seq} eq "T"){
			$downstream = substr($_, $i + 20, 20 + $tsdsize);
			next if $downstream =~ /N/;
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
		system($cmd);
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
		system($cmd);
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
	system($cmd);
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
    my ($aa, $ab, $ac, $ad, $ae, $tag, $count, $lines, $seq, $file, $outfile, $count_file);
    my $target = shift;
    &log("count : $target : making count.$tsd_size");

    $die_comment = "Can not open file.
Please inclease limit for current session.
ulimit -n 4096
e.g. ulimit -n 4096 && perl tef.pl a=sampleA,b=sampleB,ref=TAIR10\n";

    if (! -d "$wd/$target/count.$tsd_size"){
	system("mkdir $wd/$target/count.$tsd_size");
    }

    if (! -d "$wd/$target/tmp"){
	system("mkdir $wd/$target/tmp");
    }
    
    foreach $aa (@nuc){
	foreach $ab (@nuc){
	    foreach $ac (@nuc){
		foreach $ad (@nuc){
		    foreach $ae (@nuc){
			$tag = $aa . $ab. $ac . $ad . $ae;
			open($tag, "|gzip > $wd/$target/tmp/$tag.gz") || die $die_comment;
		    }
		}
	    }
	}
    }
    if ($target eq $ref){
	opendir(DIR, "$wd/$target");
	foreach $file (sort readdir(DIR)){
	    next if $file !~ /^chr/;
	    open(CHR, "$target/$file") || die $die_comment;
	    binmode(CHR);
	    &log("count : $target : $file");
	    $pos = 0;
	    while(1){
		seek(CHR, $pos++, 0);
		read(CHR, $seq, 20 + $tsd_size);
		last if length($seq) !=  20 + $tsd_size;
		next if $seq =~ /N/;
		$tag = substr($seq, 0, 5);
		print $tag "$seq\n";
		$seq = &complement($seq);
		$tag = substr($seq, 0, 5);
		print $tag "$seq\n";
	    }
	    close(CHR);
	}
    }else{
	opendir(DIR, "$wd/$target/read");
	foreach $file (sort readdir(DIR)){
	    next if $file =~ /^\./;
	    &log("count : $target : $file");
	    $file = "$wd/$target/read/" . $file;
	    if ($file =~ /gz$/){
		open(IN, "zcat $file |") || die $die_comment;;
	    }elsif ($file =~ /bz2$/){
		open(IN, "bzcat $file |") || die $die_comment;;
	    }elsif ($file =~ /xz$/){
		open(IN, "xzcat $file |") || die $die_comment;;
	    }else{
		open(IN, $file) || die $die_comment;
	    }
	    while(<IN>){
		$count = 0 if $count++ == 3;
		if ($count == 2){
		    $lines ++;
		    if ($lines % 1000000 == 0){
			&log("count : $target : $lines reads have been selected");
		    }
		    chomp;
		    $length = length($_);
		    for($i = 0; $i < $length; $i++){
			$seq = substr($_, $i, 20 + $tsd_size);
			last if length($seq) != 20 + $tsd_size;
			next if $seq =~ /N/;
			$tag = substr($seq, 0, 5);
			print $tag "$seq\n";
		    }
		    $complement = &complement($_);
		    for($i = 0; $i < $length; $i++){
			$seq = substr($complement, $i, 20 + $tsd_size);
			last if length($seq) != 20 + $tsd_size;
			next if $seq =~ /N/;
			$tag = substr($seq, 0, 5);
			print $tag "$seq\n";
		    }
		}
	    }
	}
	close(IN);
    }
    close(DIR);

    foreach $aa (@nuc){
	foreach $ab (@nuc){
	    foreach $ac (@nuc){
		foreach $ad (@nuc){
		    foreach $ae (@nuc){
			$tag = $aa . $ab. $ac . $ad . $ae;
			close($tag);
		    }
		}
	    }
	}
    }
    
    foreach $aa (@nuc){
	foreach $ab (@nuc){
	    foreach $ac (@nuc){
		foreach $ad (@nuc){
		    foreach $ae (@nuc){
			$tag = $aa . $ab. $ac . $ad . $ae;
			&log("count : $target : counting $tag"); 
			&monitorWait;
			$cmd = "perl $0 a=$a,sub=countFunc,ref=$ref,target=$target,tag=$tag,tsd_size=$tsd_size,sort_tmp=$sort_tmp &";
			system($cmd);
		    }
		}
	    }
	}
    }
    &join;
    foreach $aa (@nuc){
	foreach $ab (@nuc){
	    foreach $ac (@nuc){
		$tag = $aa . $ab . $ac;
		&log("count : $target : save $tag.gz into count.$tsd_size");
		&monitorWait;
		$cmd = "perl $0 a=$a,sub=countZipFunc,ref=$ref,target=$target,tag=$tag,tsd_size=$tsd_size,sort_tmp=$sort_tmp &";
		system($cmd);
	    }
	}
    }
    &join;
}

sub countFunc{
    if ($target eq $ref){
	$cmd = "zcat $wd/$target/tmp/$tag.gz |sort -S 1M -T $wd/$a/tmp |uniq -c |awk '{print \$2 \"\t\" \$1 * 50}'|gzip > $wd/$target/tmp/$tag.count.gz && rm $wd/$target/tmp/$tag.gz";
    }else{
	$cmd = "zcat $wd/$target/tmp/$tag.gz |sort -S 1M -T $wd/$target/tmp |uniq -c |awk '{print \$2 \"\t\" \$1}' | gzip > $wd/$target/tmp/$tag.count.gz && rm $wd/$target/tmp/$tag.gz";
    }
    system($cmd);
}

sub countZipFunc{
    $cmd = "find $wd/$target/tmp/ -name \"$tag*.count.gz\" |sort -S 1M -T $wd/$target/tmp | xargs zcat | gzip > $wd/$target/count.$tsd_size/$tag.gz && find $wd/$target/tmp/ -name \"$tag*.count.gz\" | xargs rm";
    system($cmd);

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
	system($cmd);
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
    system($cmd);
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
    
    foreach $tag (@tag){
	$fout = "$tag.$chr";
	close($fout);
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
	    system($cmd);
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
    foreach $tag (@tag){
	if ($tag{$tag}){
	    &monitorWait;
	    $cmd = "perl $0 a=$a,b=$b,ref=$ref,sub=sort20mer,tag=$tag,sort_tmp=$sort_tmp &";
	    &log($cmd);
	    system($cmd);
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
        select undef, undef, undef, $wait_time if $sub ne "countFunc";
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

sub isRepeat{
    my $seq = shift;
    my ($i, $rep, $count);
    foreach $rep ('A', 'T'){
	$count = 0;
	for ($i = 0; $i < 20; $i++){
	    if ($rep eq substr($seq, $i, 1)){
		$count++;
	    }
	}
	if ($count >= 16){
	    return 1;
	}
    }
    foreach $rep ('AC', 'AG', 'AT', 'CG', 'CT', 'GT'){
	$count = 0;
	for ($i = 0; $i < 19; $i++){
	    if ($rep eq substr($seq, $i, 2)){
		$count++;
	    }
	}
	if ($count >= 8){
	    return 1;
	}
    }
    foreach $rep ('AAT', 'ATT', 'CAA', 'CCA', 'CCT', 'CTT', 'GAA', 'GGA', 'GGT', 'GTT'){
	$count = 0;
	for ($i = 0; $i < 18; $i++){
	    if ($rep eq substr($seq, $i, 3)){
		$count++;
	    }
	}
	if ($count >= 5){
	    return 1;
	}
    }
    return 0;
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
