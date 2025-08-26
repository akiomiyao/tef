#!/usr/bin/perl

$target = $ARGV[0];
$head   = $ARGV[1];
$tail   = $ARGV[2];
$llimit = $ARGV[3];
$ulimit = $ARGV[4];

$llimit = 200 if $llimit == 0;
$ulimit = 20000 if $ulimit == 0;

if ($ARGV[2] eq ""){
    print "
Usage: perl sliceTE.pl ref/fasta_file TE_head_sequence TE_tail_sequence [Lower_limit] [Upper_limit]

Compressed FASTA files (.gz, .bz2, or .xz) are also supported.

Size limits are optional. Default: Lower_limit = 200, Upper_limit = 20000

Example usage:
    perl sliceTE.pl TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz TGATAATAGGTTTTTCAGTC TGCAAAATACACACTTATCA 100 50000

";
    exit;
}

$s = {};

if ($target =~ /gz$/){
    open(IN, "zcat $target|");
}elsif($target =~ /bz2$/){
    open(IN, "bzcat $target|");
}elsif($target =~ /xz$/){
    open(IN, "xzcat $target|");
}else{
    open(IN, $target);
}
while(<IN>){
    chomp;
    if (/^>/){
	@row = split;
	($name = $row[0]) =~ s/^>//;
    }else{
	y/a-z/A-Z/;
	$chr{$name} .= $_;
    }
}

foreach $name (sort keys %chr){
    $seq = $chr{$name};
    while(1){
	$pos = index($seq, $head, $pos + 1);
	last if $pos == -1;
	$s->{$name}{$pos} = "head";
    }
    while(1){
	$pos = index($seq, $tail, $pos + 1);
	last if $pos == -1;
	$s->{$name}{$pos} = "tail";
    }

}

foreach $name (sort keys %$s){
    foreach $pos (sort bynumber keys %{$s->{$name}}){
	$hpos = $pos if ($s->{$name}{$pos} eq "head");
	if ($s->{$name}{$pos} eq "tail"){
	    $length = $pos - $hpos + length($tail);
	    $te = substr($chr{$name}, $hpos, $length);
	    $hpos ++;
	    $pos += length($tail);
#	    print ">$name $hpos-$pos $length bp\n$te\n" if $length < 20000;
#	    print ">$name $hpos-$pos $length bp\n$te\n" if $length > 4000 and $length < 4500;
	    print ">$name $hpos-$pos $length bp\n$te\n" if $length > $llimit and $length < $ulimit;
	}
    }
}

sub bynumber{
    $a <=> $b;
}
