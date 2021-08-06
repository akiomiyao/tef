#!/usr/bin/perl

$target = $ARGV[0];
$head   = $ARGV[1];
$tail   = $ARGV[2];

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
	    print ">$name $hpos-$pos $length bp\n$te\n" if $length > 4000 and $length < 4500;
	}
    }
}

sub bynumber{
    $a <=> $b;
}
