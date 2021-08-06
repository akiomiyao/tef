#!/usr/bin/perl

open(IN, $ARGV[0]);
binmode(IN);
seek(IN, $ARGV[1]-1, 0);
read(IN, $seq, $ARGV[2] - $ARGV[1] + 1);
print "$seq\n";
