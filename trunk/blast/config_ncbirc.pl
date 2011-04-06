#!/usr/local/bin/perl

use Shell qw(uname);
$arch = uname("-m");
chomp $arch;
$gatroot = $ENV{GATHOME};
if(!$gatroot) { die "GATHOME environment variable not set.\n"; }

$ncf = "$gatroot/tools/blastall/.ncbirc";
open NCBI, ">$ncf" or die "Couldn't open $ncf.\n";
print NCBI "[NCBI]\n";
print NCBI "ROOT = $gatroot/tools/blastall\n";
print NCBI "DATA = $gatroot/tools/blastall/data\n";

