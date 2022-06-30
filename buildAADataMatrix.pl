#!/usr/bin/perl -w

use strict;

use File::Basename;
use Carp;
use Getopt::Long;
use Data::Dumper;

use Common;
use PDB;

my $prog = basename ($0);
my $verbose = 0;


GetOptions ("v"=>\$verbose);


if (@ARGV != 4)
{
	print "build data matrix per amino acid\n";
	print "Usage: $prog [options] <pdb.list> <3dna.matrix.dir> <rbs.txt> <out.txt>\n";
	print " -v : verbose\n";
	exit (1);
}


my ($pdbListFile, $AAMatrixDir, $RBSFile, $outFile) = @ARGV;

my $fin;


print "reading RBS file $RBSFile ...\n" if $verbose;

open ($fin, "<$RBSFile") || Carp::croak "cannot openf ile $RBSFile to read\n";

my %RBSHash;

while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	next if $line =~/^\#/;

	my ($acc, $site, $score, $structStr) = split ("\t", $line);
	
	my @structs = split (/\|\|/, $structStr);
	#my @rbs = ($acc, $site, $score);
	#Carp::croak Dumper (\@structs), "\n";	

	foreach my $struct (@structs)
	{
		$RBSHash{$struct}{$site}{'acc'} = $acc;
		$RBSHash{$struct}{$site}{'score'} += $score;
	}
	
}
close ($fin);


#print Dumper ($RBSHash{'2ERR'}), "\n";

open ($fin, "<$pdbListFile") || Carp::croak "cannot open file $pdbListFile to read\n";

my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
my $printHeader = 1;

while (my $pdb = <$fin>)
{
	chomp $pdb;
	next if $pdb =~/^\s*$/;

	#next unless $pdb eq '2ERR';

	print "processing $pdb ...\n" if $verbose;

	my $AAMatrixFile = $AAMatrixDir . "/" . $pdb . ".out.aa.txt";
	my $fin2;
	
	open ($fin2, "<$AAMatrixFile") || Carp::croak "cannot open file $AAMatrixFile to read\n";

	my $header = <$fin2>;
	chomp $header;

	if ($printHeader)
	{
		print $fout join ("\t", "pdb", "protein_acc", "xl_score", $header), "\n";
		$printHeader = 0;
	}

	my $keep = 0;
	my @lines;

	my $acc = "";

	while (my $line = <$fin2>)
	{
		chomp $line;
		next if $line =~/^\s*$/;

		my ($idstr, $aacode, @cols) = split ("\t", $line);
		my ($chain, $site) = split (/\./, $idstr);
		

		my $pos = substr ($site, 3);
		
		$aacode = three2one ($aacode);
		next if $aacode eq '';

		$site = $aacode . $pos;
#		print "site=$site\n";

		if (exists $RBSHash{$pdb} && exists $RBSHash{$pdb}{$site})
		{
			print "$site is crosslinked\n" if $verbose;

			$keep = 1;
			$acc = $RBSHash{$pdb}{$site}{'acc'};
			push @lines, join ("\t", $RBSHash{$pdb}{$site}{'score'}, $line);
		}
		else
		{
			push @lines, join ("\t", 0, $line);
		}		
	}
	close ($fin2);	

	if ($keep)
	{
		foreach my $line (@lines)
		{
			print $fout join ("\t", $pdb, $acc, $line), "\n";
		}
	}
}


close ($fin);
close ($fout);



