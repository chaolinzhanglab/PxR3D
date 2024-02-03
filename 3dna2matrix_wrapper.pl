#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Carp;

use Data::Dumper;
use Getopt::Long;

use Bio::Structure::IO;

use MyConfig;
use PDB;

my $prog = basename ($0);
my $progDir = dirname ($0);

my $verbose = 0;
#my $modelNum = 1;
my $keepCache = 0;

my $CMD_3dna = "~/tools/3dna/";
my $cache = MyConfig::getDefaultCache ($prog);

GetOptions (
#	"nmr:i"=>\$modelNum,
	"c:s"=>\$cache,
	"3dna:s"=>\$CMD_3dna,
	"keep-cache"=>\$keepCache,
	"v"=>\$verbose);

if (@ARGV != 3)
{
	print "$prog [options] <in.pdb> <out.rna.txt> <out.aa.txt>\n";
#	print " --nmr [int]   : specify the number of models if nmr\n";
	print " --3dna [stirng]: path to 3DNA package\n";
	print " -c     [string]: cache dir ($cache)\n";
	print " --keep-cache   : keep cache dir\n";
	print " -v             : verbose\n";
	exit (1);
}

my ($pdbFile, $outRNAFile, $outAAFile) = @ARGV;

my $CMD_dssr = "$CMD_3dna/x3dna-dssr.dms";
my $CMD_snap = "$CMD_3dna/x3dna-snap.dms";

my $ret = system ("mkdir $cache");
Carp::croak "cannot create cache dir $cache\n" if $ret != 0;


my $io = Bio::Structure::IO->new(-file => $pdbFile, -format => 'PDB');

my $struct = $io->next_structure;

my $modelNum = $struct->get_models;
print "Total model number = $modelNum\n" if $verbose;


if ($modelNum <= 1)
{
	my $dssr_prefix = "$cache/pdb_dssr.out";
	my $cmd = "$CMD_dssr -i=$pdbFile --part=nts --summary --prefix=$dssr_prefix -o=$dssr_prefix";

	print "CMD=$cmd\n" if $verbose;

	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	
	
	my $snap_prefix = "$cache/pdb_snap.out";
	$cmd = "$CMD_snap -i=$pdbFile --interface -o=$snap_prefix";
	
	print "CMD=$cmd\n" if $verbose;
	
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;



	$cmd = "perl $progDir/3dna2matrix.pl -v $pdbFile $dssr_prefix-summary.csv $snap_prefix $outRNAFile $outAAFile";

	print "CMD=$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
}
else
{

	my $ntHeaderLine = "";
	my $aaHeaderLine = "";
	my (@ntidstr, @ntcode, @ntchain, @ntdata, @aaidstr, @aacode, @aachain, @aadata);

	for (my $i = 1; $i<= $modelNum; $i++)
	{
		print "\nModel=$i ...\n" if $verbose;
		my $modelFile = "$cache/model.$i.pdb";
		my $cmd = "$CMD_dssr -i=$pdbFile --select-model=$i -o=$modelFile";
		$ret = system ($cmd);
        Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
		

		my $dssr_prefix = "$cache/pdb_dssr.$i.out";
		$cmd = "$CMD_dssr -i=$modelFile --part=nts --summary --prefix=$dssr_prefix -o=$dssr_prefix";

		print "CMD=$cmd\n" if $verbose;

		$ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;
	
	
		my $snap_prefix = "$cache/pdb_snap.$i.out";
		$cmd = "$CMD_snap -i=$modelFile --interface -o=$snap_prefix";
	
		print "CMD=$cmd\n" if $verbose;
	
		$ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;


		my $tmpOutRNAFile = "$cache/out.$i.rna.txt";
		my $tmpOutAAFile = "$cache/out.$i.aa.txt";
		$cmd = "perl $progDir/3dna2matrix.pl -v $pdbFile $dssr_prefix-summary.csv $snap_prefix $tmpOutRNAFile $tmpOutAAFile";

		print "CMD=$cmd\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

		print "loading output of model $i...\n" if $verbose;

		my $fin;
		open ($fin, "<$tmpOutRNAFile") || Carp::croak "cannot open file $tmpOutRNAFile to read\n";
		$ntHeaderLine = <$fin>;
		chomp $ntHeaderLine;

		my $ntNumModel = 0;
		while (my $line = <$fin>)
		{
			chomp $line;
			next if $line =~/^\s*$/;
			
			my ($idstr, $ntcode, $chain, @cols) = split ("\t", $line);
			if ($i == 1)
			{
				push @ntidstr, $idstr;
				push @ntcode, $ntcode;
				push @ntchain, $chain;
			}

			for (my $j = 0; $j < @cols; $j++)
			{
				$ntdata[$ntNumModel][$j] += $cols[$j];
			}
			$ntNumModel++;
		}
		close ($fin);
		my $ntNum = @ntidstr;

		Carp::croak "inconsistency of ntNumModel = $ntNumModel in model $i ($ntNum expected)\n" unless $ntNumModel == $ntNum;


		open ($fin, "<$tmpOutAAFile") || Carp::croak "cannot open file $tmpOutAAFile to read\n";
		$aaHeaderLine = <$fin>;
		chomp $aaHeaderLine;

		my $aaNumModel = 0;
		while (my $line = <$fin>)
		{
			chomp $line;
			next if $line =~/^\s*$/;
			
			my ($idstr, $aacode, $chain, @cols) = split ("\t", $line);
			if ($i == 1)
			{
				push @aaidstr, $idstr;
				push @aacode, $aacode;
				push @aachain, $chain;
			}

			for (my $j = 0; $j < @cols; $j++)
			{
				$aadata[$aaNumModel][$j] += $cols[$j];
			}
			$aaNumModel++;
		}
		close ($fin);
		my $aaNum = @aaidstr;

		Carp::croak "inconsistency of ntNumModel = $aaNumModel in model $i ($aaNum expected)\n" unless $aaNumModel == $aaNum;

	} 

	print "write output to $outRNAFile ...\n" if $verbose;
	my $fout;
	open ($fout, ">$outRNAFile") || Carp::croak "cannot open file $outRNAFile to write\n";
	print $fout $ntHeaderLine, "\n";
	for (my $i = 0; $i < @ntidstr; $i++)
	{
		print $fout join ("\t", $ntidstr[$i], $ntcode[$i], $ntchain[$i]);
		my @row = map {$_/$modelNum} @{$ntdata[$i] };
		print $fout "\t", join ("\t", @row), "\n";
	}
	close ($fout);
	
	print "write output to $outAAFile ...\n" if $verbose;
	open ($fout, ">$outAAFile") || Carp::croak "cannot open file $outAAFile to write\n";
	print $fout $aaHeaderLine, "\n";
	for (my $i = 0; $i < @aaidstr; $i++)
	{
		print $fout join ("\t", $aaidstr[$i], $aacode[$i], $aachain[$i]);
		my @row = map {$_/$modelNum} @{$aadata[$i] };
		print $fout "\t", join ("\t", @row), "\n";
	}
	close ($fout);
}

system ("rm -rf $cache") unless $keepCache;




