#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Basename;
use Carp;
use Data::Dumper;

use Bio::Structure::IO;
#use Bio::Annotation::SimpleValue;

use PDB;

my $prog = basename ($0);
my $verbose = 0;

my $summaryFile = "";
my $dirListFile = "";

GetOptions ("list:s"=>\$dirListFile,
		"summary:s"=>\$summaryFile,
	"v"=>\$verbose);

if (@ARGV != 2)
{
	print "extract (binary) complex from pdb structures\n";
	print "$prog [options] <in.dir> <out.dir>\n";
	print " -list    [string]: consider only subdirs in the list\n";
	print " -summary [string]: summary file\n";
	print " -v : verbose\n";
	exit (1);
}

my ($inDir, $outDir) = @ARGV;


my %dirHash;

my $fin;
if (-f $dirListFile)
{
	print "reading list of subdirs to process...\n" if $verbose;
	open ($fin, "<$dirListFile") || Carp::croak "cannot open file $dirListFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		$dirHash{$line} = 1;
	}
	close($fin);

	my $n = keys %dirHash;
	print "$n subdirs to be processed\n" if $verbose;
	Carp::croak "no subdirs to be processed?\n" unless $n > 0;
}

my $dh;
opendir ($dh, $inDir) || Carp::croak "cannot open $inDir to read\n";

my $fout;
if ($summaryFile ne '')
{
	open ($fout, ">$summaryFile") || Carp::croak "cannot open file $summaryFile to write\n";
}

mkdir ($outDir);

print $fout join ("\t", "#struct_id", "title", 
						"mol1_id", "mol1_chain", "mol1_dbref", "mol1_acc", "mol1_name", "mol1_desc", "mol1_type", "mol1_seq",
						"mol2_id", "mol2_chain", "mol2_dbref", "mol2_acc", "mol2_name", "mol2_desc", "mol2_type", "mol2_seq"), "\n";

while (my $d = readdir ($dh))
{
	next if $d eq '.' || $d eq '..';

	my $dir = "$inDir/$d";
	next unless -d $dir;

	if (-f $dirListFile)
	{
		next unless exists $dirHash{$d};
	}

	my $dh2;
	opendir ($dh2, $dir) || Carp::croak "cannot open dir $dir to read\n";
	while (my $f = readdir ($dh2))
	{
		next if $f eq '.' || $f eq '..';

		my $inFile = "$outDir/$f";
		$inFile =~s/\.gz$//;
		my $cmd = "zcat $dir/$f | grep -v \"^HETATM\" | grep -v \"ANISOU\" | grep -v \"SIGATM\" | grep -v \"SIGUIJ\"> $inFile";
		#removed teh HETATM, since these were not parsed by bioperl properly
		my $ret = system ($cmd);
		Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;

		print "processing $inFile ...\n" if $verbose;
		
		#check if there are still ATOM lines
		my $n = `grep \"^ATOM\" $inFile | wc -l`;
		chomp $n;
		if ($n == 0)
		{
			unlink $inFile;
			next;
		}

		my $io = Bio::Structure::IO->new(-file => $inFile, -format => 'PDB');
		my $struct = $io->next_structure;
		
		my $annot = $struct->annotation;
		my ($title) = $annot->get_Annotations('title');
		my $titleStr = $title->value;
    	$titleStr=~s/\s+$//g;
		
		my $complexInfo = extractComplexInfo ($struct);

		#system ("rm -rf $dir/$f"); #remove the file so that we do not scan next time
		if ($complexInfo == 0)
		{
			unlink $inFile;
			next;
		}

		my $complex = $complexInfo->{'chains'};
		my $molecules = $complexInfo->{'molecules'};

		#Carp::croak Dumper($molecules), "\n";		

		if (@$molecules != 2)
		{
			unlink $inFile;
			next;
		}
		
		my @chainid = sort keys %$complex;

		if ($summaryFile ne '')
		{
			print $fout join ("\t", $struct->id, $titleStr,
				$molecules->[0]{'id'}, $molecules->[0]{'chain'}, $molecules->[0]->{'db'}, $molecules->[0]->{'acc'}, $molecules->[0]->{'name'}, $molecules->[0]{'desc'}, $molecules->[0]->{'type'}, $molecules->[0]->{'seq'}, 
				$molecules->[1]{'id'}, $molecules->[1]{'chain'}, $molecules->[1]->{'db'}, $molecules->[1]->{'acc'}, $molecules->[1]->{'name'}, $molecules->[1]{'desc'}, $molecules->[1]->{'type'}, $molecules->[1]->{'seq'}), "\n";
		}
	}
	closedir ($dh2);
}

closedir ($dh);

if ($summaryFile ne '')
{
	close ($fout);
}


