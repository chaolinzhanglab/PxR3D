#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use File::Basename;
use Carp;
use Data::Dumper;

use Bio::Structure::IO;

use PDB;

my $prog = basename ($0);
my $verbose = 0;

my @aaTable = qw(ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SEC SER THR TRP TYR VAL ASX GLX);


GetOptions ("v"=>\$verbose);

if (@ARGV != 2)
{
	print "extract rna-protein interaface from pdb structure\n";
	print "$prog [options] <in.pdb> <out.txt>\n";
	print " -v : verbose\n";
	exit (1);
}

my ($inFile, $outFile) = @ARGV;

print "loading structure from $inFile...\n" if $verbose;

my $io = Bio::Structure::IO->new(-file => $inFile, -format => 'PDB');
my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

print $fout join ("\t", "#struct_id", "model_id", "chain_id", "nucl_id", "closest_aa", "closest_dist", "nucl_atom", "aa_atom", @aaTable), "\n";
while (my $struct = $io->next_structure) 
{

	my $structId = $struct->id;

	print "process structure $structId...\n" if $verbose;

	#Carp::croak "structid=$structId\n";	

	my $complexInfo = extractComplexInfo ($struct);
	my $molecules = $complexInfo->{'molecules'};

	Carp::croak "not a binary complex:", Dumper ($molecules), "\n"  unless @$molecules == 2;
	
	Carp::croak "not a protein-rna complex: ", Dumper ($molecules), "\n"
        unless ($molecules->[0]{'type'} eq 'protein' && $molecules->[1]{'type'} eq 'rna')
            || ($molecules->[0]{'type'} eq 'rna' && $molecules->[1]{'type'} eq 'protein');

	#get the list of chains for each molecule
	foreach my $mol (@$molecules)
	{
		my $chainStr = $mol->{'chain'};
		my @chains = split(/\s*\,\s*/, $mol->{'chain'});
		my %chainHash = map {$_=>1} @chains;
		$mol->{'chainHash'} = \%chainHash;
	}

	#Carp::croak "chainHash=", Dumper ($molecules), "\n";

	my ($protein, $rna) = @$molecules;
	
	if ($protein->{'type'} eq 'rna')
	{
		($protein, $rna) = ($rna, $protein);
	}

	my $iter = 0;
	foreach my $model ($struct->get_models)
	{
		my $modelId = $model->id;
		print "pocessing model $modelId...\n" if $verbose;
		
		#get the list of chains for protein and rna	
		my (@proteinChains, @rnaChains);
		foreach my $chain ($struct->get_chains ($model))
		{
			my $chainId = $chain->id;
			if (exists $protein->{'chainHash'}{$chainId})
			{
				push @proteinChains, $chain;
			}
			elsif (exists $rna->{'chainHash'}{$chainId})
			{
				push @rnaChains, $chain;
			}
			else
			{
				Carp::croak "unexpected chain $chainId\n";
			}
		}	

		foreach my $rnaChain (@rnaChains)
		{
			my $rnaChainId = $rnaChain->id;
			foreach my $nucl ($struct->get_residues ($rnaChain))
			{
				my $nuclId = $nucl->id;
			
				#Carp::croak "nuclId = $nuclId\n";

				my %distHash;
				
				foreach my $proteinChain (@proteinChains)
				{

					foreach my $aa ($struct->get_residues ($proteinChain))
					{
						my $aaId = $aa->id;
						#Carp::croak "aaId = $aaId\n";
	
						my ($res, $n) = split("-", $aaId);

						my $d2 = res2res ($struct, $nucl, $aa);
						#Carp::croak Dumper ($d2), "\n";

						if (exists $distHash{$res})
						{
							if ($d2->{'d'} < $distHash{$res}{'d'})
							{
								$d2->{'id'} = $aaId;
								$distHash{$res} = $d2;
							}
						}
						else
						{
							$d2->{'id'} = $aaId;
							$distHash{$res} = $d2;
						}
					}
				}
				
				my @dist = map {exists $distHash{$_} ? $distHash{$_}{'d'} : 'NA'} @aaTable;
				#Carp::croak Dumper (\@dist), "\n";			

				#get the closest amino acid
				my @aaList = values %distHash;
				my $partner = $aaList[0];
			
				foreach my $aa (@aaList)
				{
					$partner = $aa if $partner->{'d'} > $aa->{'d'};
				}
			
				print $fout join ("\t", $structId, $modelId, $rnaChainId, $nuclId, $partner->{'id'}, $partner->{'d'}, $partner->{'atom1'}, $partner->{'atom2'}, @dist), "\n";
			} #for $nucl
		} #for $rnaChain
	} #for model
}

close ($fout);

