#!/usr/bin/env perl

use strict;

use File::Basename;
use Carp;
use Getopt::Long;
use Data::Dumper;

use Bio::Structure::IO;

use Common;
use PDB;

my $prog = basename ($0);
my $verbose = 0;

my @aaTable = qw (ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL);
my @ntTable = qw (A C G U);

my %aaHash = map {$aaTable[$_] => $_} (0..($#aaTable));

my %aaGroups = (
	polar=>[qw(GLN ASN HIS SER THR TYR CYS TRP)],
	positive=>[qw(LYS ARG HIS)],
	negative=>[qw(ASP GLU)],
	hydrophobic=>[qw(ALA ILE LEU MET PHE VAL PRO GLY)],
	aromatic=>[qw(PHE TYR TRP HIS)],
	aliphatic=>[qw(ALA LEU VAL PRO ILE)]
);

my @aaGroupNames = qw(polar positive negative hydrophobic aromatic aliphatic);

my @hbondType = qw(po4:sidechain po4:backbone sugar:sidechain sugar:backbone base:sidechain base:backbone);
my %hbondTypeHash = map {$_=>1} @hbondType;


GetOptions ("v"=>\$verbose);


if (@ARGV != 5)
{
	print "parse 3dna output to matrix format\n";
	print "Usage: $prog [options] <in.pdb> <dssr.summary.csv> <snap.out> <out.rna.mat> <out.aa.mat>\n";
	print " -v : verbose\n";
	exit (1);
}


my ($pdbFile, $dssrFile, $snapFile, $outRNAFile, $outAAFile) = @ARGV;

my $fin;

open ($fin, "<$dssrFile") || Carp::croak "cannot openf ile $dssrFile to read\n";

#my @dssrKey = qw(modified pseudoknotted turn break u-turn k-turn sugar_conf puckering BI
#	non-stack canonical non-canonical non-pair-contact 
#	helix index stem index coaxial-stack multiplet index hairpin-loop index bulge index internal-loop index junction-loop index 
#	ss-non-loop A-minor G-tetrad i-motif ribose-zipper kissing-loop cap-donor cap-acceptor phosphate splayed-apart);


my %rnaDataHash;
my $headerLine = <$fin>;

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	my @cols = split (',', $line);

	my ($idstr, $ntcode, $chain, $singleStranded,  $sugarConf, $puckering, $BI, $ntstacking) = ($cols[0], $cols[1], $cols[3], $cols[5], $cols[12], $cols[13], $cols[14], $cols[15]);	
	#$nt = substr($nt, 2);
	
	my $sugarConfCoded = {anti=>0, syn=>0, other=>0};
	$sugarConf = "other" if $sugarConf eq '.';
	
	if (exists $sugarConfCoded->{$sugarConf})
	{
		$sugarConfCoded->{$sugarConf} = 1;
	}
	else
	{
		Carp::croak "sugarConf=$sugarConf is undefined\n";
	}

	my $puckeringCoded = {"~C3'-endo"=>0, "~C2'-endo"=>0, other=>0};
	$puckering = "other" if $puckering eq '.';

	if (exists $puckeringCoded->{$puckering})
	{
		$puckeringCoded->{$puckering} = 1;
	}
	else
	{
		Carp::croak "puckering=$puckering is undefined\n";
	}

	my $BICoded = {BI=>0, BII=>0, other=>0};
	$BI = "other" if $BI eq ".";

	if (exists $BICoded->{$BI})
	{
		$BICoded->{$BI} = 1;
	}
	else
	{
		Carp::croak "BI=$BI is undefined\n";
	}

	my $ntstackingCoded = {upper=>0, down=>0};
	if ($ntstacking eq 'upper:')
	{
		$ntstackingCoded->{'upper'} = 1;
	}
	elsif ($ntstacking eq ':down')
	{
		$ntstackingCoded->{'down'} = 1;
	}
	elsif ($ntstacking eq 'upper:down')
	{
		$ntstackingCoded->{'upper'} = 1;
		$ntstackingCoded->{'down'} = 1;
	}
	
	$rnaDataHash{$idstr} = {
			iter=>$iter,
			nt=>$ntcode,
			chain=>$chain,
			single_stranded=> $singleStranded eq 'single-stranded' ? 1 : 0,
			sugar_conf=>$sugarConfCoded,
			puckering=>$puckeringCoded,
			BI=>$BICoded,
			ntstacking=>$ntstackingCoded
						
	};
	$iter++;	
}

close ($fin);

#Carp::croak Dumper (\%dataHash), "\n";

my %aaDataHash;
open ($fin, "<$snapFile") || Carp::croak "cannot open file $snapFile to read\n";

while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;

	#print $line, "\n";
	if ($line=~/^List of (\d+) nucleotide\/amino-acid interactions/g)
	{

		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntaa, $res_nt, $res_aa, @junc) = split (/\s+/, $line);
			

			#$nt=substr($nt, 2);
			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);
			
			
			#Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$rnaDataHash{$res_nt}{'interact'}{$aa}+= 1;
			$aaDataHash{$res_aa}{'interact'}{$nt} += 1;
		}
	}
	elsif ($line=~/^List of (\d+) phosphate\/amino-acid H-bonds/g)
	{
		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntatom, $aaatom, $dist, $type) = split (/\s+/, $line);
			
			my ($junc, $res_nt) = split("\@", $ntatom);

			my ($junc2, $res_aa) = split("\@", $aaatom);

			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);

			Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$type =~s/:salt-bridge$//g;
			#salt-bridge is treated as h-bond here
			Carp::croak "hbondType = $type is undefined\n" unless exists $hbondTypeHash{$type};
			
			$rnaDataHash{$res_nt}{'hbond'}{$type}{$aa}+= 1;
			$aaDataHash{$res_aa}{'hbond'}{$type}{$nt} += 1;

		}
		#Carp::croak Dumper ($dataHash{'hbond'});
	}
	elsif ($line=~/^List of (\d+) sugar\/amino-acid H-bonds/g)
	{
		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntatom, $aaatom, $dist, $type) = split (/\s+/, $line);
			
			my ($junc, $res_nt) = split("\@", $ntatom);

			my ($junc2, $res_aa) = split("\@", $aaatom);
			
			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);

			Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$type =~s/:salt-bridge$//g;
			#salt-bridge is treated as h-bond here
			Carp::croak "hbondType = $type is undefined\n" unless exists $hbondTypeHash{$type};
		
			$rnaDataHash{$res_nt}{'hbond'}{$type}{$aa}+= 1;
			$aaDataHash{$res_aa}{'hbond'}{$type}{$nt}+= 1;
		}
		#Carp::croak Dumper ($dataHash{'hbond'});
	}
	elsif ($line=~/^List of (\d+) base\/amino-acid H-bonds/g)
	{
		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntatom, $aaatom, $dist, $type) = split (/\s+/, $line);
			
			my ($junc, $res_nt) = split("\@", $ntatom);

			my ($junc2, $res_aa) = split("\@", $aaatom);
			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);

			#Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$type =~s/:salt-bridge$//g;
			#salt-bridge is treated as h-bond here
			Carp::croak "hbondType = $type is undefined\n" unless exists $hbondTypeHash{$type};

			$rnaDataHash{$res_nt}{'hbond'}{$type}{$aa}+= 1;
			$aaDataHash{$res_aa}{'hbond'}{$type}{$nt}+= 1;
		}
		#Carp::croak Dumper ($dataHash{'hbond'});
	}
	elsif ($line=~/^List of (\d+) base\/amino-acid pairs/g)
	{
		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntaa, $res_nt, $res_aa, $dist, $angle) = split (/\s+/, $line);
			
			#$nt = substr($nt, 2);
			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);

			Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$rnaDataHash{$res_nt}{'base_aa_pair'}{$aa}+= 1;
			$aaDataHash{$res_aa}{'base_aa_pair'}{$nt} += 1;
		}
		#Carp::croak Dumper ($dataHash{'base_aa_pair'});
	}
	elsif ($line=~/^List of (\d+) base\/amino-acid stacks/g)
	{
		my $head = <$fin>;
		while (my $line=<$fin>)
		{
			chomp $line;
			last if $line=~/^\s*$/;

			$line=~s/^\s+//g;

			my ($iter, $id, $ntaa, $res_nt, $res_aa, $dist, $angle) = split (/\s+/, $line);
			
			#$nt = substr($nt, 2);
			my $aa = cleanAA (substr($res_aa, 2, 3));
			my $nt = substr($res_nt, 2, 1);

			Carp::croak "aa=$aa is undefined\n" unless exists $aaHash{$aa};

			$rnaDataHash{$res_nt}{'base_aa_stack'}{$aa}+= 1;
			$aaDataHash{$res_aa}{'base_aa_stack'}{$nt}+= 1;
		}
		#Carp::croak Dumper ($dataHash{'base_aa_stack'});
	}
}

close ($fin);

#Carp::croak Dumper (\%dataHash), "\n";


my $fout;
open ($fout, ">$outRNAFile") || Carp::croak "cannot open file $outRNAFile to write\n";

print $fout join ("\t", "idstr", 'ntcode', 'chain', 'single_stranded', 'sugar_conf_anti', 'sugar_conf_syn', 'sugar_conf_other', 'puckering_C3_endo', 'puckering_C2_endo', 'puckering_other', 'BI_BI', "BI_BII", "BI_other", "ntstacking_upper", "ntstacking_down");

my @interact = map {"interact_" . $_} (@aaTable, @aaGroupNames);
print $fout "\t", join ("\t", @interact);

for my $type (@hbondType)
{
	my @hbond = map {"hbond_" . $type . "_" . $_} (@aaTable, @aaGroupNames);
	print $fout "\t", join ("\t", @hbond);
}

my @basepair = map {"base_aa_pair_" . $_} (@aaTable, @aaGroupNames);
print $fout "\t", join ("\t", @basepair); 

my @basestack = map {"base_aa_stack_" . $_} (@aaTable, @aaGroupNames);
print $fout "\t", join ("\t", @basestack); 
print $fout "\n";


foreach my $idstr (sort {$rnaDataHash{$a}->{'iter'} <=> $rnaDataHash{$b}->{'iter'}} keys %rnaDataHash)
{
	my $ntinfo = $rnaDataHash{$idstr};

	print $fout join ("\t", $idstr, 
					$ntinfo->{'nt'},
					$ntinfo->{'chain'},
					$ntinfo->{'single_stranded'},
					$ntinfo->{'sugar_conf'}{'anti'}, $ntinfo->{'sugar_conf'}{'syn'}, $ntinfo->{'sugar_conf'}{'other'},
					$ntinfo->{'puckering'}{"~C3'-endo"}, $ntinfo->{'puckering'}{"~C2'-endo"}, $ntinfo->{'puckering'}{"other"},
					$ntinfo->{'BI'}{'BI'}, $ntinfo->{'BI'}{'BII'}, $ntinfo->{'BI'}{'other'},
					$ntinfo->{'ntstacking'}{'upper'}, $ntinfo->{'ntstacking'}{'down'});

	my @interact = map {exists $ntinfo->{'interact'} && exists $ntinfo->{'interact'}{$_} ? $ntinfo->{'interact'}{$_} : 0} @aaTable;

	my @interact2 = aggregateAAGroups (\@interact, \@aaGroupNames);

	print $fout "\t", join ("\t", @interact, @interact2);

	foreach my $type (@hbondType)
	{
		my @hbond = map {exists $ntinfo->{'hbond'} && exists $ntinfo->{'hbond'}{$type} && exists $ntinfo->{'hbond'}{$type}{$_} ? $ntinfo->{'hbond'}{$type}{$_} : 0} @aaTable;
		my @hbond2 = aggregateAAGroups (\@hbond, \@aaGroupNames);

		print $fout "\t", join ("\t", @hbond, @hbond2);
	}
	
	my @basepair = map {exists $ntinfo->{'base_aa_pair'} && exists $ntinfo->{'base_aa_pair'}{$_} ? $ntinfo->{'base_aa_pair'}{$_} : 0} @aaTable;
	my @basepair2 = aggregateAAGroups (\@basepair, \@aaGroupNames);

	print $fout "\t", join ("\t", @basepair, @basepair2);

	my @basestack = map {exists $ntinfo->{'base_aa_stack'} && exists $ntinfo->{'base_aa_stack'}{$_} ? $ntinfo->{'base_aa_stack'}{$_} : 0} @aaTable;
    my @basestack2 = aggregateAAGroups (\@basestack, \@aaGroupNames);;

	print $fout "\t", join ("\t", @basestack, @basestack2);

	print $fout "\n";
}

close ($fout);


my $io = Bio::Structure::IO->new(-file => $pdbFile, -format => 'PDB');

my $struct = $io->next_structure;

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

my $model;

for $model ($struct->get_models)
{
	last;
}

my @proteinChains;
foreach my $chain ($struct->get_chains ($model))
{
	my $chainId = $chain->id;
    if (exists $protein->{'chainHash'}{$chainId})
    {
    	push @proteinChains, $chain;
	}
}


open ($fout, ">$outAAFile") || Carp::croak "cannot open file $outAAFile to write\n";

print $fout join ("\t", "idstr", 'aacode', 'chain');

@interact = map {"interact_" . $_} @ntTable;
print $fout "\t", join ("\t", @interact);

for my $type (@hbondType)
{
	my @hbond = map {"hbond_" . $type . "_" . $_} @ntTable;
	print $fout "\t", join ("\t", @hbond);
}

@basepair = map {"base_aa_pair_" . $_} @ntTable;
print $fout "\t", join ("\t", @basepair); 

@basestack = map {"base_aa_stack_" . $_} @ntTable;
print $fout "\t", join ("\t", @basestack); 
print $fout "\n";


foreach my $proteinChain (@proteinChains)
{
	my $proteinChainId = $proteinChain->id;
	foreach my $aa ($struct->get_residues ($proteinChain))
	{
		my $aaId = $aa->id;
		
		my ($res, $n) = split ("-", $aaId);
		next if $res eq 'HOH';

		my $idstr = $proteinChainId . "." .$res . $n;
	
		my $aainfo = $aaDataHash{$idstr};

		print $fout join ("\t", $idstr, 
					$res,
					$proteinChainId);

		my @interact = map {exists $aaDataHash{$idstr} && exists $aainfo->{'interact'} && exists $aainfo->{'interact'}{$_} ? $aainfo->{'interact'}{$_} : 0} @ntTable;

		print $fout "\t", join ("\t", @interact);

		foreach my $type (@hbondType)
		{
			my @hbond = map {exists $aaDataHash{$idstr} && exists $aainfo->{'hbond'} && exists $aainfo->{'hbond'}{$type} && exists $aainfo->{'hbond'}{$type}{$_} ? $aainfo->{'hbond'}{$type}{$_} : 0} @ntTable;

			print $fout "\t", join ("\t", @hbond);
		}
	
		my @basepair = map {exists $aaDataHash{$idstr} && exists $aainfo->{'base_aa_pair'} && exists $aainfo->{'base_aa_pair'}{$_} ? $aainfo->{'base_aa_pair'}{$_} : 0} @ntTable;

		print $fout "\t", join ("\t", @basepair);

		my @basestack = map {exists $aaDataHash{$idstr} && exists $aainfo->{'base_aa_stack'} && exists $aainfo->{'base_aa_stack'}{$_} ? $aainfo->{'base_aa_stack'}{$_} : 0} @ntTable;

		print $fout "\t", join ("\t", @basestack);

		print $fout "\n";
	}
}
close ($fout);




sub cleanAA
{
	my ($aa) = @_;
	if ($aa eq 'MSE')
	{
		$aa = 'MET';
	}
	return $aa;
}

sub aggregateAAGroups
{
	my ($val, $grps) = @_;
	my @val2;
    foreach my $grp (@$grps)
    {
        my $aaList = $aaGroups{$grp};
        my @aaIdx = map {$aaHash{$_}} @$aaList;
        my @tmp = map {$val->[$_]} @aaIdx;
        push @val2,  sum(\@tmp);
    }
	return @val2;
}


