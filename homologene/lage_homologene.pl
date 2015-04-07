#!/usr/bin/perl
#perl lage_homologene.pl homologene_aged.txt homologene.data.txt
#
local $|=1;
use strict;
use warnings;
use List::Util qw(max min);
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::NodeI;
use Bio::Tree::DistanceFactory;
use Bio::Annotation::Collection;
use XML::Bare;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Bio::Tools::Run::Phylo::Phyml;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Muscle;

my $fam_age_file = shift;
my $homologene_data = shift;
my %gi_geneid;
my %fam_gis;
my %human_gis;
my %len_ages;
my %rates;
my %geneid_famage;

open(my $fam_age_fh,"<",$fam_age_file) or die $!;
while(<$fam_age_fh>){
	next if $.==1;
	chomp;
	(my $gene_id, my $family_id, my $age) = (split("\t",$_))[0,1,6];
	$geneid_famage{$gene_id}=$age;
}
close $fam_age_fh;

open(my $data_fh,"<",$homologene_data) or die $!;
while(<$data_fh>){
	chomp;
	(my $famid, my $taxid, my $gene_id, my $prot_gi) = (split("\t",$_))[0,1,2,4];
	if (exists $fam_gis{$famid}){
		$fam_gis{$famid}.=';'.$prot_gi;
	}
	else{
		$fam_gis{$famid} = $prot_gi;
	}
	$gi_geneid{$prot_gi} = $gene_id;
	$human_gis{$prot_gi} = 1 if $taxid == 9606;
}
close $data_fh;

my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                       -db      => 'protein',
                                       -rettype => 'fasta',
                                       -email   => 'ajorr1@asu.edu');
my $counter = 1;
FAMILY: for my $fam (sort keys %fam_gis){
	$counter++;
	print "\rOn family: $counter of ".scalar(keys %fam_gis) if $counter%100==0;
	eval{{
		my @ids = split(';',$fam_gis{$fam});
		$factory->reset_parameters(-id => \@ids,
								   -eutil => 'efetch',
								   -db => 'protein',
								   -rettype => 'fasta',
								   -email => 'ajorr1@asu.edu');
		my $fasta = $factory->get_Response()->content();
		open(my $fastafh, "<", \$fasta);
		my $seqiofac = Bio::SeqIO->new(-fh => $fastafh,
									   -format => 'fasta',
									   -alphabet => 'protein');
		
		my @seqs;
		while(my $seq = $seqiofac->next_seq){
			push @seqs, $seq;
		}
		next if scalar(@seqs) lt 3;
		close $fastafh;
		my $musclefac = Bio::Tools::Run::Alignment::Muscle->new(QUIET=>1);
		my $alignment = $musclefac->align(\@seqs);
		###PHYML STUFF HERE
		my $alignfile = 'tempaln';
		my $treefile = $alignfile.'_phyml_tree.txt';
		my @phymlfiles = ($treefile, $alignfile.'phyml_stats.txt',
			$alignfile.'phyml_boot_trees.txt', $alignfile.'phyml_boot_stats.txt',
			$alignfile.'phyml_rand_trees.txt');
		my $alignout = Bio::AlignIO->new(-file => ">$alignfile", -format => 'phylip');
		$alignout->write_aln($alignment);
		
		system("phyml -i $alignfile -d aa --no_memory_check >nul")==0 or die;
		
		my $treefac = Bio::TreeIO->new(-file => $treefile, -format => 'newick');
		my $tree = $treefac->next_tree;
		###PHYML STUFF ENDS
		
		
		
		
		
		###MEGA STUFF HERE
		# my $alignfile = 'tempaln.fasta';
		# my $treefile = $alignfile.'_phyml_tree.txt';
		# my @phymlfiles = ($treefile, $alignfile.'phyml_stats.txt',
			# $alignfile.'phyml_boot_trees.txt', $alignfile.'phyml_boot_stats.txt',
			# $alignfile.'phyml_rand_trees.txt');
		# my $alignout = Bio::AlignIO->new(-file => ">$alignfile", -format => 'fasta');
		# $alignout->write_aln($alignment);
		# `M6CC -a M6CC.mao -d tempaln.fasta -o temptree -s`;
		# my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $treefile);
		# my $tree = $treeio->next_tree;
		
		###MEGA STUFF ENDS
		my @leaves = $tree->get_leaf_nodes;
		my @human_leaves;
		my $root_node = $tree->get_root_node;
		my $root_length = $root_node->height();
		next if $root_length==0;
		my @lengths;
		
		for my $leaf (@leaves){
			my $id = $leaf->id;
			my $depth = $leaf->depth;
			(my $field,$id) = (split(/\|/,$id))[0,1];
			warn "\n$id is not a gi, it's a $field\n" if $field ne 'gi';
			push @lengths, $depth;
			push @human_leaves, $leaf if exists $human_gis{$id};
		}
		my $min_length = min(@lengths);
		for my $leaf (@human_leaves){
			my $id = $leaf -> id();
			$id = (split(/\|/,$id))[1];	
			next unless(exists $gi_geneid{$id});
			my $geneid = $gi_geneid{$id};
			next unless(exists $geneid_famage{$geneid});
			my $fam_age = $geneid_famage{$geneid};
			my $parent_node = $leaf->ancestor;
			my $depth = $parent_node->depth;
			$len_ages{$geneid} = (1 - ($depth/$root_length))*($fam_age);
			$rates{$geneid} = ($root_length-$min_length)/($fam_age);
			#print(join("\t",$id,$fam_age,$len_ages{$uniprot_id}),"\n");
			# print "uni_id:$uniprot_id\n";
			# print "family_age:$fam_age\n";
			# print "len_age:$len_ages{$uniprot_id}";
			# exit;			
		}
	}};
	if ($@){
		warn "\nI had a problem with fam:$fam"."\n".$@;
		}
}
open($fam_age_fh,"<",$fam_age_file) or die $!;
open(my $newfh,">>","rates_len_".$fam_age_file) or die $!;
print{$newfh}(join("\t",qw(gene_id refseq fam_age len_age rate)),"\n");

while(<$fam_age_fh>){
	next if $.==1;
	chomp;
	(my $id, my $refseq, my $age) = (split("\t",$_))[0,3,6];
	$len_ages{$id} = "NA" unless exists $len_ages{$id};
	$rates{$id} = "NA" unless exists $rates{$id};
	print {$newfh} (join("\t",$id,$refseq,$age,$len_ages{$id},$rates{$id}),"\n");
	}
close $fam_age_fh;
close $newfh;
exit;
