#!usr/bin/perl
#perl lage_hogenom.pl hogenom_aged.txt acc_id5_id6
#

use strict;
use warnings;
use List::Util qw(max min);
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::NodeI;
use Bio::Annotation::Collection;
use XML::Bare;

my $agesfile = shift;
my $uni_hog_mapfile = shift;
open(my $agesfh,"<",$agesfile) or die $!;
my %fam_ages;
my %len_ages;
my %rates;
my %prefix_txnid;
my %uni_hog_map;
while(<$agesfh>){
	next if $.==1;
	chomp;
	(my $uniprot_id, my $refseq, my $fam_age) = split("\t");
	$fam_ages{$uniprot_id}=$fam_age;
}
close $agesfh;
print "Got family ages!\n";

open(my $accessionsfh, "<", $uni_hog_mapfile) or die $!;
while(<$accessionsfh>){
	chomp;
	next unless $_;
	(my $uni_id, my $hog_id) = (split(/\s+/,$_))[0,3];
	$uni_hog_map{$hog_id}=$uni_id;
}
close $accessionsfh;
print "Got the uniprot -> hog map!\n";

# my $treefolder = 'ALN_TREE';
# opendir(my $trees, $treefolder) or die $!;
# my @treefiles = grep{!/^\./ && /\.phb$/} readdir($trees);
# closedir $trees;

open(my $treefh,"<",'hogenom6.phyml.txt') or die $!;
my @treefilearray;
while(<$treefh>){
	chomp;
	push @treefilearray, (split(m/\s+/,$_))[1];
}
close $treefh;
open(my $logfh,">>",'distance_log.txt') or die $!;

my $counter = 0;
for my $treefileline (@treefilearray){
	$counter++;
	next unless $treefileline;
	next if $treefileline eq "(NO";
	open(my $treelinefh,"<",\$treefileline) or die $!;
	my $treeio = Bio::TreeIO -> new(-format => 'newick',
									-fh => $treelinefh );
	eval{
		while( my $tree = $treeio->next_tree ){
			my @leaves = $tree->get_leaf_nodes;
			my @human_leaves;
			my $root_node = $tree->get_root_node;
			my $root_length = $root_node->height();
			next if $root_length==0;
			my @lengths;
			
			for my $leaf (@leaves){
				my $id = $leaf->id;
				my $depth = $leaf->depth;
				(my $other,$id) = split(/\|/,$id);
				print {$logfh} ($other,'|',$id,"\n") if $other ne $id;
				my $prefix = (split('_',$id))[0];
				# print "id:$id\n";
				# print "prefix:$prefix\n";
				# print "OK!\n";
				# exit;
				push @lengths, $depth;
				push @human_leaves, $leaf if $prefix =~ m/^HS.*/;
			}
			my $min_length = min(@lengths);
			
			for my $leaf (@human_leaves){
				
				my $id = $leaf -> id();
				$id = (split(/\|/,$id))[1];	
				
				next unless(exists $uni_hog_map{$id});
				my $uniprot_id = $uni_hog_map{$id};
				next unless(exists $fam_ages{$uniprot_id});
				my $fam_age = $fam_ages{$uniprot_id};

				my $parent_node = $leaf->ancestor;
				my $depth = $parent_node->depth;
				$len_ages{$uniprot_id} = (1 - ($depth/$root_length))*($fam_age);
				$rates{$uniprot_id} = ($root_length-$min_length)/($fam_age);
				#print(join("\t",$id,$fam_age,$len_ages{$uniprot_id}),"\n");
				# print "uni_id:$uniprot_id\n";
				# print "family_age:$fam_age\n";
				# print "len_age:$len_ages{$uniprot_id}";
				# exit;			
				
			}

			
		}
	};
	print $@."\nError processing line $counter!\n" if $@;
}
close $logfh;
open($agesfh,"<",$agesfile) or die $!;
open(my $newfh,">>","rates_len_".$agesfile) or die $!;
print{$newfh}(join("\t",qw(uniprot_id refseq_id fam_age len_age rate)),"\n");

while(<$agesfh>){
	next if $.==1;
	chomp;
	my $id = (split("\t",$_))[0];
	$len_ages{$id} = "NA" unless exists $len_ages{$id};
	$rates{$id} = "NA" unless exists $rates{$id};
	print {$newfh} (join("\t",$_,$len_ages{$id},$rates{$id}),"\n");
	}
close $agesfh;
close $newfh;
exit;
