#!usr/bin/perl
#
#perl calibration.pl Compara.phyloxml_aa_trees.22 member.txt gene_tree_node.txt gene_tree_root.txt gene_tree_node_attr.txt calibration_info.txt
#
#Takes the tree archive from release 22 (compara release 75), member.txt, gene_tree_node.txt, and gene_tree_root.txt tables
#and gets age & length-based age & metrics describing the error of the lineage
# see ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql
#

use strict;
use warnings;

use strict;
use warnings;
use feature qw(say);
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::NodeI;
use Bio::Annotation::Collection;
use XML::Bare;
use LWP::Simple;
use List::Util qw(min max sum);
use Bio::EnsEMBL::Registry;

my $tree_folder = shift;
my $member_file_name = shift;
my $gene_tree_node_name = shift;
my $gene_tree_root_name = shift;
my $gene_tree_node_attr_name = shift;
my $calibration_file = shift;
my $outfile = 'calibration_out.txt';

#print count_trees($tree_folder);
#exit;
(my $id_nodeid_file_ref, my $node_boot_parentid_ref) = load_tables($member_file_name, $gene_tree_root_name, $gene_tree_node_name, $gene_tree_node_attr_name);
say('Loaded tables');
load_file($calibration_file,$id_nodeid_file_ref,$node_boot_parentid_ref,$tree_folder,$outfile);


#load_file($filename, \%gene_file);
#loads file $filename and processes it
sub load_file{
	my $file_name = shift;
	my $id_hash_ref = shift;
	my $node_hash_ref = shift;
	my $tree_folder = shift;
	my $outfile = shift;
	open(my $infh , "<", $file_name) or die $!;
	process_file($infh, $id_hash_ref,$node_hash_ref, $tree_folder, $outfile);
	close $infh;
}

#process_file($filehandle)
#processes calibration info given in $filehandle
#
sub process_file{
	my $filehandle = shift;
	my $id_hash_ref = shift;
	my $node_hash_ref = shift;
	my $tree_folder = shift;
	my $outfile = shift;
	while(<$filehandle>){
		next if $.==1;
		chomp;
		(my $gene, my $species_list) = split("\t",$_);
		my @species = split(',',$species_list);
		my @species_lc = map {lc($_)} @species;
		$gene = map_hugo($gene);
		my @max_cutoffs = calibrate_gene($gene, $id_hash_ref, $node_hash_ref, $tree_folder, \@species_lc);
		next unless @max_cutoffs;
		my %out_hash;
		for my $i (0..$#max_cutoffs){
			$out_hash{$i/$#max_cutoffs} = $max_cutoffs[$i];
		}
		print_hash(\%out_hash,$outfile);
	}
}

#print_hash($hashref,$filename)
#prints a hash (tab delimited)
sub print_hash{
	my $hashref = shift;
	my $outfile = shift;
	my %hash = %$hashref;
	open(my $outfh,">>",$outfile) or die $!;
	map { print {$outfh} (join("\t",$_,$hash{$_}),"\n") } sort keys %hash;
	close $outfh;
}


#translate_gene
#given a gene, translates to ensembl id
sub map_hugo{
	my $hugo_gene_id = shift;
	my $ensembl_registry = 'Bio::EnsEMBL::Registry';
	$ensembl_registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');
	my $gene_adaptor = $ensembl_registry->get_adaptor('Human', 'Core', 'Gene');
	my $gene = $gene_adaptor->fetch_by_display_label($hugo_gene_id);
	die "Couldn't find $hugo_gene_id" unless $gene;
	return $gene->stable_id();
}

#calibrate_gene
#given a gene and list of species, returns maximum bootstrap cutoffs to get there
sub calibrate_gene{
	my $gene_id = shift;
	my $id_hash_ref = shift;
	my $node_hash_ref = shift;
	my $tree_folder = shift;
	my $speciesref = shift;
	my @species = @$speciesref;
	@species = sort(@species);
	my @unwanted_species = grep (/^!.+/,@species);
	@unwanted_species = map { $_ = substr($_,1) } @unwanted_species;
	@species = grep (!/^!.+/,@species);
	my @descendent_undesired;
	
	my @bootstrap_values;
	my $file = $id_hash_ref->{$gene_id}->{'file'};
	unless ($file){
		warn "Couldn't find the file for $gene_id";
		return;
	}
	my $full_file = find_full_path($tree_folder, $file);
	unless ($full_file){
		warn "Couldn't find the file for $gene_id : $file";
		return;
	}
	#say("Found the file for $gene_id ($full_file)");								#DEBUG
	my $treeio = Bio::TreeIO -> new(-format => 'phyloxml',
									-file => $full_file);
	while( my $tree = $treeio->next_tree ){
		for my $leaf ($tree->get_leaf_nodes){
			next unless lc($leaf->id) eq lc($gene_id);
			my $node = $leaf;
			my $node_id = $id_hash_ref->{$gene_id}->{'nodeid'};
			while($node->ancestor){
				$node = $node->ancestor;
				$node_id = $node_hash_ref->{$node_id}->{'parentid'};
				my $bootstrap_val = $node_hash_ref->{$node_id}->{'boot'};
				push @bootstrap_values, $bootstrap_val;
				my @descendent_leaves = grep {$_->is_Leaf} $node->get_all_Descendents;
				my @descendent_species = map{ lc(get_species($_)) } @descendent_leaves;
				my @descendent_desired = grep{ $_ ~~ @species } @descendent_species;
				@descendent_undesired = grep {$_ ~~ @unwanted_species} @descendent_species;
				my %unique_descendents;
				@unique_descendents{@descendent_desired} = 1;
				@descendent_desired = sort keys %unique_descendents;
				if (@descendent_desired ~~ @species){
					#say("Target species are: " . join("\t",@species));
					#say("I have: " . join("\t",@descendent_desired));
					last;
				}
				if ( ! $node->ancestor ){
					warn("There weren't enough species in $gene_id in $full_file");
					say("Target species are: " . join("\t",@species));
					say("I have: " . join("\t",@descendent_desired) . "\n");
					undef @bootstrap_values;
				}
			}
			say("I found some species that you didn't expect to be in this lineage!") if @descendent_undesired;
			if (@bootstrap_values){
				return @bootstrap_values;
			}
		}
	}
	
return;	
	
}

#get_species($leaf)
#gets the species name in $leaf
sub get_species{
	my $leaf = shift;
	my $string_leaf = '<xml>' . $leaf->to_string() . '</xml>';
	my $xml_ob = new XML::Bare( text=>$string_leaf );
	my $root = $xml_ob->parse();
	my $speciesname = $root->{xml}->{clade}->{taxonomy}->{scientific_name}->{value};
	return $speciesname;
}

#load_tables($member_file_name, $gene_tree_root_file_name, $gene_tree_node_file_name, $gene_tree_node_attr_file)
#loads tables in files. these are from pan ensembl release 22 (based on ensembl release 75)
#returns 2 hashes which map: stableids to node ids + file names, node ids to bootstrap values + parent ids.
sub load_tables{
	
	my $member_file = shift;
	my $gene_tree_root_file = shift;
	my $gene_tree_node_file = shift;
	my $gene_tree_node_attr_file = shift;
	
	my %id_member;
	my %root_type;
	my %member_node;
	my %node_name;
	my %id_node_id;
	my %id_nodeid_file;
	my %node_boot_parentid;
	my %node_boot;
	my %node_parent;
	my %node_root;
	
	open(my $memberfh,"<",$member_file) or die $!;
	while(<$memberfh>){
		chomp;
		(my $stab_id, my $txn_id, my $canonical_id) = (split("\t",$_))[1,4,8];
		next unless $txn_id==9606;
		$id_member{$stab_id}=$canonical_id;
	}
	close $memberfh;
	
	open(my $gene_tree_root_fh,"<",$gene_tree_root_file) or die $!;
	while(<$gene_tree_root_fh>){
		chomp;
		(my $root_id, my $type, my $clusterset_id, my $tree_name) = (split("\t",$_))[0,2,3,7];
		next unless $clusterset_id eq 'default';
		next unless $type eq 'tree';
		$root_type{$root_id} = $clusterset_id; 
		$node_name{$root_id} = $tree_name;
	}
	close $gene_tree_root_fh;
	  
	open(my $gene_tree_node_fh,"<",$gene_tree_node_file) or die $!;
	while(<$gene_tree_node_fh>){
		chomp;
		(my $node_id, my $parent_id, my $root_id, my $member_id) = (split("\t",$_))[0,1,2,6];
		next if $node_id eq '\N';
		next unless $root_type{$root_id};
		$node_root{$node_id} = $root_id;
		$member_node{$member_id}=$node_id;
		$node_parent{$node_id} = $parent_id;
		}
	close $gene_tree_node_fh;
	
	open(my $gene_tree_node_attr_fh,"<",$gene_tree_node_attr_file) or die $!;
	while(<$gene_tree_node_attr_fh>){
		chomp;
		(my $node_id, my $bootstrap) = (split("\t",$_))[0,3];
		next if $bootstrap eq '\N';
		$node_boot{$node_id} = $bootstrap;
	}
	close $gene_tree_node_attr_fh;
	map{$id_node_id{$_} = $member_node{$id_member{$_}}} sort keys %id_member;
	for my $id (sort keys %id_member){
		next unless $id;
		next unless $id_node_id{$id};
		next unless exists $node_root{$id_node_id{$id}};
		next unless exists $node_name{$node_root{$id_node_id{$id}}};
		my $node_id = $id_node_id{$id};
		my $filename = $node_name{$node_root{$node_id}} . '.aa.xml'; #$tree_folder_path
		$id_nodeid_file{$id}{'nodeid'} = $node_id;
		$id_nodeid_file{$id}{'file'} = $filename;
	}
	
	for my $nodeid (sort keys %node_parent){
		my $boot = $node_boot{$nodeid};
		my $parent = $node_parent{$nodeid};
		$node_boot_parentid{$nodeid}{'boot'} = $boot;
		$node_boot_parentid{$nodeid}{'parentid'} = $parent;
	}
	
	
	
	# map{$id_filename{$_} = $tree_folder_path . "/" . $node_name{$id_node_id{$_}}} sort keys %id_node_id;
	
	undef %id_member;
	undef %root_type;
	undef %member_node;
	undef %node_name;
	undef %id_node_id;
	undef %node_boot;
	undef %node_parent;
	undef %node_root;
	
	return (\%id_nodeid_file, \%node_boot_parentid);
	
}

sub find_full_path{
	my $tree_folder_path = shift;
	my $filename = shift;
	my $desired;
	opendir(my $dirh, $tree_folder_path) or die $!;
	while(readdir $dirh){
		next if $_ =~ m/^\.+/;
		my $outer = $_;
		opendir(my $dirh2, $tree_folder_path . "/" . $outer) or die $!;
		while(readdir $dirh2){
			next if $_ =~ m/^\.+/;
			$desired = $outer . "/" . $_ if $_ eq $filename;
			last if $desired;
		}
		closedir $dirh2;
		last if $desired;
	}
	closedir $dirh;
	return $tree_folder_path . "/" . $desired if $desired;
	return;
}

sub count_trees{
	my $counter;
	my $tree_folder_path = shift;
	opendir(my $dirh, $tree_folder_path) or die $!;
	while(readdir $dirh){
		next if $_ =~ m/^\.+/;
		my $outer = $_;
		opendir(my $dirh2, $tree_folder_path . "/" . $outer) or die $!;
		while(readdir $dirh2){
			next if $_ =~ m/^\.+/;
			$counter++;
		}
		closedir $dirh2;
	}
	closedir $dirh;
	return $counter;
}

sub find_full_path2{
	my $tree_folder_path = shift;
	my $gene_id = shift;
	my $desired;
	opendir(my $dirh, $tree_folder_path) or die $!;
	while(readdir $dirh){
		next if $_ =~ m/^\.+/;
		opendir(my $dirh2, $tree_folder_path . "/" . $_) or die $!;
		my $outer = $_;
		while(readdir $dirh2){
			next if $_ =~ m/^\.+/;
			my $treeio = Bio::TreeIO -> new(-format => 'phyloxml',
											-file => $tree_folder_path."/".$outer."/".$_);
			while( my $tree = $treeio->next_tree ){
				my @leaves = $tree->get_leaf_nodes;
				for my $leaf (@leaves){
					$desired = $tree_folder_path."/".$outer."/".$_ if $leaf->id eq $gene_id;
				}
			}
			last if $desired;
		}
		closedir $dirh2;
		last if $desired;
	}
	closedir $dirh;
	return $desired if $desired;
	return;
}

#use gene_tree_root to translate tree ID to filename
#my %id_node_id; stable id to node id
#my %node_root; node id to root id
#gene_tree_root[0] is node id
#gene_tree_root[7] is the tree id + filename
#gene_tree_root[6] is the canonical id (it will be \N if this tree is the canonical one!)