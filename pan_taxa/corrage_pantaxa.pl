#!usr/bin/perl
#
#perl error_ager.pl InputFolder member.txt gene_tree_node.txt gene_tree_node_attr.txt gene_tree_root.txt
#
#Takes the tree archive from release 22 (compara release 75), member.txt, gene_tree_node.txt, gene_tree_node_attr.txt, and gene_tree_root.txt tables
#and gets age & length-based age & metrics describing the error of the lineage
# see ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql
#

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

my $outfilename = 'pantaxa_erraged.txt';
my $deleted_filename = 'problem_genes.txt';
my $infolder = shift;
my $member_file = shift;
my $gene_tree_node_file=shift;
my $gene_tree_node_attr_file=shift;
my $gene_tree_root_file=shift;
my $counter = 0;
my $counter_two = 0;
my $counter_three = 0;
my $counter_four = 0;
my $counter_dupes = 0;
my @problem_ids;
my %time_machine;
my %id_member;		#stable id => member id			###ONLY USED TO GET NODE ID
my %member_node;	#member id => node id			###ONLY USED TO GET NODE ID
my %id_node_id;
my %node_parent;	#node id => parent node id		Parent ID can be \N if the node is a root
my %node_attr;		#node id => attributes
my %node_boot;		#node id => bootstrap			Below 70 often regarded as unreliable SEE http://sysbio.oxfordjournals.org/content/42/2/182.short
my %id_boot;
my %corrected_age;	#id => fixed age
my %corrected_boot;	#id => fixed root bootstrap
my %corrected_len_age;
my %corrected_rate;
my %root_type;
my %node_root;		#node id => root node id
my %fam_ages;		#								The table has 6903769 non-null bootstrap values.
my %len_ages;		#								3581244 of these are over 70.
my %rates;			#
my $min_cutoff = 60;
my $max_cutoff = 70;
my $target_species = 'homo sapiens';


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
	(my $root_id, my $clusterset_id) = (split("\t",$_))[0,3];
	next unless $clusterset_id eq 'default';
	$root_type{$root_id} = $clusterset_id;
}
close $gene_tree_root_fh;

open(my $gene_tree_node_fh,"<",$gene_tree_node_file) or die $!;
while(<$gene_tree_node_fh>){
	chomp;
	(my $node_id, my $parent_id, my $root_id, my $member_id) = (split("\t",$_))[0,1,2,6];
	next if $node_id eq '\N';
	next unless $root_type{$root_id};
	$member_node{$member_id}=$node_id;
	$node_parent{$node_id}=$parent_id;
	$node_root{$node_id} = $root_id;
}
close $gene_tree_node_fh;

for my $id (sort keys %id_member){
	my $node_id = $member_node{$id_member{$id}};
	$id_node_id{$id} = $node_id;
}
undef %id_member;
undef %member_node;

open(my $gene_tree_node_attr_fh,"<",$gene_tree_node_attr_file) or die $!;
while(<$gene_tree_node_attr_fh>){
	chomp;
	(my $node_id, my $node_type, my $bootstrap) = (split("\t",$_))[0,1,3];
	$node_attr{$node_id} = $node_type;
	next if $bootstrap eq '\N';
	$node_boot{$node_id} = $bootstrap;
}
close $gene_tree_node_attr_fh;


opendir(my $outerdirh,$infolder) or die $!;
while(readdir $outerdirh){
	next if $_ =~ m/^\./;
	my $directory = "$infolder/$_";
	next unless -d "$directory";
	opendir(my $innerdirh,$directory);
	while(readdir $innerdirh){
		next if $_ =~ m/^\./;
		my $treeio = Bio::TreeIO -> new(-format => 'phyloxml',
										-file => $directory."/".$_);
		while( my $original_tree = $treeio->next_tree ){
			my $tree_clone = $original_tree->clone();
			my @result = age_tree($tree_clone, $target_species);
			next unless $result[0];
			for (my $i = 0; $i<$#result; $i+=2){
				my $tree = $original_tree->clone();
				my @result = age_tree($tree, $target_species);
				next unless $result[0];
				my $node = $result[$i];
				next unless $node;
				my @cutoffs = get_linear_cutoffs($tree,$node,$min_cutoff,$max_cutoff);
				(my $id, my $age, my $len_age, my $rate) = split('/',$result[$i+1]);
				$fam_ages{$id} = $age;
				$len_ages{$id} = $len_age;
				$rates{$id} = $rate;
				my $node_id = $id_node_id{$id};
				my $root_id = $node_root{$node_id};
				my $rootstrap = $node_boot{$root_id};
				$id_boot{$id} = $rootstrap;
				my $cutoff_ref = \@cutoffs;
				my @error_id = find_error($node_id, $node, $cutoff_ref, 0);
				my $fixed_age = 'ER';
				my $fixed_len_age = 'ER';
				my $fixed_rate = 'ER';
				my $fixed_boot = 'ER';
				if ($error_id[0] == -1){
					$fixed_age = 'NA';
					$fixed_len_age = 'NA';
					$fixed_rate = 'NA';
					$fixed_boot = 'NA';
				}
				elsif($error_id[0] == 0){
					$fixed_age = $age;
					$fixed_len_age = $len_age;
					$fixed_rate = $rate;
					$fixed_boot = $rootstrap;
				}
				else{
					$counter++;
					(my $corrected_root_id, my $corrected_root_node) = @error_id;
					$fixed_boot = $node_boot{$corrected_root_id};
					$corrected_root_node -> ancestor(undef);
					my $subtree = Bio::Tree::Tree->new(-root => $corrected_root_node);
					####
					#print tree to see if it doesn't work
					####
					if (!$subtree){
						$counter_three++;
						next;
					}
					my @fixed_result = age_tree($subtree,$target_species);
					if ($fixed_result[0]==0){
						$counter_four++;
						$fixed_age = 0;
						$fixed_len_age = 0;
						$fixed_rate = 0;
						$fixed_boot = 0;
						push @problem_ids, $id;
					}
					$counter_two++;
					for (my $j = 0; $j<$#fixed_result; $j+=2){
						(my $corr_id, my $corr_age, my $corr_len_age, my $corr_rate) = split('/', $fixed_result[$j+1]);
						next unless $id eq $corr_id;
						$fixed_age = $corr_age;
						$fixed_len_age = $corr_len_age;
						$fixed_rate = $corr_rate;
					}
				}

				$corrected_age{$id} = $fixed_age;
				$corrected_boot{$id} = $fixed_boot;
				$corrected_len_age{$id} = $fixed_len_age;
				$corrected_rate{$id} = $fixed_rate;
			}
		}
	}
	closedir $innerdirh;
}
closedir $outerdirh;
open(my $delfh,">>",$deleted_filename) or die $!;
for my $id (@problem_ids){
	print {$delfh} $id , "\n";
	}
close $delfh;
	
open(my $newfh,">>",$outfilename) or die $!;
print{$newfh}(join("\t",qw(ens_id fam_age len_age rate root_boot corr_age corr_len_age corr_rate corr_root_boot)),"\n");
for my $id (sort keys %fam_ages){
	$corrected_boot{$id} = '' unless exists $corrected_boot{$id};
	$id_boot{$id} = '' unless exists $id_boot{$id};
	print {$newfh} (join("\t",$id,$fam_ages{$id},$len_ages{$id},$rates{$id},$id_boot{$id},$corrected_age{$id},$corrected_len_age{$id},$corrected_rate{$id},$corrected_boot{$id}),"\n");
	}
close $newfh;
say("I've found " . $counter . ' ages that need correcting.');
say("I was able to repair " . $counter_two . ".");
say("There were " . $counter_three . " skipped due to problems with the subtree.");
say("There were " . $counter_four . " (".scalar(@problem_ids).") given 0's due to problems aging the subtree.");
say("I made exceptions for " . $counter_dupes . " nodes because they were duplications.");
#

#

#

#recursive sub climbs the tree to find any bad parents.
#returns -1 if the initial node is bad, 0 if all nodes are good, and an array with the node id of the good node before the 1st bad one + the node itself if a bad one is found.
#find_error($leaf_node_id, \$node, @cutoffs,0);
sub find_error {
	local $|=1;
	my $node_id = shift;
	my $node = shift;
	my $cutoffs_ref = shift;
	my @cutoffs = @{$cutoffs_ref};
	my $indx = shift;
	my $boot_cutoff = $cutoffs[$indx];
	if ($node_parent{$node_id} eq '\N' || $node_root{$node_id} == $node_id){
		return(0);
	}
	if (! exists($node_parent{$node_id})){
		die "Couldn't find a parent for $node_id!"
	}
	if (! $node->ancestor){
		die "Node: " . $node->id . "\t" . $node_id . " encountered a problem.\nThis node has no ancestor, but the gene_tree_node table says there should be one!\n(And it should be: " . $node_parent{$node_id} . '.)';
	}
	return(-1) if is_bad($node_id,$boot_cutoff);		#This will only happen if the initial node is bad.
	
	my $parent_id = $node_parent{$node_id};
	if (is_bad($parent_id,$boot_cutoff)){
		return($node_id,$node);
		}
	else{
		my $parent_node = $node->ancestor;
		$indx++;
		return find_error($parent_id,$parent_node,$cutoffs_ref,$indx);
	}	
}

#sub checks if a node is dubious or has a bootstrap below the cutoff
#is_bad($node_id, $cutoff);
sub is_bad{
	my $node_id = shift;
	my $boot_cutoff = shift;
	if (((exists $node_attr{$node_id}) && ($node_attr{$node_id} eq 'dubious')) || ((exists $node_boot{$node_id}) && ($node_boot{$node_id} < $boot_cutoff))) {
		if ((exists $node_attr{$node_id}) && ($node_attr{$node_id} eq 'duplication')){
			$counter_dupes++;
			return 0;
		}
		else{
			return 1;
		}
	}
	else{
		return 0;
	}
}

#sub for determining variable cutoffs
#get_linear_cutoffs($tree, $node, $min, $max)
sub get_linear_cutoffs {
	my $tree = shift;
	my $node = shift;
	my $min = shift;
	my $max = shift;
	
	my @ancestors = $tree->get_lineage_nodes($node);
	my $num_nodes = scalar(@ancestors) + 1;
	my @cutoffs;
	my $increment = ($max-$min)/$num_nodes;
	push @cutoffs, $max;
	for my $indx (0..($num_nodes - 1)){
		push @cutoffs, $cutoffs[$indx] - $increment
	}
	return @cutoffs;
}

#sub for getting the age of genes of a particular species.
#age_tree($tree, 'homo sapiens');
#takes a BIO::Tree object and a species name
#returns key/values of the node followed by stats, where each value in the array is a different gene: gene_id/age/length-based age/rate
#or 0 if there are no genes.
sub age_tree {
	my $tree = shift;
	my $base_species = shift;
	my @result;
	
	my @leaves = $tree->get_leaf_nodes;
	my $root_node = $tree->get_root_node;
	my $root_length = $root_node->height();
	return 0 if $root_length == 0;
	my @ages;
	my @lengths;
	my @species;
	my @target_leaves;
	
	for my $leaf ( @leaves ) {
		push @lengths, override_depth($leaf,$root_node);
		my $string_leaf = '<xml>' . $leaf->to_string() . '</xml>';
		my $xml_ob = new XML::Bare( text=>$string_leaf );
		my $root = $xml_ob->parse();
		my $speciesname = $root->{xml}->{clade}->{taxonomy}->{scientific_name}->{value};
		unless (exists $time_machine{$speciesname}){
			$time_machine{$speciesname} = get_time($speciesname,$base_species);
		}
		my $age = $time_machine{$speciesname};
		push @ages, $age unless $age==-1;
		push @target_leaves, $leaf if lc($speciesname) eq lc($base_species);
	}
	return 0 if scalar(@target_leaves)==0;
	my $min_length = min(@lengths);
	my $fam_age = max(@ages);
	for my $leaf (@target_leaves){
		my $id = $leaf -> id();
		my $len = $leaf->branch_length();
		my $len_age = (1 - ($len/$root_length))*($fam_age);
		my $rate = 'NA';
		$rate = ($root_length-$min_length)/($fam_age) if $fam_age;
		my @partial_result = ($id, $fam_age, $len_age, $rate);
		push @result, $leaf;
		push @result, join('/',@partial_result);
	}
	return @result;
}

#returns the length from this node to the root.
#fixed it, tho.
#override_depth($node, $root_node);
sub override_depth{
	my $node = shift;
	my $root_node = shift;
	my $depth = 0;
	return 0 if $node == $root_node;
	do {
		$depth += $node->branch_length;
		$node = $node->ancestor;
	} while( $node != $root_node);
	return $depth;
}


#the sub to get times. Takes the species names.
#returns -1 on an error, else returns the time.
sub get_time {
	local $| = 1;
	my $species_one = shift @_;
	my $species_two = shift @_;
	return 0 if lc($species_one) eq lc($species_two);
	my $website = 'http://timetree.org/mobile_query.php?taxon_a=' . $species_one . '&taxon_b=' . $species_two;
	#print "\rOPENING website for $species\t\t";
	my $output =get($website) or die "Unable to access $website!";
	#print "\rGOT site for $species !\t\t";
	
	my $xml = new XML::Bare( text=>$output);
	my $root = $xml->parse();
	return -1 if exists $root->{result}->{error}->{value};
	my $avg = $root->{result}->{mean}->{avg}->{value};
	my $num_papers = $root->{result}->{all_total}->{value};
	my @times;
	if ($num_papers==1){
		push @times, $root->{result}->{time_estimates}->{node}->{time}->{value};
	}
	else{
		for my $paper (0..($num_papers-1)){
			push @times, $root->{result}->{time_estimates}->{node}->[$paper]->{time}->{value};
		}
	}
	#Chauvenet's Criterion is bugged for the web api
	if ($avg != sum(@times)/@times){
			my $stdev = 0;
			my $sigma = 0;
			if (@times > 5){
				my $cconstant = .9969 + (.4040 * log(scalar(@times)));
				for my $value (@times){
					$sigma += ($value - (sum(@times)/@times))**2;
				}
				$sigma /= scalar(@times);
				$stdev = $sigma ** .5;
				my $dmax = $cconstant * $stdev;
				for my $index (0..$#times){
					next unless $times[$index];
					if (($times[$index] < (sum(@times)/@times - $dmax)) || ($times[$index] > (sum(@times)/@times + $dmax))){
						splice(@times,$index,1);
						$index--;
					}
				}
			$avg = sum(@times)/@times;
			}
		}
	return $avg;
}	
exit;

__END__

Use Bio::Tree->clone($internal_node) to get a subtree to use when dubious nodes are detected.
Dubious are annotated as dubious OR bootstrap < 70

Make a function to do the normal tree aging
Then a function to do the error detection & aging; if error is detected, use clone to get a subtree and use the normal function on that

NOTE: We automatically accept null bootstrap values!!