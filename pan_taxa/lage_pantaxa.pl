#!usr/bin/perl
#
#perl error_ager.pl InputFolder outfile_name.txt
#
#Takes the tree archive from release 22 (compara release 75)
#and gets age & length-based age & rate
#see ftp://ftp.ensemblgenomes.org/pub/pan_ensembl/release-22/mysql
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

my $infolder = shift;
my $outfilename = shift;
my %time_machine;
my %fam_ages;
my %len_ages;
my %rates;
my $target_species = 'homo sapiens';

#Navigate the directory structure
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
			#Take a copy of the tree so we don't break the original
			my $tree_clone = $original_tree->clone();
			my @result = age_tree($tree_clone, $target_species);
			next unless $result[0];
			for (my $i = 0; $i<$#result; $i+=2){
				my $node = $result[$i];
				next unless $node;
				(my $id, my $age, my $len_age, my $rate) = split('/',$result[$i+1]);
				$fam_ages{$id} = $age;
				$len_ages{$id} = $len_age;
				$rates{$id} = $rate;
    		}
		}
	closedir $innerdirh;
	}
}

#Print to output
open(my $newfh,">>",$outfilename) or die $!;
print{$newfh}(join("\t",qw(ens_id fam_age len_age rate)),"\n");
for my $id (sort keys %fam_ages){
	print {$newfh} (join("\t",$id,$fam_ages{$id},$len_ages{$id},$rates{$id}),"\n");
	}
close $newfh;

#

#

#

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
#it was broken for our purposes
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

#the sub to get divergence times. Takes the species names.
#returns -1 on an error, else returns the time.
sub get_time {
	local $| = 1;
	my $species_one = shift @_;
	my $species_two = shift @_;
	return 0 if lc($species_one) eq lc($species_two);
	my $website = 'http://timetree.org/mobile_query.php?taxon_a=' . $species_one . '&taxon_b=' . $species_two;
	my $output =get($website) or die "Unable to access $website!";
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
	#Chauvenet's Criterion is bugged for the web api, so let's fix it.
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

