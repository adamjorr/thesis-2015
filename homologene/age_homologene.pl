#!usr/bin/perl;
#perl age_homologene.pl

use strict;
use warnings;
use LWP::Simple;
use List::Util qw(sum max);
use Bio::DB::EUtilities;
#use feature qw(say);

my %taxa;
my %genes;
my %ages;
my $error = 0;
my $desired = 0;
my $speciestimedfile = 'timed_taxid_taxname';
my $outfile = 'homologene_aged.txt';
my %timer;

#Get IDs for Homo Sapien genes
my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                    -email => 'ajorr1@asu.edu',
                                    -db    => 'homologene',
									-term  => '"Homo sapiens"[Organism]',
									-retmax => 133548,
									-usehistory => 'y');
my $count = $factory ->get_count;
my $hist = $factory ->next_History || die 'No history returned!';
$factory -> set_parameters(-eutil => 'esummary',
                           -history => $hist);
 my ($retmax, $retstart) = (500,0);
 my $retry = 0;
#Get the docsum files
#to get taxa
RETRIEVE_IDS:
while($retstart < $count){
	$factory -> set_parameters(-retmax => $retmax,
							  -retstart => $retstart);
	eval{
		$factory -> get_Response(-cb => sub{
			while(my $docsum = $factory->next_DocSum()){
				my $id = $docsum->get_id;
				#Ignore the first header
				my $item = $docsum -> next_Item();
				my @taxids;
				my @geneids;
				my @genesyms;
				#Ignore the second headers
				while(my $subitem = $item -> next_Item('flattened')){
					#Get the third content if it's taxid or geneid
					while(my $subsubitem = $subitem -> next_Item){
						#print "Name: ",$subsubitem -> get_name, "\tContent: ",$subsubitem->get_content,"\n";
						push @taxids, $subsubitem -> get_content if $subsubitem -> get_name eq 'TaxId';
						push @geneids, $subsubitem -> get_content if $subsubitem -> get_name eq 'GeneID';
					}
				}
				#Add to hash
				$taxa{$id} = \@taxids;
				$genes{$id} = \@geneids;
			}
		});
	};
	if ($@) {
		die "Server error: $@. Try again\n" if $retry == 5;
		print STDERR "Server error, redo #$retry\n";
		$retry++ && redo RETRIEVE_IDS;
	}
	#say "Retrieved $retstart";
	$retstart += $retmax;
	#say "Size of hash: ",scalar(keys(%taxa));
}
#Get ages of each fam
open(my $famfh,"<",$speciestimedfile) or die $!;
while(<$famfh>){
	chomp;
	(my $taxonid, my $name, my $taxonage) = split("\t");
	$timer{$taxonid} = $taxonage;
	}
close $famfh;

for my $fam (sort {$a <=> $b} keys(%taxa)){
	my @famtaxa = @{$taxa{$fam}};
	my @famages;
	for my $eachtaxa (@famtaxa){
		push @famages, $timer{$eachtaxa};
	}
	$ages{$fam} = max(@famages);
}
#Get all human genes so we can compare later.
$factory -> reset_parameters(-eutil => 'esearch',
							-email => 'ajorr1@asu.edu',
							-db    => 'gene',
							-term  =>'"Homo sapiens"[Organism]',
							-retmax => 220000);
my @humgeneids = $factory -> get_ids;
my %humgenes;
@humgenes{@humgeneids} = ();						 
my %refseq;
	
#Find out which human genes we need
my %hum_fam_needed;			#famid  => @humangeneids
my %hum_needed_names;		#geneid => genename
my %hum_needed_aliases;		#geneid => @genealiases
#%ages						#famid  => age
my %refseq_acc;				#refseq_gi => refseq_acc

for my $fam (sort {$a <=> $b} keys(%ages)){
	my @famgenes = @{$genes{$fam}};
	my @humfamgenes = grep {exists $humgenes{$_}} @famgenes;
	$hum_fam_needed{$fam} = \@humfamgenes;
	@hum_needed_names{@humfamgenes} = ();
}

#Fill hum_needed_names and hum_needed_aliases;
my @hum_ids = sort(keys %hum_needed_names);
$retry = 0;
(my $length,my $arystart) = (400,0);
$count = scalar(@hum_ids);

RETRIEVE_NAMES:
while($arystart < $count){
	my @hum_segment = @hum_ids[$arystart..($arystart+$length-1)];
	until (defined($hum_segment[$#hum_segment])){pop @hum_segment;} #so we have an array of only defined values
	$factory ->reset_parameters(-eutil => 'elink',
								-email=> 'ajorr1@asu.edu',
								-db    => 'nuccore',
								-dbfrom=> 'gene',
								-id => \@hum_segment,
								-correspondence=>1);
	while (my $ds = $factory->next_LinkSet) {
		if ($ds->get_link_name eq 'gene_nuccore_refseqrna'){
			my @rsid = $ds->get_submitted_ids;
			die "There's a problem with the number of submitted ids" if (scalar(@rsid) != 1);
			my $rsacc = join('/',$ds->get_ids);
			$refseq{pop(@rsid)}=$rsacc;
		}
	}	
	
	$factory -> reset_parameters( -eutil => 'esummary',
								-email => 'ajorr1@asu.edu',
								-db => 'gene',	
								-id => \@hum_segment);
	while (my $ds = $factory->next_DocSum) {
		my $id = $ds->get_id;
		while (my $item = $ds->next_Item('flattened'))  {
			#print $item->get_name,"\t",$item->get_content,"\n";
			$hum_needed_names{$id} = $item->get_content if $item->get_name eq 'Name';
			$hum_needed_aliases{$id} = $item->get_content if $item->get_name eq 'OtherAliases';
			$hum_needed_aliases{$id} = join('/',split(/,\s/,$hum_needed_aliases{$id})) if $hum_needed_aliases{$id}; #
			#if the alias field doesn't work, the if statement in the line above is why.
		}
		#print "My name is: ",$hum_needed_names{$id},"\n";
		#print "My aliases are: ",$hum_needed_aliases{$id},"\n";
		#sleep 10
	}

	$arystart += $length;
	#say("Retrieved $arystart");
	#say("Refseq length: ",scalar(keys(%refseq)));
	#my $ex = (keys(%refseq))[int(rand(scalar(keys(%refseq))))];
	#say("Here's an example: ",$ex,"\t",$refseq{$ex});
	#say("Key: ",(keys(%refseq))[0],"\tVal: ",$refseq{(keys(%refseq))[0]});
}

$retry = 0;
($length,$arystart) = (100,0);

$factory -> reset_parameters(-eutil => 'esummary',
						-email  => 'ajorr1@asu.edu',
						-db     => 'nuccore');
my @all_refseq;

for my $genes (sort(keys %refseq)){
	push @all_refseq, split('/',$refseq{$genes})# if defined $refseq{$genes};
}
#say("I have ",scalar(@all_refseq)," refseq GIs!");
$count = scalar(@all_refseq);
RETRIEVE_REFSEQ:				
while($arystart < $count){
		my @hum_part = @all_refseq[$arystart..($arystart+$length-1)];
		until (defined($hum_part[$#hum_part])){pop @hum_part;}
		$factory -> reset_parameters( -eutil => 'esummary',
									-email=> 'ajorr1@asu.edu',
									-db => 'nuccore',	
									-id => \@hum_part);
		while (my $ds = $factory->next_DocSum) {
			my $id = $ds->get_id;
			while (my $item = $ds->next_Item('flattened'))  {
				$refseq_acc{$id} = $item->get_content if $item->get_name eq 'Caption';
			}
		}
		$arystart += $length;
		#say("Retrieved $arystart");
		#say("Hash length: ",scalar(keys(%refseq_acc)));
		#my $ex = (keys(%refseq_acc))[int(rand(scalar(keys(%refseq_acc))))];
		#say("Here's an example: ",$ex,"\t",$refseq_acc{$ex});
	}	
	
#Print the gene ids, fam, name, aliases, and ages.
open(my $outfh,">>",$outfile) or die $!;
print{$outfh}(join("\t",qw(gene_id family_id refseq_acc name aliases age)),"\n");
for my $fam (sort {$a <=> $b} keys(%ages)){
	no warnings "uninitialized";
	for my $humgene (@{$hum_fam_needed{$fam}}){
	if ($hum_needed_aliases{$humgene}){
		print {$outfh}(join("\t",$humgene,$fam,join('/',@refseq_acc{split('/',$refseq{$humgene})}),$hum_needed_names{$humgene},$hum_needed_aliases{$humgene},$ages{$fam}),"\n");
		}
	else{
		print {$outfh}(join("\t",$humgene,$fam,join('/',@refseq_acc{split('/',$refseq{$humgene})}),$hum_needed_names{$humgene},'',$ages{$fam}),"\n");
		}
	}
}
close $outfh;
exit;
