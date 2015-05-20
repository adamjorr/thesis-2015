#!usr/bin/perl
#perl age_hogenom.pl
#
#Searches a list of HOGENOM families and
#appends times to each human gene from
#a given list of species divergence times.

use strict;
use warnings;
use List::Util qw(max);
use File::Slurp;
use LWP::UserAgent;

local $\ = "\n";
local $| = 1;
my $timesfile = 'timed_species_code';
my $fam_hogenomid = 'FAM_HOGENOMID';
my $species_code = 'species_code';
my $acc_id5_id6 = 'acc_id5_id6';
my $outfile = 'hogenom_aged.txt';
my %timemachine; 	#species name => species divergence time

open(TIME,"<",$timesfile) or die $!;
#Get times and give to hash
#Build time machine
while(<TIME>){
	chomp;
	next unless $_;
	(my $name, my $time) = split("\t");
	$timemachine{$name} = $time;
	}
close TIME;
	
#Build HogenomID => Family (With FAM_HOGENOMID)
open(my $fam_id,"<",$fam_hogenomid) or die $!;
my %id_fam;
while(<$fam_id>){
	chomp;
	next unless $_;
	s/\s+$//;
	(my $fam, my $id) = split(/\s+/,$_,2);
	$id_fam{$id} = $fam;
}
close $fam_id;

#Build Family => \@HogenomIDs (With FAM_HOGENOMID)
open($fam_id,"<",$fam_hogenomid) or die $!;
my %fam_ids;
while(<$fam_id>){
	chomp;
	next unless $_;
	(my $fam, my $id) = (split(/\s+/,$_,3))[0,1];
	if (exists $fam_ids{$fam}){
		$fam_ids{$fam} .= ";$id";
	}
	else{
		$fam_ids{$fam} = $id;
	}
}
close $fam_id;
for my $fam (sort keys %fam_ids){
	my @ids = split(';',$fam_ids{$fam});
	$fam_ids{$fam} = \@ids;
}

#Build HogenomID => Species (With species_code)
open(my $id_specs_fh,"<",$species_code) or die $!;
my %id_specs;
while(<$id_specs_fh>){
	chomp;
	next unless m/^\|/;
	s/^\|//;
	s/\s+$//;
	s/\*//;
	(my $id, my $species) = (split('\|'))[0,3];
	$id_specs{$id} = $species;
}
close $id_specs_fh;

#Build Family =>@speciesNames => Age using time machine
my %family_age;
#for my $family (sort keys %fam_ids){
	#print "\nfamily: $family";
	#print "ids: ", join(',',@{$fam_ids{$family}});
	#print "prefixes :",join(',',(map({(split('_',$_))[0]} @{$fam_ids{$family}})));
	#print "species: ", join(',',@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})});
	#print "times: ",join(',',@timemachine{@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})}});
	#print "max: ",max(@timemachine{@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})}});
	#sleep 2;
#}


for my $family (sort keys %fam_ids){
	eval{
		local $SIG{__WARN__} = sub{};
		my @prefixes = map({(split('_',$_))[0]} @{$fam_ids{$family}});
		my @lengthy = grep {!exists $id_specs{$_}} @prefixes;
		my @ok = grep{exists $id_specs{$_}} @prefixes;
		for my $longone (@lengthy){
			my @shortened = grep {$longone =~ m/$_.*/} sort keys %id_specs;
			if (scalar(@shortened) == 1){
				push @ok, @shortened;
			}
			else{
				push @ok, $shortened[0];
			}
		}
		$family_age{$family} = max(@timemachine{@id_specs{@ok}});
	};
	#if ($@){
	#	print "\n$@";
	#	next unless (scalar(@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})}) < 6);
	#	print "family: $family";
	#	print "ids: ", join(',',@{$fam_ids{$family}});
	#	print "prefixes :",join(',',(map({(split('_',$_))[0]} @{$fam_ids{$family}})));
	#	print "species: ", join("\n",@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})});
	#	print "times: ",join(',',@timemachine{@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})}});
	#	print "max: ",max(@timemachine{@id_specs{map({(split('_',$_))[0]} @{$fam_ids{$family}})}});
	#	sleep 5;
	#}
}

#Build Hogenom ID => Age
my %id_age;
my @all_ids = sort keys %id_fam;
my @human_ids = grep {m/^HS\d+_/} @all_ids;
@id_age{@human_ids} = @family_age{@id_fam{@human_ids}};

#Build UniProt => HogenomID with acc_id5_id6
open(my $acc_id5_id6_fh,"<",$acc_id5_id6);
my %uni_id;
while(<$acc_id5_id6_fh>){
	chomp;
	(my $uni, my $id) = (split(/\s+/,$_))[0,3];
	$uni_id{$uni} = $id;
}
my @human_unis = grep {exists($id_age{$uni_id{$_}})} sort keys %uni_id;

#UniProt Conversion
my %uni_refseq;
my $base = 'http://www.uniprot.org';
my $tool = 'mapping';
my $contact = 'ajorr1@asu.edu';
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';
GETTING:
for my $gene (@human_unis){
	next unless $gene;
	my $params = { 
		from => 'ACC',
		to => 'REFSEQ_NT_ID',
		format => 'list',
		query => $gene
	};
	next unless $params;
	my $response = $agent->post("$base/$tool/", $params);
	while (my $wait = $response->header('Retry-After')) {
	  print STDERR "Waiting ($wait)...\n";
	  sleep $wait;
	  $response = $agent->get($response->base);
	}
	my $retry = 0;
	my $got_data;
	
	eval{
	$response->is_success ?
		$got_data = $response->content :
		die 'Failed, got ' . $response->status_line .
			' for ' . $response->request->uri . "\n";
	};
	
	if ($@){
		print "Something went wrong getting accessions for $gene. Trying again...\n$@";
		$retry++;
		die "I couldn't do it! $@" if $retry == 8;
		redo GETTING;
	}
	chomp $got_data;
	unless($got_data){$uni_refseq{$gene}="NA";}
	else{$uni_refseq{$gene}=join('/',split("\n",$got_data));}
}

#Display UniProt, UniProt => RefSeq , UniProt => HogenomID => Age
open(my $outfh,">>",$outfile) or die $!;
print{$outfh}(join("\t",qw(uniprot_id refseq_id age)));
for my $geneid (@human_unis){
	no warnings "uninitialized";
	print{$outfh}(join("\t",$geneid,$uni_refseq{$geneid},$id_age{$uni_id{$geneid}}));
}

exit;
__END__
