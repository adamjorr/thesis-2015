#!usr/bin/perl
#time_species_code.pl species_code

use warnings;
use LWP::Simple;
use List::Util qw(sum);
local $| = 1;
my $specfile = shift @ARGV;
my $desired;
my $error;

open(SPECIES,"<",$specfile) or die $!;
open(OUT,">>","timed_".$specfile) or die $!;
while (<SPECIES>){
	chomp;
	next unless m/^\|/;
	my $species = (split(/\|/,$_))[4];
	$species =~ s/\s+$//;
	my $website = 'http://timetree.org/mobile_query.php?taxon_a='."$species".'&taxon_b=homo%20sapiens';
	#print "\rOPENING website for $species\t\t";
	my $output =get($website) or die "Unable to access $website!";
	#print "\rGOT site for $species !\t\t";
	open(TTXML,"<",\$output) or die $!;
	my $avg;
	my $numpapers;
	my @times;
	while(<TTXML>){
		chomp;
		if ($_ =~ m/<mean>/){
			$desired = 1;
			next;
			}
		if ($_ =~ m/<error>/){
			$avg = $';
			$error = 1;
			#print "\rEncountered an error for $species...\t\t\n";
			last;
			}
		if ($_ =~ m/<avg>/ && $desired){
			my $after = $';
			if ($after =~ m|</avg>|){
				$avg = $`;
				#print "\rFound time for $species!\t\t\n";
				next;
				}
			}
		if ($_ =~ m|</mean>|){
			$desired = 0;
			next;
			}
		if ($_ =~ m/<all_total>/){
			my $after = $';
			if ($after =~ m|</all_total>|){
				$numpapers=$`;
				next;
				}
			}
		if ($_ =~ m/<time>/g){
			my $after = $';
			if ($after =~ m|</time>|){
				push @times, $`;
				next;
				}
			}
		}
	next unless @times;
	close TTXML;
	unless ($error){
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
	}
	$error = 0;
	print {OUT} (join("\t",$species,$avg),"\n");
}
close SPECIES;
close OUT;
exit;