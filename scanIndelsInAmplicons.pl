#!/usr/bin/perl -w
##############################################################
# Scan INDEL in amplicons from a custom SAM alignment file
# -uses site info for amplicons (tab separated text) with cols:
#	ampliconName,startInReference,cutSite,endInReference
# -outputs indel frequency in tab delimited text with cols:
#	ampliconName,#indelsInCutSite,#readsWithIndelsInCutSite,
#	#perfectMatch,#totalReadsAligned
##############################################################

use strict;
use warnings;
use List::Util qw(min max);
#use List::Util;

## CHANGE HERE IF WANT TO CHANGE THE RANGE ##
my $cutSiteRange=10;
print "\nCut site range is set to:\t$cutSiteRange\n";

if((!defined $ARGV[0])||(!defined $ARGV[1])||(!defined $ARGV[2])||(!defined $ARGV[3]))
	{die "\nUsage: $0 [BAM_AlignmentFile] [ampliconCutSiteReferenceFile] [outputIndelFreqFile] [referenceName(barcodeReferencePrefix/faidx-index-filename)]\n";}

open(AMP,$ARGV[1]) or die $!;
my %cutSites;
while(<AMP>) #read cutsites in Amplicons
	{
	next unless(/^\S+\t\d+\t\d+\t\d+\t\d+$/);
	chomp;
	my @line=split(/\t/);
	$cutSites{$line[0]}{amp}=$line[1];
	$cutSites{$line[0]}{ref}=$line[3];
	}
close(AMP);

if($ARGV[0]=~/sam$/) #convert to BAM if input is SAM
	{
	warn "\nYour input alignment looks like SAM\nConverting to BAM first...\n";
	my $outf=(split(/\./,$ARGV[0]))[0];
	$outf.='.bam';
	print `samtools view -bS -o $outf $ARGV[0]`;
	$ARGV[0]=$outf;
	}
unless(-e $ARGV[0].'.bai') #create index of the sam
	{
	warn "\nLooks like you don't have BAM index here needed for creating sub-alignments\nCreating index...\n";
	print `samtools index $ARGV[0]`;
	}

my $customAlignmentDir=$ARGV[3].'_AmpliconAlignments';
print "\nCreating directory for SAM files: $customAlignmentDir...";
print `mkdir $customAlignmentDir`;
print "\nCreating data sheet for indel frequency...";
open(OUT,">$ARGV[2]") or die $!;
print OUT "AmpliconName\t#TotalIndels\t#ReadsWithIndels\t%ReadsWithIndels\t#ReadsWithMatch\t#ReadsAlignedTotal\n";
foreach my $amplicon(sort(keys(%cutSites)))
	{
	#create custom alignment range (+-15 base of cutSite)
	my $begin=$cutSites{$amplicon}{ref}-$cutSiteRange;
	my $end=$cutSites{$amplicon}{ref}+$cutSiteRange;
	#create custm alignment SAM file
	my $customAlignmentFile="$customAlignmentDir\/$ARGV[3]\_$amplicon\_$begin-$end\.sam";
	my $indelAlignmentFile="$customAlignmentDir\/$ARGV[3]\_$amplicon\_$begin-$end\_indels\.sam";
	print "\nCreating custom alignment for $amplicon\:$begin-$end...";
	print `samtools view -ho $customAlignmentFile $ARGV[0] $ARGV[3]:$begin-$end`;
	
	###parsing sam input to check indels in cut site region###
	my $indelCount=0; my $indelReadCount=0; my $readCount=0;
	print "\nCreating SAM file for $amplicon indels...";
	open(SAMIN,$customAlignmentFile) or die $!;
	open(SAMOUT,">$indelAlignmentFile") or die $!; #create SAM with only indels
	while(<SAMIN>)
		{
		if(/^@/) #skip SAM headers to indels file
			{print SAMOUT $_;next;}
		my $cigar=(split(/\t/))[5]; #extract cigar code
		next if($cigar=~/H/); #skip all split alignments
		$readCount++; my $indelReadFlag=0;
		if($cigar=~/[DI]/) #look if cigar has Insertions or Deletions
			{
			## Processing CIGAR for counting if there is an indel in cut site ##
			my @nt=split(/\D/,$cigar); #note all nucleotides
			if(!defined $nt[$#nt])
				{pop @nt;}
			#print "!!@nt!\t";
			my @var=($cigar=~/(\D)/g); #note all variations
			my $offset=(split(/\t/))[3]; #extract offset position of read
			my $start=$offset;
			foreach my $feature(0..$#var) #span all features of the CIGAR code
				{
				my $stop=$start+$nt[$feature]-1; #stop for current feature
				if($var[$feature]=~/[DI]/) #look if current feature is an indel
					{
					my $overlap=min($stop,$end)-max($start,$begin); #Calculate overlap of the feature
					if($overlap>=0) #look if current indel is in cutsite
						{
						$indelCount++;
						$indelReadFlag=1; #flag the read as indel +ve
						}
					}
				$start=$stop+1 if($var[$feature]!~/I/); #start for next feature (will not change if current feature was an insertion- insert length adjustment)
				}
			}
		if($indelReadFlag) #if current read was indel +ve
			{
			$indelReadCount++;
			print SAMOUT $_;
			}
		}
	close(SAMIN);
	close(SAMOUT);
	my $matchCount=$readCount-$indelReadCount;
	my $indelReadPercent=($indelReadCount/$readCount)*100;
	print OUT "$amplicon\t$indelCount\t$indelReadCount\t";
	printf OUT ("%.3f",$indelReadPercent);
	print OUT "%\t$matchCount\t$readCount\n";
	}
close(OUT);
print "\n\nIndel quantitation finished...\n";
