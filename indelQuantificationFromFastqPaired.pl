#!/usr/bin/perl -w
################################################################################
# Quantitate InDels from Paired End sequenced (Illumina) reads.
# Steps:
#	1. Trim reads to exclude the adapters on the other end using Prinseq-Lite
#	2. Quality control (quality based trimming) using TrimGalore
#	3. Merge paired end reads using FLASH
#	4. Generate Alignment using BWA (SAM format)
#	5. Prepare (convert SAM->BAM, sort) alignment (BAM) files using Samtools
#	6. Calculate indels in the cut-site; export custom alignment
#		for reads having indels in the cut-site, in SAM format
#
# Author: Piyush Ranjan (piyuranjan@gatech.edu)
# Date created: May 08, 2014
# Last updated: ************
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday tv_interval);
use Sort::Naturally;

##All function definitions
#sub GetLoggingTime;

sub Steps
	{
	print '
################################################################################
# Quantitate InDels from Paired End sequenced (Illumina) reads.
# Steps:
#	1. Trim reads to exclude the adapters on the other end using Prinseq-Lite
#	2. Quality control (quality based trimming) using TrimGalore
#	3. Merge paired end reads using FLASH
#	4. Generate Alignment using Burrows-Wheeler Aligner (SAM format)
#	5. Prepare (convert SAM->BAM, sort) alignment (BAM) files using Samtools
#	6. Calculate indels in the cut-site; export custom alignment
#		for reads having indels in the cut-site, in SAM format
#
# Author: Piyush Ranjan (piyuranjan@gatech.edu)
# Date created: May 08, 2014
# Last updated: ************
################################################################################
	';
	}
sub Usage
	{
	print "
Usage:\n$0 [options: provide all flags separately]
\nPackage requirements (in PATH):\nPrinseq-Lite(prinseq-lite.pl); TrimGalore(trim_galore); FLASH(flash); Burrows-Wheeler Aligner(bwa); Samtools(samtools)
\nOptions:
 -d|dataDir	provide directory path with paired end read files
			This option assumes paired end files to be ordered
			next to each other.
			Default: pwd
 -w|workDir	provide path to create all intermediate/result files
			Default: pwd/IndelAnalysis
 -dW|delWorkDir	Delete previous workDir if exists (has same name)
			Default: off (Quit with error)
 -h|help	Print detailed help with steps involved
	\n";
	}

#Steps; Usage;
#print GetLoggingTime();

my ($dataDir,$workDir,$delWorkDir);
my ($help)=(0) x 1; #all 0 valued scalars
if(!GetOptions('d|dataDir=s' => \$dataDir,
				'w|workDir=s' => \$workDir,
				'dW|delWorkDir' => \$delWorkDir,
				'h|help' => \$help)
				||(!defined $dataDir))
	{Usage; exit 1;}

my %pairedFqFiles=%{ScanSeqFiles($dataDir)};
foreach my $key(nsort keys %pairedFqFiles)
	{print "$key\t$pairedFqFiles{$key}\n";}
##Write later: Input method by a config file

#####################################
######### Subroutines Below #########
#####################################

sub GetLoggingTime #find and return current logging timestamp
	{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $niceTimestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
	return $niceTimestamp;
	}

sub ScanSeqFiles #Scan and guess paired end files from a given location
	{
	my $dataDir=$_[0]; #record path for data files
	my %pairedFqFiles; #will record key->value as pair1File->pair2File
	opendir(DIR,$dataDir) or die $!;
	my @seqFiles=nsort(readdir(DIR)); #sorted aphanumerically
	my @fqFiles;
	foreach my $seqFile(@seqFiles)
		{
		next unless(-f "$dataDir/$seqFile"); #skip any directory/hidden files
		next unless($seqFile=~/\.f(ast)?q$/i); #skip files not with extension fastq/fq
		#print $seqFile."\n";
		push(@fqFiles,$seqFile);
		}
	for(my $i=0;$i<=$#fqFiles;$i+=2)
		{
		$pairedFqFiles{$fqFiles[$i]}=$fqFiles[$i+1];
		}
	return \%pairedFqFiles;
	}
