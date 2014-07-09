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
# Last updated: Check on https://github.com/piyuranjan/NucleaseIndelActivityScript
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday tv_interval);
use Cwd;
use Sort::Naturally;

##All function definitions (for previous perl versions, prior to v5.10)
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
# Last updated: Check on https://github.com/piyuranjan/NucleaseIndelActivityScript
################################################################################
	';
	}
sub Usage
	{
	print "
Usage:\n$0 [options: provide all flags separately]
\nPackage requirements (in system PATH variable):\nPrinseq-Lite(prinseq-lite.pl); TrimGalore(trim_galore); FLASH(flash); Burrows-Wheeler Aligner(bwa); Samtools(samtools)
\nOptions:
 -c|config	[string:required] provide config file (with location)
			(pwd assumed as location if not given)
 -d|dataDir	[string] provide directory path which contains paired-end Read files,
			Reference sequence file(s) and tab-delimited Amplicon cut-site files.
			Program will find the filenames specified in the config, in the
			complete directory tree under given path.
			Default: pwd
 -w|workDir	[string] provide path to create all intermediate/result files
			result files will be created in path/IndelAnalysis.\$timestamp
			Default: pwd
 -v|verbose	Print all logging steps (on STDOUT)
 -debug		Turn on Debug mode with multiple logging values. Turn this on only
			if you are the developer/contributor to the program.
 -h|help	Print detailed help with steps involved
	\n";
	}

our @entryParameters; #this hash stores all the parameters given in config file
my ($configFile); #all undef vars
my $pwd=cwd();
our ($dataDir,$workDir)=($pwd) x 2; #default for dataDir and workDir
my ($help,$verbose,$debug)=(0) x 3; #all 0 valued scalars
#scan all command line options
if(!GetOptions('c|config=s' => \$configFile,
				'd|dataDir=s' => \$dataDir,
				'w|workDir=s' => \$workDir,
				'v|verbose' => \$verbose,
				'debug' => \$debug,
				'h|help' => \$help)
				||(!defined $configFile))
	{Usage; exit 1;}
if($help) #quit with help
	{Steps;Usage;exit 0;}
print "\n...Debug mode on...\n" if $debug;
print "\nTime before beginning: ".GetLoggingTime()."\n" if ($verbose||$debug);
$dataDir=~s/\/$//; #remove trailing slash if present
die "$dataDir: Directory absent / not readable\n" unless((-d $dataDir)&&(-r $dataDir));
$workDir=~s/\/$//; #remove trailing slash if present
$workDir.='/IndelAnalysis.'.GetLoggingTime(); #default for workDir

###display all configuration options###


##Deprecated: Scan fastq files on a location
# my %pairedFqFiles=%{ScanSeqFiles($dataDir)};
# foreach my $key(nsort keys %pairedFqFiles)
	# {print "$key\t$pairedFqFiles{$key}\n" if $debug;}

###Scan all parameters from the configuration file###
ScanParametersFromConfig($configFile);
die "\nError: No valid entries in the config file: $configFile\n" if($#entryParameters<0); #record check
mkdir $workDir or die $!; #creating work directory if all checks passed

###Processing each entry in config for Indel quantification###
foreach my $entry(0..$#entryParameters)
	{
	##create sub directories for each sample
	my $fqFile=(split(/\//,${$entryParameters[$entry]}{ReadPair1}))[-1];
	$fqFile=~/(\S+)\.f.*q/;
	my $entryDir="$workDir/$1";
	#print "$entryDir\n" if $debug;
	mkdir $entryDir or die $!;
	
	##trimming reads by length=avgReadLength-minAmpliconSize
	my $trimLength=${$entryParameters[$entry]}{AvgReadLength}-${$entryParameters[$entry]}{MinAmpliconLength};
	my $trimDir=$entryDir."/1.trimmedReads";
	mkdir $trimDir or die $!;
	#invoke function
	}


#####################################
######### Subroutines Below #########
#####################################

sub GetLoggingTime #find and return current logging timestamp
	{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $niceTimestamp = sprintf ( "%04d%02d%02d-%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
	return $niceTimestamp;
	}

## ScanSeqFiles deprecated now
sub ScanSeqFiles #Scan and guess paired end files from a given location
	{
	my $dataDir=$_[0]; #record path for data files
	my %pairedFqFiles; #will record key->value as pair1File->pair2File
	opendir(DIR,$dataDir) or die $!;
	my @seqFiles=nsort(readdir(DIR)); #sorted aphanumerically
	closedir(DIR);
	my @fqFiles;
	foreach my $seqFile(@seqFiles)
		{
		next unless(-f "$dataDir/$seqFile"); #skip any directory/hidden files
		next unless($seqFile=~/\.f(ast)?q$/i); #skip files not with extension fastq/fq
		#print $seqFile."\n" if $debug;
		push(@fqFiles,$seqFile);
		}
	for(my $i=0;$i<=$#fqFiles;$i+=2)
		{
		$pairedFqFiles{$fqFiles[$i]}=$fqFiles[$i+1];
		}
	return \%pairedFqFiles;
	}

sub ScanParametersFromConfig #Scans the configuration file for all the parameters and returns a hash with all parameter values
	{
	my $configFile=$_[0];
	die "\nError: Config file empty." if(-z $configFile);
	my @paramOrder; #array to store order of the values
	my %paramLabels=qw(ReadPair1 0 ReadPair2 0 AvgReadLength 0 MinAmpliconLength 0 ForwardAdapter 0 ReverseAdapter 0 ReferenceSeqFile 0 AmpliconCutSites 0); #hash to confirm all labels and store their location
	open(CONF,$configFile) or die "$configFile: $!\n";
	my $entryCounter=0;
	while(<CONF>)
		{
		next if(/(#|^$)/); #skip any comments/blank lines
		chomp();
		if(/^LABELS:/) #scan, confirm integrity, store location of all labels
			{
			@paramOrder=split(/\t/,$_);
			shift(@paramOrder); #remove "LABELS:" tag from param array
			die "\nError: All parameter labels not present in configFile: $configFile at line $.\n" if($#paramOrder<7); #integrity check
			for(my $param=0;$param<=$#paramOrder;$param++)
				{
				if(defined $paramLabels{$paramOrder[$param]}) #if name matches, store location
					{$paramLabels{$paramOrder[$param]}=$param;}
				else
					{die "\nError: parameter mismatch: no parameter matches label \"$paramOrder[$param]\" in configFile: $configFile at line $.\n";}
				}
			next;
			}
		my @paramValues=split(/\t/,$_);
		# foreach my $val(@paramValues)
			# {print "!$val\n" if $debug;}
		# print "\n@paramValues\n$#paramValues" if $debug;
		# print scalar @paramValues if $debug;
		die "\nError: All parameter values not present for entry at line $. in configFile: $configFile\n" if($#paramValues<7); #integrity check for all values in configFile
		#scan values and complete file paths
		my $readPair1=FindFilePath($paramValues[$paramLabels{ReadPair1}]);
		${$entryParameters[$entryCounter]}{ReadPair1}=$readPair1;
		#print ${$entryParameters[$entryCounter]}{ReadPair1}."!!\n" if $debug;
		my $readPair2=FindFilePath($paramValues[$paramLabels{ReadPair2}]);
		${$entryParameters[$entryCounter]}{ReadPair2}=$readPair2;
		my $referenceSeqFile=FindFilePath($paramValues[$paramLabels{ReferenceSeqFile}]);
		${$entryParameters[$entryCounter]}{ReferenceSeqFile}=$referenceSeqFile;
		my $ampliconCutSites=FindFilePath($paramValues[$paramLabels{AmpliconCutSites}]);
		${$entryParameters[$entryCounter]}{AmpliconCutSites}=$ampliconCutSites;
		#scan parameters
		${$entryParameters[$entryCounter]}{AvgReadLength}=$paramValues[$paramLabels{AvgReadLength}];
		${$entryParameters[$entryCounter]}{MinAmpliconLength}=$paramValues[$paramLabels{MinAmpliconLength}];
		${$entryParameters[$entryCounter]}{ForwardAdapter}=$paramValues[$paramLabels{ForwardAdapter}];
		${$entryParameters[$entryCounter]}{ReverseAdapter}=$paramValues[$paramLabels{ReverseAdapter}];
		$entryCounter++;
		}
	}

sub FindFilePath #looks for complete path of a file in dataDir
	{
	my $fileName=$_[0];
	my $command=`find $dataDir -name $fileName`;
	#print "xx$command\n" if $debug;
	my @results=split(/\n/,$command);
	if($#results>0)
		{
		warn "\nWarning: More than one file found by name: $fileName in $dataDir\n";
		foreach my $result(@results)
			{warn "$result\n";}
		warn "\nOnly $results[0] will be used for processing\n";
		}
	elsif($#results<0)
		{die "\nError: File: $fileName not found in directory tree under: $dataDir\n";}
	#print "$results[0]\n" if $debug;
	return $results[0];
	}

sub TrimReadsByLength #Trims reads by a given length criteria using prinseq-lite
	{
	my $outDir=$_[0]; #output files will be generated here
	my $trimLength=$_[1]; #reads will be trimmed by this length from 3' end
	my $entryNumber=$_[2]; #this record number in the config file will be processed
	
	}


