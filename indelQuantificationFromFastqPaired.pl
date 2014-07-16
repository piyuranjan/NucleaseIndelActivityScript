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
use List::MoreUtils qw(uniq);

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
 -r|cutRange	[int] InDels in this range (+/-) of the cut site (provided in config)
			will be calculated
			Default: 15
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
my $cutSiteRange=15; #default cut site range (+/-) of the given cut site NT
#scan all command line options
if(!GetOptions('c|config=s' => \$configFile,
				'd|dataDir=s' => \$dataDir,
				'w|workDir=s' => \$workDir,
				'r|cutRange=i' => \$cutSiteRange,
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

###Create Indexes for BWA###
my $indexDir="$workDir/ReferenceIndexes";
mkdir $indexDir or die $!; #creating directory for all reference indexes for BWA
CreateBWAIndex($indexDir);


###Processing each entry in config for Indel quantification###
foreach my $entry(0..$#entryParameters)
	{
	print "\n","=" x 50 if($verbose||$debug);
	print "\nProcessing ${$entryParameters[$entry]}{SampleName}\n" if($verbose||$debug);
	##create sub directories for each sample
	my $fqFile=(split(/\//,${$entryParameters[$entry]}{ReadPair1}))[-1];
	$fqFile=~/(\S+)\.f(ast)?q$/;
	my $entryDir="$workDir/${$entryParameters[$entry]}{SampleName}";
	#print "$entryDir\n" if $debug;
	mkdir $entryDir or die $!;
	
	##trimming reads by length=avgReadLength-minAmpliconSize
	my $trimLength=${$entryParameters[$entry]}{AvgReadLength}-${$entryParameters[$entry]}{MinAmpliconLength};
	my $trimDir=$entryDir."/1.trimmedReads";
	mkdir $trimDir or die $!;
	my @trimmedFiles=TrimReadsByLength($entry,$trimLength,$trimDir); #invoke function
	#print "\nFiles generated: @trimmedFiles\n" if $debug;
	
	##quality control and adapter trimming
	my $qualTrimDir=$entryDir."/2.highQualityReads";
	mkdir $qualTrimDir or die $!;
	my @qualTrimmedFiles=TrimReadsByQuality($entry,@trimmedFiles,$qualTrimDir);
	#print "\nFiles generated: @qualTrimmedFiles\n" if $debug;
	
	##merge read pairs to get longer reads
	my $mergeDir=$entryDir."/3.mergedReads";
	mkdir $mergeDir or die $!;
	my $mergedFile=MergePairedReads($entry,@qualTrimmedFiles,$mergeDir);
	
	##create alignments of the reads to amplicon sequences
	my $alignDir=$entryDir."/4.readAlignment";
	mkdir $alignDir or die $!;
	my $alignFile=CreateBWAAlignment($entry,$mergedFile,$alignDir);
	
	##quantify indels in the alignments generated for each amplicon
	my $indelDir=$entryDir."/indelQuantification";
	mkdir $indelDir or die $!;
	ScanIndels($entry,$alignFile,$indelDir,$cutSiteRange);
	
	print "\nProcessing for ${$entryParameters[$entry]}{SampleName} finished successfully\n" if($verbose||$debug);
	print "=" x 50,"\n" if($verbose||$debug);
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
	my %paramLabels=qw(SampleName 0 ReadPair1 0 ReadPair2 0 AvgReadLength 0 MinAmpliconLength 0 ForwardAdapter 0 ReverseAdapter 0 ReferenceSeqFile 0 AmpliconCutSites 0); #hash to confirm all labels and store their location
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
			die "\nError: All parameter labels not present in configFile: $configFile at line $.\n" if($#paramOrder<8); #integrity check
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
		die "\nError: All parameter values not present for entry at line $. in configFile: $configFile\n" if($#paramValues<8); #integrity check for all values in configFile
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
		${$entryParameters[$entryCounter]}{SampleName}=$paramValues[$paramLabels{SampleName}];
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
	my $entryCounter=$_[0]; #this record number in the config file will be processed
	my $trimLength=$_[1]; #reads will be trimmed by this length from 3' end
	my $outDir=$_[2]; #output files will be generated here
	my $prinSeqOut=`prinseq-lite.pl -fastq ${$entryParameters[$entryCounter]}{ReadPair1} -fastq2 ${$entryParameters[$entryCounter]}{ReadPair2} -trim_right $trimLength -out_good $outDir/${$entryParameters[$entryCounter]}{SampleName} -out_bad null 2>&1`; #command for trimming
	#print "return code = ", $?, "\n" if $debug;
	if($?) {die "$prinSeqOut\n$!";} #sanity check
	print "\nTrimming using prinseq-lite by $trimLength NT from 3' for: ${$entryParameters[$entryCounter]}{SampleName}\n$prinSeqOut\nTrimming done...\n" if $verbose;
	my @trimmedFiles=("$outDir/${$entryParameters[$entryCounter]}{SampleName}_1.fastq","$outDir/${$entryParameters[$entryCounter]}{SampleName}_2.fastq");
	if((-f $trimmedFiles[0])&&(-f $trimmedFiles[1]))
		{return @trimmedFiles;} #return names of the files generated
	else
		{die "Files not formed: @trimmedFiles\nCheck if any error in prinseq-lite\n";} #integrity check
	}

sub TrimReadsByQuality #Trims reads from 3' end by phred score 20 and removes adapters using TrimGalore
	{
	my $entryCounter=$_[0]; #this record number in the config file will be processed
	my $forFile=$_[1]; my $revFile=$_[2]; #sequences in these files will be trimmed from 3' end
	my $outDir=$_[3]; #output files will be generated here
	#print "\nentry:$entryCounter\nTrimmed Files: $forFile\t$revFile\noutDir:$outDir\n" if $debug;
	my $trimGaloreOut=`trim_galore -q 20 --paired -a ${$entryParameters[$entryCounter]}{ForwardAdapter} -a2 ${$entryParameters[$entryCounter]}{ReverseAdapter} -o $outDir --retain_unpaired $forFile $revFile 2>&1`; #command for quality trimming
	#print "return code = ", $?, "\n" if $debug;
	if($?) {die "$trimGaloreOut\n$!";} #sanity check
	print "\nQuality and adapter trimming using Trim_Galore by phred score 20 from 3' for: ${$entryParameters[$entryCounter]}{SampleName}\n$trimGaloreOut\nTrimming done...\n" if $verbose;
	my @qualTrimmedFiles=($forFile,$revFile);
	my $filePair=0;
	foreach my $file(@qualTrimmedFiles) #getting names of the newer files generated to return
		{
		$filePair++;
		$file=(split(/\//,$file))[-1];
		$file=~/(\S+)\.f(ast)?q$/;
		$file="$outDir/$1_val_$filePair.fq";
		}
	if((-f $qualTrimmedFiles[0])&&(-f $qualTrimmedFiles[1]))
		{return @qualTrimmedFiles;} #return names of the files generated
	else
		{die "Files not formed: @qualTrimmedFiles\nCheck if any error in trim_galore\n";} #integrity check
	}

sub MergePairedReads #merges read pairs in to a single read using FLASH
	{
	my $entryCounter=$_[0]; #this record number in the config file will be processed
	my $forFile=$_[1]; my $revFile=$_[2]; #sequences in these files will be trimmed from 3' end
	my $outDir=$_[3]; #output files will be generated here
	my $minOverlapToMerge=30; #this sets the min NT overlap to be matched for successful merging
	my $mergeOut=`flash -m $minOverlapToMerge -o ${$entryParameters[$entryCounter]}{SampleName} -d $outDir $forFile $revFile 2>&1`; #command for merging read pairs
	#print "return code = ", $?, "\n" if $debug;
	if($?) {die "$mergeOut\n$!";} #sanity check
	print "\nMergining paired reads using FLASH with min overlap of $minOverlapToMerge NT for: ${$entryParameters[$entryCounter]}{SampleName}\n$mergeOut\nMerging done...\n" if $verbose;
	my $mergedFile="$outDir/${$entryParameters[$entryCounter]}{SampleName}.extendedFrags.fastq";
	if(-f $mergedFile)
		{return $mergedFile;} #return name of the file generated
	else
		{die "File not formed: $mergedFile\nCheck if any error in FLASH\n";} #integrity check
	}

sub CreateBWAIndex #creates indexes for the reference sequences provided and updates entry array with index prefix from sequence fileNames
	{
	my $indexDir=$_[0]; #location where indexes will be created
	my @refSeqs;
	foreach my $entryCounter(0..$#entryParameters)
		{
		push(@refSeqs,${$entryParameters[$entryCounter]}{ReferenceSeqFile}); #store refSeq fileName
		my $refIndexPrefix=(split(/\//,${$entryParameters[$entryCounter]}{ReferenceSeqFile}))[-1];
		$refIndexPrefix="$indexDir/$refIndexPrefix";
		${$entryParameters[$entryCounter]}{ReferenceSeqFile}=$refIndexPrefix; #update RefSeqFile element with BWA index prefix to be used later during alignment
		}
	my @uniqRefSeqs=uniq(@refSeqs);
	foreach my $uniqRefSeq(@uniqRefSeqs)
		{
		my $refIndexPrefix=(split(/\//,$uniqRefSeq))[-1];
		$refIndexPrefix="$indexDir/$refIndexPrefix";
		my $indexOut=`bwa index -p $refIndexPrefix -a is $uniqRefSeq 2>&1`; #command for making BWA index
		#print "return code = ", $?, "\n" if $debug;
		if($?) {die "$indexOut\n$!";} #sanity check
		print "\nPreparing BWA index for $uniqRefSeq\n$indexOut\nIndex Prepared...\n" if $verbose;
		};
	}

sub CreateBWAAlignment #creates alignment of reads to amplicon sequences, converts SAM->BAM and sorts BAM file
	{
	##Aligning reads to reference
	my $entryCounter=$_[0]; #this record number in the config file will be processed
	my $readFile=$_[1]; #sequences in this file will be aligned to corresponding reference
	my $outDir=$_[2]; #output files will be generated here
	my $alignFile="$outDir/${$entryParameters[$entryCounter]}{SampleName}.sam";
	my $errFile="$outDir/stderr.log";
	#bwa mem -t 4 $path/ReferenceIndexes/Barcode-01-ref.fa $path2/3.mergedReads/Barcode-01.extendedFrags.fastq >$path3/4.readAlignment/Barcode-01.sam
	`bwa mem -t 4 ${$entryParameters[$entryCounter]}{ReferenceSeqFile} $readFile >$alignFile 2>$errFile`; #command for generating alignment
	#print "return code = ", $?, "\n" if $debug;
	open(BWAERR,"$errFile") or die $!;
	my @stderrOut=<BWAERR>; #record STDERR for BWA
	close(BWAERR);
	unlink $errFile;
	if($?) {die "@stderrOut\n$!";} #sanity check
	print "\nCreating alignment of reads on amplicon refs using BWA for ${$entryParameters[$entryCounter]}{SampleName}\n@stderrOut\nAlignment finished...\n" if $verbose;
	
	##Converting SAM -> BAM
	my $alignBamFile=$alignFile;
	$alignBamFile=~s/\.sam$/.bam/;
	#print "\nBAM file: $alignBamFile\n" if $debug;
	#samtools view -bSh -o Barcode-01.bam Barcode-01.sam
	my $samtoolsOut=`samtools view -bSh -o $alignBamFile $alignFile 2>&1`; #command to convert SAM to BAM
	#print "return code = ", $?, "\n" if $debug;
	if($?) {die "$samtoolsOut\n$!";} #sanity check
	print "\nConverting alignment from SAM to BAM for ${$entryParameters[$entryCounter]}{SampleName}\n$samtoolsOut\nConversion finished...\n" if $verbose;
	
	##Sorting BAM
	my $alignSortedBamPrefix=$alignBamFile;
	$alignSortedBamPrefix=~s/\.bam/.Sorted/;
	#samtools sort Barcode-01.bam Barcode-01.Sorted
	my $bamSortOut=`samtools sort $alignBamFile $alignSortedBamPrefix 2>&1`; #command to sort BAM
	#print "return code = ", $?, "\n" if $debug;
	if($?) {die "$bamSortOut\n$!";} #sanity check
	print "\nSorting alignment by chromosomal coordinates in BAM for ${$entryParameters[$entryCounter]}{SampleName}\n$bamSortOut\nSorting finished...\n" if $verbose;
	
	my $alignSortedBamFile=$alignSortedBamPrefix.".bam";
	if(-f $alignSortedBamFile)
		{return $alignSortedBamFile;} #return name of the file generated
	else
		{die "File not formed: $alignSortedBamFile\nCheck if any error in BWA/Samtools\n";} #integrity check
	}

sub ScanIndels #scans indels in amplicons using the alignment files
	{
	my $entryCounter=$_[0]; #this record number in the config file will be processed
	my $alignFile=$_[1]; #alignment in this file will be used for indel calculation
	my $outDir=$_[2]; #output files will be generated here
	my $cutSiteRange=$_[3]; #indels in this range of the cut site will be calculated
	
	##Scan all cut sites
	my %cutSites;
	open(SITE,${$entryParameters[$entryCounter]}{AmpliconCutSites}) or die $!;
	while(<SITE>)
		{
		next if(/(#|^$)/); #skip comments/blank lines
		chomp();
		my @entry=split(/\t/);
		$cutSites{$entry[0]}=$entry[1];
		}
	close(SITE);
	
	##Create index of the alignment file
	unless(-e $alignFile.'bai')
		{
		my $samtoolsOut=`samtools index $alignFile 2>&1`; #command to create index of the alignment
		#print "return code = ", $?, "\n" if $debug;
		if($?) {die "$samtoolsOut\n$!";} #sanity check
		print "\nCreating BAM alignment index for ${$entryParameters[$entryCounter]}{SampleName}\n$samtoolsOut\nIndexing finished...\n" if $verbose;
		}
	}