#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) pseudo-it-nextflow.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A script to iteratively call the NextFlow pseudo-it-nextflow.nf
##      workflow.  The NextFlow framework is a powerful framework for
##      describing dependencies in a serial and parallel fashion.  
##      Unfortunately it does not handle iteration very well. This
##      script simply sets up the data files for a NextFlow run calls
##      the workflow in an iterative fashion.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
=head1 NAME

pseudo-it-nextflow.pl - Iteratively run the pseudo-it-nextflow workflow

=head1 SYNOPSIS

  pseudo-it-nextflow.pl [-reference <fasta>] -PE1 <paired_end_R1.fastq>
                        -PE2 <paired_end_R2.fastq> [-SE <unpaired.fastq>]
                        [-iteration <#>] -outputDir <directory_for_resulsts>

 Run 7 iterations of the pseudo-it-nextflow.nf workflow and save the results
 of each round to "-outputDir". Optionally allow recovery of a failed run
 by specifying the "-iteration" parameter and omiting the "-reference" 
 parameter ( in this case reference is assumed to be the output of the
 previous (iteration-1) consensus file output.

 For example:

 Hrothgar:
   ssh hrothgar
   qlogin -pe sm 10 -q Yoda -P communitycluster
   cd /lustre/scratch/username
   nohup /lustre/work/daray/software/pseudo-it-nextflow-v0.4/pseudo-it-nextflow.pl \
        -reference /lustre/work/foobar/foo.fa \
        -PE1 /lustre/work/foobar/R1.fastq \
        -PE2 /lustre/work/foobar/R2.fastq \
        -SE  /lustre/work/foobar/SE.fastq \
        -cluster hrothgar \
        -outputDir /lustre/scratch/username/results-genome1 >& run-genome1.log &

                             

 Quanah:
   ssh quanah
   qlogin -pe sm 12 -q omni -P quanah
   cd /lustre/scratch/username
   nohup /lustre/work/daray/software/pseudo-it-nextflow-v0.4/pseudo-it-nextflow.pl \
        -reference /lustre/work/foobar/foo.fa \
        -PE1 /lustre/work/foobar/R1.fastq \
        -PE2 /lustre/work/foobar/R2.fastq \
        -SE  /lustre/work/foobar/SE.fastq \
        -cluster quanah \
        -outputDir /lustre/scratch/username/results-genome1 >& run-genome1.log &


  If a run fails you will see an indication of that in the *.log file.  Also
  NextFlow should "qdel" all the running or queued jobs automatically.  If it
  looks like it was "Killed" rather than error'd out due to an data problem,
  then you should just restart the run using the previous iterations output.
  To do this simply alter the command line by removing the "-reference" 
  parameters and including a new "-iteration #" parameter.  The iteration
  number is the one you want to re-run.  I.e if it died while working iteration
  5 simply start:

   nohup /lustre/work/daray/software/pseudo-it-nextflow-v0.4/pseudo-it-nextflow.pl \
        -PE1 /lustre/work/foobar/R1.fastq \
        -PE2 /lustre/work/foobar/R2.fastq \
        -SE  /lustre/work/foobar/SE.fastq \
        -cluster quanah \
        -iteration 5 \
        -outputDir /lustre/scratch/username/results-genome1 >& run-genome1.log &



=head1 DESCRIPTION

A script to iteratively run the NextFlow port of the Pseudo-it project.  
The Pseudo-it project, originally written by Brice Sarver 
( https://github.com/bricesarver/pseudo-it ), is an approach that 
"iteratively generates pseudoreferences, incorporating sample-specific 
variation and reducing mapping biases.". 

By porting it to NextFlow aspects of the workflow can be massively 
parallelized ( particularlly index building, mapping, and variant 
calling ) while other aspects of the flow can remain serial.  In 
addition the parallel aspects may be run on a single node with many
processors or on many nodes through a job-scheduler.

=head1 Workflow Overview

  [serial]   Split FASTQ files into batches
  [parallel] Index starting reference FASTA file 
              - make dictionary
              - make faIdx
              - make BWT
  [parallel] Map Read Batches ( bwa - both paired and optional single ended reads )
               - Added sort to each job to speed up merging
  [serial]   Merge BAM files
  [serial]   Generate read groups ( picard AddOrReplaceReadGroups )
  [serial]   Dedup ( picard MarkDuplicates )
  [serial]   Index ( samtools index )
              - No need to run gatk RealignerTargetCreator ( only needed by 
                UnifiedGenotyper )
              - No need to run gatk IndelRealigner ( only needed by 
                UnifiedGenotyper )
  [serial]   Create BED file batches for HaplotypeCaller
  [parallel] Variant calling ( gatk HaplotypeCaller )
               - Latest version of gatk doesn't support the "-stand_emit_conf 10"
                 option anymore. From my reading of the docs it may be that this
                 was overriden by the "-stand_call_conf 30" parameter anyway.
  [serial]   Merge VCF files
  [serial]   Select SNPs ( gatk SelectVariants, gatk VariantFiltration )
  [serial]   Call consensus ( gatk FastaAlternateReferenceMaker )

=head1 NextFlow Environment

NextFlow maintains an execution state database that needs to be accessible
to the node running the workflow.  In addition this database needs to
stored on a filesystem that support file locking.  Since /lustre does not
support file locking I choose to create a wrapper script for NextFlow to
set the NextFlow database and run directory to the /tmp/<username> directory.
If the directory doesn't already exist it's created before NextFlow is started
up.  In this way a locking file system is guaranteed.  Unfortunately this means
that it's important to remember which system NextFlow was initiated on in order
to use features like "-resume". The NextFlow working directories on the other
hand need to be located on a shared filesytem ( file locking not necessary ).
The wrapper currently sets this to the NXF_WORK subdirectory in the current
working directory.  The protocol for running on Hrothgar would then be:

   # Login to the cluster head node
   ssh hrothgar

   # Pick an adequately large queue to run the NextFlow environment from
   qlogin -pe sm 10 -q Yoda -P communitycluster

   # Change to a shared filesystem to run from
   cd /lustre/scratch/myworkdir

   # Run NextFlow wrapper on a workflow
   /lustre/work/daray/software/Nextflow/workflow run my-workflow.nf 

   # Now...you can finde a new NXF_WORK directory containinig all the
   # NextFlow intermediate files:
   ls NXF_WORK

   # And you can look at the NextFlow database and report files:
   ls /tmp/myusername


The options are:

=over 4

=item -reference

The starting reference to align reads to

=item -bwaBatchSize <#>

Optional parameter to control the size of the mapping batches.  Currently
set to 5000000.  By increasing this, jobs will run slower on the cluster ( typically
not a good way to batch ) but the number of jobs submitted to the queue will
be smaller  This is important if the queue admins have a limit on the number
of jobs that can be submitted at once.

=item -haploBatchSize <#>

Optional parameter to control the size of the HaplotypeCaller batches.  Currently
defaults to 20000000.

=back

=head1 DEPENDENCIES

Java 1.8u155 JRE, GATK, Picard, Samtools, BWA

=head1 SEE ALSO
 
=head1 COPYRIGHT

Copyright 2017 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

# 
# Modules
#
use strict;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use warnings;
use strict;
use File::Spec;
use FindBin;


#
# Where NextFlow and the pseudo-it-nextflow.nf file exist:
#
my $NextFlow = "/lustre/work/daray/software/Nextflow/nextflow";
my $pseudoItWorkflow = "$FindBin::RealBin/pseudo-it-nextflow.nf";
my $Version = "0.1";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    "-reference=s",
    "-PE1=s",
    "-PE2=s",
    "-SE=s",
    "-name=s",
    "-cluster=s",
    "-iteration=i",
    "-outputDir=s",
    "-bwaBatchSize=i",
    "-haploBatchSize=i"
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}


sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0|less";
  exit;
}

if ( ! keys( %options ) ||
     ( ! ( $options{'PE1'} && $options{'PE2'} && 
           $options{'cluster'} && $options{'outputDir'} ) ) )
{
  usage;
}

my $iteration = 1;
$iteration = $options{'iteration'} if ( exists $options{'iteration'} );

my $username = getpwuid($<);
my $outputDir = $options{'outputDir'};

mkdir($outputDir) if ( ! -d $outputDir);
if ( ! -d "$outputDir/iter-0" && exists $options{'reference'} )
{
  mkdir("$outputDir/iter-0");
  system("ln -s $options{'reference'} $outputDir/iter-0/consensus.fa" );
}
mkdir( "$outputDir/iter-$iteration" ) if ( ! -d "$outputDir/iter-$iteration" );

if ( ! exists $options{'bwaBatchSize'} )
{
  $options{'bwaBatchSize'} = 5000000;
}

if ( ! exists $options{'haploBatchSize'} )
{
  $options{'haploBatchSize'} = 20000000;
}

while ( $iteration < 8 && $iteration > 0 )
{
  $options{'outputDir'} = $outputDir . "/iter-$iteration";
  $options{'reference'} = "$outputDir/iter-" . ( $iteration-1 ) . "/consensus.fa";
  
  my $optionStr = "";
  foreach my $option ( keys %options ) 
  {
    $optionStr .= " --$option " . $options{$option};
  }

  my $reportFile = "nf-$$-iter-$iteration-report";
  my $timelineFile = "nf-$$-iter-$iteration-timeline";
  my $dagFile = "nf-$$-iter-$iteration-dag";
  my $iterCmd = "$NextFlow run $pseudoItWorkflow $optionStr -with-report $reportFile -with-timeline $timelineFile -with-dag $dagFile ";
  
  # Run
  print "Running iteration $iteration\n";
  print "cmd: $iterCmd\n";
  system($iterCmd);
  
  if ( -s "$outputDir/iter-$iteration/consensus.fa" &&
       -s "$outputDir/iter-$iteration/filteredSNPs.vcf" &&
       -s "/tmp/$username/$reportFile" )
  {
    # Copy report and trace.txt files to shared folders from the /tmp/user file system
    system("mv /tmp/$username/$reportFile $outputDir/iter-$iteration");
    system("mv /tmp/$username/$timelineFile $outputDir/iter-$iteration");
    system("mv /tmp/$username/$dagFile $outputDir/iter-$iteration");
    $iteration++;
  }else {
    $iteration = 0;
  }
}

1;
