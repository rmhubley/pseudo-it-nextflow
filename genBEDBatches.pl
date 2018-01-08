#!/usr/bin/env perl
use strict;
use Data::Dumper;

# Where this script can find liftUp, twoBitInfo and twoBitToFa 
my $ucscBinDir = "/lustre/work/daray/software";


my $fastaFile = $ARGV[0];
my $maxSize = $ARGV[1];


my $batches = getBatches( $fastaFile, $maxSize );

my $batchNum = 1;
#open LST,">batchFileList.txt" or die;
foreach my $batch ( @{$batches} ) {
    my $fileName = "batch-$batchNum.bed";
    open OUT,">$fileName" or die;
    foreach my $seq ( @{$batch} )
    {
      print OUT "".join("\t",@{$seq}) . "\n";
    }
    close OUT;
    #print LST "batch-$batchNum.bed\n";
    print "batch-$batchNum.bed\n";
    $batchNum++;
}
#close LST;

sub getBatches {
  my $fastaFile = shift;
  my $maxSeqSize = shift;
  my $numBatches = shift;
  
  ## Open up the fasta file and tabulate the sequence sizes
  # BUG: This doesn't work.  Strange thing to find in a Jim Kent project.
  #my $cmd = "$ucscBinDir/faToTwoBit $fastaFile stdout | $ucscBinDir/twoBitInfo stdin stdout |";
  #if ( ! -s "$fastaFile".".2bit" )
  #{
  #  system("faToTwoBit $fastaFile $fastaFile" . ".2bit");
  #}
  #my $cmd = "$ucscBinDir/twoBitInfo $fastaFile".".2bit stdout |";
  #open(P, $cmd)
  #  || die "Couldn't open pipe ($cmd): $_\n";
  my $totalSize = 0;
  my %seqSizes = ();
  my @batches = ();
  my $batchSize = 0;
  my @seqList = ();
  #while (<P>) {
  #  chomp;
  #  my ($seq, $seqSize) = split("\t");
  #  $totalSize += $seqSize;
  #  $seqSizes{$seq} = $seqSize;
  #}
  #close P;
  #
  open IN,"$fastaFile" or die;
  my $seq;
  my $id;
  while (<IN>){
    if ( /^>(\S+)/ ){
      my $tmpId = $1;
      if ( $seq ne "" ) 
      {
        $seqSizes{$id} = length($seq);
      }
      $seq = "";
      $id = $tmpId;
      next;
    }
    s/[\n\r\s]+//g;
    $seq .= $_;
  }
  if ( $seq )
  {
    $seqSizes{$id} = length($seq);
  }
  $seq = "";
  $id = "";
  close IN;

  if ( defined $numBatches && $numBatches > 0 )
  {
    $maxSeqSize = int($totalSize / $numBatches) + 1;
  }

  foreach my $seq ( keys(%seqSizes ) )
  {
    my $seqSize = $seqSizes{$seq};
    if ( ($batchSize + $seqSize) <= $maxSeqSize )
    {
      # Add it and move on 
      push @seqList, [ $seq, 0, $seqSize ];
      $batchSize += $seqSize;
      next;
    }

    # First chunk
    my $firstChunkSize = $maxSeqSize - $batchSize;
    push @seqList, [$seq, 0, $firstChunkSize];
    push @batches, [ @seqList ];
    $batchSize = 0;
    @seqList = ();


    # Middle chunks
    my $start = $firstChunkSize;
    for ( my $i = 0; $i < int(($seqSize-$firstChunkSize)/$maxSeqSize); $i++ )
    {
      push @seqList, [$seq, $start, $start + $maxSeqSize];
      push @batches, [ @seqList ];
      $batchSize = 0;
      @seqList = ();
      $start += $maxSeqSize;
    }
        
    # Last chunk
    if ( $start < $seqSize ) 
    {
       push @seqList, [$seq, $start, $seqSize];
       $batchSize += $seqSize - $start;
    }
  }
  if ( @seqList )
  {
    push @batches, [ @seqList ];
  }
  return \@batches;
} # getBatches

