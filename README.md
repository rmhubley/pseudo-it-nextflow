
### Pseudo-It-Nextflow : A NextFlow port of the Pseudo-it project

Pseudo-it is a method developed by Brice Sarver( https://github.com/bricesarver/pseudo-it ).
which generates a pseudoreference through a process of iterative mapping and reference calling
steps.  This project was designed to run on the TTU clusters but could easily be generalized 
for other environments.

Run pseudo-it-nextflow.pl without any parameters to see documentation.

TL;DR

  ssh hrothgar
  qlogin -pe sm 10 -q Yoda -P communitycluster
  cd /lustre/scratch/my_project
  nohup /lustre/work/daray/software/pseudo-it-nextflow-v0.4/pseudo-it-nextflow.pl \
        -reference /lustre/work/foobar/foo.fa \
        -PE1 /lustre/work/foobar/R1.fastq \
        -PE2 /lustre/work/foobar/R2.fastq \
        -SE  /lustre/work/foobar/SE.fastq \
        -cluster hrothgar \
        -outputDir /lustre/scratch/username/results-genome1 >& run-genome1.log &


Release Notes:

0.4 : Initial release
0.5 : Added bwaBatchSize and haploBatchSize as command line parameters.


Robert Hubley, 2017
