#!/usr/bin/env nextflow
/*
vim: syntax=groovy

Pseudo-It-Nextflow : A NextFlow port of the Pseudo-it project by 
                     Brice Sarver ( https://github.com/bricesarver/pseudo-it ).
                     Pseudo-it is an approach that "iteratively generates 
                     pseudoreferences, incorporating sample-specific 
                     variation and reducing mapping biases."

Robert Hubley, 2017


Workflow Overview:
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
  
*/

////////////////////// LOCAL CONFIGURATION AND DEFAULTS /////////////////////////
// Default batch sizes for parallel jobs
params.bwaBatchSize = 5000000
params.haploBatchSize = 20000000

// The cluster to run on or "local"
params.cluster = "local"

// Java heap size in gb.
params.heapSize = 14

// Default output directory
params.outputDir = "results";

//
// Software dependencies
//
// The "module load java" isn't up to the task @TTU.  They have an older version
// on their clusters which has some bugs when used with the latest
// bioinformatics packages.
java = "/lustre/work/daray/software/jre1.8.0_151/bin/java -Xmx${params.heapSize}g -Djava.io.tmpdir=tmp"
picardJar = "/lustre/work/daray/software/picard-tools-2.15.0/picard.jar"
samtools = "/lustre/work/daray/software/samtools-1.3/samtools"
bwa = "/lustre/work/daray/software/bwa-0.7.12/bwa"
gatkJar = "/lustre/work/daray/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

//
// Setup executor for different environments
//
if ( params.cluster == 'local' ) {
  thisExecutor = 'local'
  thisQueue = ''
  thisOptions = ''
}else if ( params.cluster == 'hrothgar' ){
  // Hrothgar Cluster
  //   Yoda:  504GB,  1 node,  40 CPUs 
  //   R2D2:   32GB,  1 nodes, 20 CPUs
  //   Chewie: 16GB,  3 nodes, 20 CPUs each
  //                  3 nodes, 40 CPUs each
  //
  // Run NextFlow from: 
  //   qlogin -pe sm 10 -q Yoda -P communitycluster
  // 
  proc = 10
  //java = "/opt/apps/nfs/west/java/jre1.8.0_151/bin/java -Xmx${params.heapSize}g -Djava.io.tmpdir=tmp"
  thisExecutor = 'sge'
  thisQueue = 'Chewie,Yoda,R2D2'
  thisOptions = "-pe sm ${proc} -P communitycluster -S /bin/bash"
}else if ( params.cluster == 'quanah' ) {
  // Quanah Cluster
  //   omni:  188GB, 467 nodes, 36 CPUs each
  //
  // Run NextFlow from:
  //   qlogin -pe sm 24 -q omni -P quanah
  //
  proc = 12
  thisExecutor = 'sge'
  thisQueue = 'omni'
  thisOptions = "-pe sm ${proc} -P quanah -S /bin/bash"
}else {
  println "Unknown environment chosen!"
}

///////////////////////////////////////////////////////////////////////////////////

// Example dataset
params.reference = "${workflow.projectDir}/sample/600-ref.fa"
params.PE1 = "${workflow.projectDir}/sample/600_1.fastq"
params.PE2 = "${workflow.projectDir}/sample/600_2.fastq"
params.SE = "${workflow.projectDir}/sample/600_1.fastq"
if ( params.reference == "${workflow.projectDir}/sample/600-ref.fa" )
{
  // Lower default batch size....to have meaningful effect
  params.bwaBatchSize = 50000
  params.haploBatchSize = 30000
}


references = Channel.fromPath( params.reference )
refFile = file(params.reference)
refFileBase = refFile.getBaseName()
pe1File = file(params.PE1)
pe2File = file(params.PE2)

// Setup SE fastq batches if requested
if ( params.SE ) {
  seFile = file(params.SE)
  bwaSEBatchChan = Channel
               .from( seFile )
               .splitFastq( by: params.bwaBatchSize, file:true )
}

// Split PE fastq into batches for mapping
bwaPEBatchChan = Channel
               .from('foo', [ pe1File, pe2File ])
               .splitFastq( by: params.bwaBatchSize, pe:true, file:true )

// Say hello
log.info "\n"
log.info "Pseudo-It-NextFlow : Pseudoreference Generator- ver 0.4"
log.info "===================================================================="
log.info "Working Directory   : " + workflow.workDir
log.info "Ouptut Directory    : " + params.outputDir
log.info "Cluster             : " + params.cluster
log.info "Reference           : " + params.reference
log.info "PE1                 : " + params.PE1
log.info "PE2                 : " + params.PE2
log.info "SE                  : " + params.SE
log.info "bwaBatchSize        : " + params.bwaBatchSize + " reads"
log.info "haploBatchSize      : " + params.haploBatchSize + " records"
log.info "heapSize            : " + params.heapSize + " gb"
log.info "\n"


//////////////////////// W O R K F L O W ///////////////////////////////

process makeDictionary {
  input: 
  file refFile

  output:
  file "${refFileBase}.dict" into dictChan, dictChan2, dictChan3, dictChan4, dictChan5, dictChan6

  script:
  """
  ${java} -jar ${picardJar} CreateSequenceDictionary R=${refFile} O=${refFileBase}.dict
  """
}

process makeFaIdx {
  input: 
  file refFile

  output:
  file "${refFile}.fai" into faIdxChan, faIdxChan2, faIdxChan3, faIdxChan4, faIdxChan5

  script:
  """
  ${samtools} faidx ${refFile}
  """
}

process makeBWT {

  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions 

  input: 
  file refFile

  output:
  set file("${refFile}.bwt"), file("${refFile}.sa"), file("${refFile}.ann"), file("${refFile}.pac"), file("${refFile}.amb") into bwtChan, bwtChan2

  // Generates *.amb *.ann *.bwt *.pac *.sa

  script:
  """
  # Consider using -b 250000000 when ref is 2GB ( ie. refsize/8 )
  ${bwa} index ${refFile}
  """
}

process mapPE {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions 

  input:
  file foo from dictChan
  file bar from faIdxChan
  set file(bwtFile), file(saFile), file(annFile), file(pacFile), file(ambFile) from bwtChan
  file refFile
  set file(pe1File), file(pe2File) from bwaPEBatchChan

  output:
  file "${pe1File}${pe2File}.bam" into peBamChan
  
  script:
  """
  ${bwa} mem -M -t ${proc} ${refFile} ${pe1File} ${pe2File} | ${samtools} view -Sb - | ${samtools} sort - > ${pe1File}${pe2File}.bam 2> mapToBam.stderr
  """  
}

// Single Ended reads are optional to this workflow.  If specified
// slightly modify the workflow to include them
if ( params.SE ) {
  
  // Additional map SE reads
  process mapSE {
    executor = thisExecutor
    queue = thisQueue
    clusterOptions = thisOptions 

    input:
    file foo from dictChan2
    file bar from faIdxChan2
    set file(bwtFile), file(saFile), file(annFile), file(pacFile), file(ambFile) from bwtChan2
    file refFile
    file seFile from bwaSEBatchChan
  
    output:
    file "${seFile}.bam" into seBamChan
    
    script:
    """
    ${bwa} mem -M -t ${proc} ${refFile} ${seFile} | ${samtools} view -Sb - | ${samtools} sort - > ${seFile}.bam 2> mapToBam.stderr
    """  
  }

  // Merge both PE and SE bam files
  process mergePESEBamFiles {
    input:
    file(bamList) from peBamChan.toList()
    file(bamList) from seBamChan.toList()
  
    output:
    file 'merged.bam' into mergedBamChan
  
    script:
    """
    ${samtools} merge tmpMerge.out *.bam
    mv tmpMerge.out merged.bam
    """
  }

}else{

  // Just merge PE bam files
  process mergePEBamFiles {
    input:
    file(bamList) from peBamChan.toList()
  
    output:
    file 'merged.bam' into mergedBamChan
  
    script:
    """
    ${samtools} merge tmpMerge.out *.bam
    mv tmpMerge.out merged.bam
    """
  }
}

// For some reason the original pseudo-it sorted the bam after merging but only 
// in iteration 2 and above...?

process genReadGroups {
  input:
  file 'merged.bam' from mergedBamChan

  output:
  file 'readGrpd.bam' into readGrpdChan

  script:
  """
  ${java} -jar ${picardJar} AddOrReplaceReadGroups I=merged.bam O=readGrpd.bam SO=coordinate LB=spret_exome PL=illumina PU=misc SM=readGrpd VALIDATION_STRINGENCY=LENIENT 
  """
}

process dedupBAM {
  publishDir "${params.outputDir}", mode: 'copy'

  input:
  file 'readGrpd.bam' from readGrpdChan

  output:
  file 'dedup.bam' into dedupBamChan, dedupBamChan2

  script:
  """
  ${java} -jar ${picardJar} MarkDuplicates I=readGrpd.bam O=dedup.bam VALIDATION_STRINGENCY=LENIENT M=iteration1.dup_metrics
  """
}

process indexBAM {
  input:
  file 'dedup.bam' from dedupBamChan

  output:
  file "dedup.bam.bai" into idxdBamChan, idxdBamChan2

  script:
  """
  ${samtools} index dedup.bam
  """
}

process getBatches {

  input:
  file refFile
  file foo from idxdBamChan

  output:
  file "batch*.bed" into batchChan mode flatten

  script:
  """
  ${workflow.projectDir}/genBEDBatches.pl ${refFile} ${params.haploBatchSize}
  """
}

process HaplotypeCaller {
    executor = thisExecutor
    queue = thisQueue
    clusterOptions = thisOptions 

    input: 
    file refFile
    file dedupFile from dedupBamChan2
    file idxdFile from idxdBamChan2
    file dictFile from dictChan3
    file faIdx from faIdxChan3
    file batch_file from batchChan 

    output:
    file "${batch_file}.vcf" into vcfsChan
    
    script:
    """
    #
    # The original pseudo-it.py script used the now depreciated "-stand_emit_conf 10" 
    # option. Here we omit that but keep the "-stand_call_conf 30".
    #
    ${java} -jar ${gatkJar} -T HaplotypeCaller -R ${refFile} -I ${dedupFile} --genotyping_mode DISCOVERY -stand_call_conf 30 -o ${batch_file}.vcf -nct ${proc} -L ${batch_file} >& ${batch_file}.log
    """
}

process mergeVCFs {
  input:
  file reference_dict from dictChan6
  file(vcfList) from vcfsChan.toSortedList()

  output:
  file "concat.vcf" into mergedVCFChan

  script:
  """
  # Borrowed from CAW workflow:
  # first make a header from one of the VCF intervals
  # get rid of interval information only from the GATK command-line, but leave the rest
  FIRSTVCF=\$(ls batch*.vcf | sort -V | head -n 1)
  sed -n '/^[^#]/q;p' \$FIRSTVCF | \
  awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if(\$i!~/intervals=/ && \$i !~ /out=/){printf("%s ",\$i)}}printf("\\n")}' \
  > header
  (
    cat header
    for vcf in \$(ls batch*.vcf); do
      egrep -v "^#" \${vcf}
    done
  ) > unsorted-merge.vcf
  # Sort vcf file
  ${java} -jar ${picardJar} SortVcf I=unsorted-merge.vcf O=concat.vcf SEQUENCE_DICTIONARY=${reference_dict}
  """
}

process selectSNPs {
  publishDir "${params.outputDir}", mode: 'copy'
  input: 
  file refFile
  file faIdx from faIdxChan4
  file foo from dictChan4
  file vcfFile from mergedVCFChan
  
  output:
  file "filteredSNPs.vcf" into snpsChan

  script:
  """
  ${java} -jar ${gatkJar} -T SelectVariants -R ${refFile} -V ${vcfFile} -o ${vcfFile}.rawsnps --selectTypeToInclude SNP
  ${java} -jar ${gatkJar} -T VariantFiltration -R ${refFile} -V ${vcfFile}.rawsnps --filterName "mq30-5dp60" --filterExpression "MQ < 30.0 || DP < 5 || DP > 60" -o filteredSNPs.vcf 
  """
}

process callConsensus {
  publishDir "${params.outputDir}", mode: 'copy'

  input: 
  file "old_consensus.fa" from refFile
  file "old_consensus.fa.fai" from faIdxChan5
  file "old_consensus.dict" from dictChan5
  file snpFile from snpsChan
  
  output:
  file "consensus.fa" into result

  script:
  """
  ${java} -jar ${gatkJar} -T FastaAlternateReferenceMaker -R old_consensus.fa -o consensus.fa -V ${snpFile} 
  """
}

workflow.onComplete {
	    log.info """
	Pipeline execution summary
	---------------------------
	Completed at: ${workflow.complete}
	Duration    : ${workflow.duration}
	Success     : ${workflow.success}
	workDir     : ${workflow.workDir}
	exit status : ${workflow.exitStatus}
	Error report: ${workflow.errorReport ?: '-'}
	"""
	}
