#!/bin/bash -e

#  Processing the RNA-seq data 
# Files generated by this pipeline will be save in the directory specified in the config file (in the variable RESLOC)

usage(){
    NAME=$(basename $0)
    cat <<EOF
Usage:
  ${NAME}
  
Wrapper script for RNA-Seq analysis protocol used for the analysis of The Norfolk ADaPt trail's RNA-seq data.
This wrapper is based on the HISAT2/StringTie protocol (Pertea et al., Nat. Prot., 2016), 
with the adding of a QC step with TrimGalore! (Felix Krueger, The Babraham Institute).

This is a modified and extended version of the wrapper provided by the authors
of the protocol (to see the original visit ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol).

In order to configure the pipeline options (input/output files etc.)
please copy and edit a file rnaseq_pipeline-adapt.config.sh which must be
placed in the current (working) directory where this script is being launched.

Output directories "alignment-hisat2", "qc-trimgalore", "counts-prepDE" will be created in the 
current working directory or, if provided, in the given [output] (which will be created if it does not exist).


EOF
}


# if [[ "$1" ]]; then
 if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 1
 fi
# fi


## load variables
if [[ ! -f ./rnaseq-processing-adapt.config.sh ]]; then
 usage
 echo "Error: configuration file (rnaseq-processing-adapt.config.s) missing!"
 exit 1
fi

source ./rnaseq-processing-adapt.config.sh
WRKDIR=$(pwd -P)
errprog=""

if [[ ! -x $SAMTOOLS ]]; then
    errprog="samtools"
fi
if [[ ! -x $HISAT2 ]]; then
    errprog="hisat2"
fi
if [[ ! -x $STRINGTIE ]]; then
    errprog="stringtie"
fi
if [[ ! -x $TRIMGALORE ]]; then
    errprog="trimgalore"
fi
if [[ ! -x $GFFCOMPARE ]]; then
    errprog="gffcompare"
fi
if [[ ! -x $PREPDE ]]; then
    errprog="prepDE.py"
fi

if [[ "$errprog" ]]; then
  echo "ERROR: $errprog program not found, please edit the configuration script."
  exit 1
fi


#Create output directory
mkdir -p $RESLOC
cd $RESLOC


SCRIPTARGS="$@"
TRIMGALORELOC=./qc-trimgalore
ALIGNLOC=./alignment-hisat2
COUNTSLOC=./counts-prepDE

LOGFILE=./run.log

for d in "$TEMPLOC" "$TRIMGALORELOC" "$ALIGNLOC" "$QUANTIFKNOWNLOC" "$COUNTSKNOWNLOC" "$ASSEMBLYLOC" "$MERGELOC" "$MERGELOC" "$COMPARELOC" "$QUANTIFLOC" "$QUANTIFLOC" "$COUNTSLOC" ; do
 if [ ! -d $d ]; then
    mkdir -p $d
 fi
done

# main script block
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS
echo
echo "Found ${#reads1[@]} samples: "
echo "${reads1[@]%%.*}"

## QC using Trim Galore
echo
echo
echo "************ * Quality control and adapter trimming ( Trim Galore )"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    SAMPLE="${reads1[$i]%%.*}"
    SAMPLE="${SAMPLE%_*}"
    stime=`date +"%y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $SAMPLE"
    
    ## QC and adaptor sequence removal
    echo [`date +"%y-%m-%d %H:%M:%S"`] " starting qc ... "
    $TRIMGALORE -q $QUALITY --phred33 --stringency $STRINGENCY -length $LENGTH --paired -o $TRIMGALORELOC --fastqc ${FASTQLOC}/${SAMPLE}_${R1}.${FILETYPE}.gz ${FASTQLOC}/${SAMPLE}_${R2}.${FILETYPE}.gz
    echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"
done


## Alignment using HISAT2
echo 
echo 
echo "************ * Alignment of reads to genome (HISAT2)"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    SAMPLE="${reads1[$i]%%.*}"
    SAMPLE="${SAMPLE%_*}"
    stime=`date +"%y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $SAMPLE"

    echo [`date +"%y-%m-%d %H:%M:%S"`] " starting alignment ... "
    $HISAT2 -p $NUMCPUS -q --dta-cufflinks -x ${GENOMEIDX} \
      -1 ${TRIMGALORELOC}/${SAMPLE}_1_val_1.fq.gz \
      -2 ${TRIMGALORELOC}/${SAMPLE}_2_val_2.fq.gz \
      -S ${TEMPLOC}/${SAMPLE}.sam \
  	  --summary-file ${ALIGNLOC}/${SAMPLE}.txt --new-summary
	 
    
	echo [`date +"%y-%m-%d %H:%M:%S"`] " converting SAM to BAM ... "
    $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${SAMPLE}.bam ${TEMPLOC}/${SAMPLE}.sam

    echo [`date +"%y-%m-%d %H:%M:%S"`] " removing intermediate files"
    rm ${TEMPLOC}/${SAMPLE}.sam
    echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"
done




## Assembly using StringTie
echo 
echo 
echo "************ * Assembly of Transcripts (StringTie)"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    SAMPLE="${reads1[$i]%%.*}"
    SAMPLE="${SAMPLE%_*}"
    stime=`date +"%y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $SAMPLE"

    echo [`date +"%y-%m-%d %H:%M:%S"`] " starting assembly ... "
    $STRINGTIE -p $NUMCPUS -g ${GTFFILE} -o ${ASSEMBLYLOC}/${SAMPLE}.gtf \
      -l ${SAMPLE} ${ALIGNLOC}/${SAMPLE}.bam
    echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"

done



## Merge transcript files
echo
echo
echo "************ * Merge all Assembled Transcripts (StringTie)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
ls -1 ${ASSEMBLYLOC}/*.gtf > ${ASSEMBLYLOC}/mergelist.txt

$STRINGTIE --merge -p $NUMCPUS -g  ${GTFFILE} \
    -o ${MERGELOC}/stringtie_merged.gtf ${ASSEMBLYLOC}/mergelist.txt
echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"




## Compare Transcripts to the reference annotation
echo
echo
echo "************ * Compare Transcripts to the Reference Annotation (gffcompare)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
$GFFCOMPARE -r ${GTFFILE} -G -o ${COMPARELOC}/gffcompare_merged ${MERGELOC}/stringtie_merged.gtf
echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"


## Estimate Transcript Abundance
echo
echo
echo "************ * Estimate transcript abundance  (StringTie)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    SAMPLE="${reads1[$i]%%.*}"
    SAMPLE="${SAMPLE%_*}"
    stime=`date +"%y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $SAMPLE"
	
    if [ ! -d ${QUANTIFLOC}/${SAMPLE} ]; then
       mkdir -p ${QUANTIFLOC}/${SAMPLE}
    fi
    $STRINGTIE -e -B -p $NUMCPUS -G ${COMPARELOC}/gffcompare_merged.annotated.gtf \
      -o ${QUANTIFLOC}/${SAMPLE}/${SAMPLE}.gtf ${ALIGNLOC}/${SAMPLE}.bam \
      -A ${QUANTIFLOC}/${SAMPLE}/${SAMPLE}.tab

	echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"
done




## Create (unnormalised) counts table for differential expression with DESeq2 or edgeR
echo
echo
echo "************ * Estimate unnormalised count tables (prepDE.py)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
python $PREPDE -i $QUANTIFLOC -g ${COUNTSLOC}/gene_count_matrix.csv -t ${COUNTSLOC}/transcript_count_matrix.csv -l $READLENGTH
echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"



## Quantification (only known transcripts)
echo
echo
echo "************ * Estimate abundance for only known transcripts (StringTie)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    SAMPLE="${reads1[$i]%%.*}"
    SAMPLE="${SAMPLE%_*}"
    stime=`date +"%y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $SAMPLE"
	
	if [ ! -d ${QUANTIFKNOWNLOC}/${SAMPLE} ]; then
	   mkdir -p ${QUANTIFKNOWNLOC}/${SAMPLE}
	fi

	$STRINGTIE -e -B -p $NUMCPUS -G $GTFFILE \
		-o ${QUANTIFKNOWNLOC}/${SAMPLE}/${SAMPLE}.gtf ${ALIGNLOC}/${SAMPLE}.bam \
		-A ${QUANTIFKNOWNLOC}/${SAMPLE}/${SAMPLE}.tab

	echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"
done


## Create (unnormalised) counts table for differential expression with DESeq2 or edgeR (only known genes/transcripts)
echo
echo
echo "************ * Estimate unnormalised count tables for only known transcripts (prepDE.py)"
echo [`date +"%y-%m-%d %H:%M:%S"`] " starting ..."
python $PREPDE -i $QUANTIFKNOWNLOC -g ${COUNTSKNOWNLOC}/gene_count_matrix.csv -t ${COUNTSKNOWNLOC}/transcript_count_matrix.csv -l $READLENGTH
echo [`date +"%y-%m-%d %H:%M:%S"`] " ... done"




echo
echo
echo [`date +"%y-%m-%d %H:%M:%S"`] "#> DONE!!  End of Processing"




} #pipeline end


pipeline 2>&1 | tee $LOGFILE
