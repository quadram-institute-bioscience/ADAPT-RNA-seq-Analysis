## Configuration file to process RNA-seq raw data using the wrapper: rnaseq-pipeline.sh

## This file has been adapted from: ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/rnaseq_pipeline.config.sh
##
## Place this script in a working directory and edit it accordingly.
##
## The default configuration assumes that the user has downloaded the raw 
## data in the data/raw-data directory.
##
##
## File paths for input data
##  Full absolute paths are strongly recommended here.
##  Warning: if using relatives paths here, these will be interpreted 
##  relative to the  chosen output directory (which is specified in the
##  RESLOC parameter)
##
##
## Alignment was done using hisat2 indexes based on the Homo sapiens GRCh38.97 genome. It's assumed the file is in the local directory (it can be changed editing variable GENOMEIDX): 
## 			<PROJECT-DIR>/data/resources/HISAT2/indexes/homo-sapiens/GRCh38.97/genome
##
##
## Reference annotation was used from ensemble: GRCh38, release 97 (file:  Homo_sapiens.GRCh38.97.gtf). Annotation can be downloaded directly from ensembl (http://www.ensembl.org/). This script assumes the gtf file is saved to local directory (it can be changed editing variable GTFFILE):
## 			<PROJECT-DIR>/data/resources/ensembl/97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf



#how many CPUs to use on the current machine?
NUMCPUS=4

#  Base directory, if most of the input files have a common path
BASEDIR=".."    ## path the project's directory
# Results Directory
RESLOC=${BASEDIR}/results
# Log Files
LOGS=${RESLOC}/log-files


############################################
###  PROJECT CONFIGURATION   ###
############################################
#Read lenght
READLENGTH=151

FASTQLOC="$BASEDIR/data/raw-data"          # location of the raw RNA-seq data
R1="1"
R2="2"
FILETYPE="fastq"
## For temporary files use the following
TEMPLOC="$BASEDIR/tmp"


############################################
### PROGRAMS CONFIGURATION  ###
############################################

# Path to shared databases, genomes and annotation
RESOURCES="$BASEDIR/data/resources"

# If these programs are not in any PATH directories, please edit accordingly:
HISAT2=$(which hisat2)
STRINGTIE=$(which stringtie)
SAMTOOLS=$(which samtools)
TRIMGALORE=$(which trim_galore)
GFFCOMPARE=$(which gffcompare)
PREPDE=$(which prepDE.py)


############################################
###   QC: TRIMGALORE     ###
############################################
# Trim Galore
QUALITY=30
STRINGENCY=5
LENGTH=60
TRIMGALORELOC="$RESLOC/qc-trimgalore"


############################################
###   ALIGNMENT: HISAT2   ###
############################################
# HISAT2 outputs
ALIGNLOC=${RESLOC}/alignment-hisat2

## Reference genome, annotation (ensembl) file
## HISAT2 INDEXES and annotation files
#  Human, GRCh38.97
GENOMEIDX=${RESOURCES}/HISAT2/indexes/homo-sapiens/GRCh38.97/genome
GTFFILE=${RESOURCES}/ensembl/97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf


##################################################
### ASSEMBLY ONLY KNOWN TRANSCRIPTS: StringTie ###
##################################################
QUANTIFKNOWNLOC=${RESLOC}/quantification-known-transcripts-stringtie
COUNTSKNOWNLOC=${RESLOC}/counts-known-transcripts-prepDE


############################################
###  ASSEMBLY: StringTie   ###
############################################
ASSEMBLYLOC=${RESLOC}/assembly-stringtie


############################################
###  MERGE TRANSCRIPTS: StringTie   ###
############################################
MERGELOC=${RESLOC}/merge-transcript-stringtie
COMPARELOC=${RESLOC}/compare-to-reference-gffcompare


############################################
###  QUANTIFY TRANSCRIPTS: StringTie   ###
############################################
QUANTIFLOC=${RESLOC}/quantification-stringtie
COUNTSLOC=${RESLOC}/counts-prepDE



############################################
## List of Samples ##
############################################
reads1=(${FASTQLOC}/*_1.fastq.gz)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1.fastq.gz/_2.fastq.gz}")