# ADAPT-RNA-seq-Analysis
Processing of RNA-seq data and Differential Expression Analysis for The Norfolk ADaPt trial.



## Overview

The processing of raw RNA-sequencing data follows the protocol by Pertea et al. [1]. The main steps are as follows:

* QC, removal of adaptor sequences and low quality reads using Trim Galore

* Alignment of reads to the reference genome using HISAT2

* Assembly of the alignments into full-length transcripts using StringTie

* Quantification of the expression levels using StringTie

* Calculate differential expression using DESeq2

* Gene Set Enrichment Analysis (not included here)


## Requirements

- [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (>= 0.6.5)
- [HISAT2](https://daehwankimlab.github.io/hisat2/) (>= 2.1.0)
- [StringTie](https://ccb.jhu.edu/software/stringtie/) (>= 1.3.5)
- [Samtools](http://www.htslib.org/) (>= 1.4.1)
- [GFFCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) (>= 0.9.8)



## Processing RNA-seq data

To start the processing, move to the source code directory (src/), and run the wrapper script:

```
$ cd src/
$ ./rnaseq-processing.sh
```



## Differential Expression Analysis with DESeq2

Differential expression analysis is performed using DESeq2 (src/DE-analysis.Rmd). Defaults assume that the gene counts and metadata are saved as csv files in the data directory.

NOTE: RNA-seq data and processed data will be made available upon publication.



#### Contact

For queries contact Perla (dot) Rey (at) quadram (dot) ac (dot) uk


### References

[1] Pertea M., Kim D., Pertea G.M., Leek J.T., Salzberg S.L., Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown, Nature Protocols, 2016.
