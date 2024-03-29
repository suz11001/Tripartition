# Trippd (TRi-Partition based Pgt Detection)

Trippd is a prototype implementation of a simple proof-of-concept approach for detecting the presence of partial gene transfer (PGT) (i.e., horizontal transfer of a fragment of a gene) in a given gene family. trippd takes as input a multiple sequence alignment for the gene family under consideration, partitions the sites/columns of the alignment into three roughly equal parts, computes ML trees and bootstrap replicates for each partition, and compares these trees with each other to determine if that gene family has been affected by significant partial gene transfer. Further methodological details are described in the associated RECOMB-CG 2022 paper cited below. Trippd can be used to easily identify gene families whose gene trees may have been impacted by the presence of significant partial gene transfer.
 

## Dependencies
- Python 2.7.15
  - DendroPy
  - Scipy
  - Numpy
- RAxML 8.2.11


## Overview

Given a gene family multiple sequence alignment, partition the alignment in three approximately equal sized alignments, create maximum likelihood trees, create bootstrap trees, and conduct a statistical analysis for presence/absence of sub-gene/partial gene transfer event.

1. 01_splitFasta.py - Tripartition the multiple sequence alignment:  

   command: `python 01_splitFasta.py path/to/fasta/file`  
 
2. 02_raxml.sh - Generate maximum likelihood trees with RAxML:  
   
   command: `bash 02_raxml.sh /path/to/tripartitioned/alignment/directory/ model`

3. 03_raxml_bs.sh - Generate bootstrap tree (not rapid bootstraps) with RAxML:

   command: `bash 03_raxml_bs.sh /path/to/tripartitioned/alignment/directory/ model `

4. 04_hist_intersection.py - Conducts the statistical analysis for determining whether gene family is contains sub-gene transfer events:
   
   command: `python 04_hist_intersection.py --thresh 0.5 --bs_sample /path/to/tripartitioned/bootstrap/trees/ --ml_sample /path/to/tripartitioned/maximum/likelihood/trees/`

## Data Set Availability

Simulated data sets and the scripts used to generate such data are available in the sim_data directory.


## Citing trippd

Trippd can be cited as follows:

<a href="https://compbio.engr.uconn.edu/wp-content/uploads/sites/2447/2022/03/PartialGeneTransfer_RECOMBCG_2022.pdf">On Partial Gene Transfer and its Impact on Gene Tree Reconstruction</a><br>
Sumaira Zaman, Mukul S. Bansal<br>
RECOMB Comparative Genomics Conference (RECOMB-CG) 2022; LNBI 13234: 168–186
