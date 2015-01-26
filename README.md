# Introduction

**hiHMM (hierarchically-linked infinite hidden Markov model)** is a new Bayesian non-parametric method to jointly infer chromatin state maps in multiple genomes (different cell types, developmental stages, even multiple species) using genome-wide histone modification data. 

# Instruction

## Directory Structure

Each analysis or experiment should have its input and output files contained within a single folder. Note that the chromosome names in the input files should be the same as those in the mappability files and chromosome lengths files. A separate output folder for each model can also be used instead, e.g. `output_model1` and `output_model2`.

* hiHMM
 * `driv_hiHMM.m` (script to run hiHMM)
 * R_scripts
  * sample_analysis
    * Chromosome_Lengths
    * Mappability_files
    * Mapped_input_files (output from `hiHMM_Pre1.R` in **Step 2**)
    * Unmapped_input_files
    * Output (output from `driv_hiHMM.m` in **Step 3**)

The R scripts use a number of packages, which are automatically detected and installed by the scripts. These include:
* `hiHMM_Pre1.r`: [IRanges](http://bioconductor.org/packages/2.3/bioc/html/IRanges.html)
* `hiHMM_Post1.r`: [gplots](http://cran.r-project.org/web/packages/gplots/index.html), [RColorBrewer](http://cran.r-project.org/web/packages/RColorBrewer/index.html)
* `hiHMM_Post2.r`: [IRanges](http://bioconductor.org/packages/2.3/bioc/html/IRanges.html)

## Step 1: Create hiHMM input files

There should be one file for each condition and each chromosome, and should be named as `condition_chromosome.txt`. Each file should contain all the ChIP-seq tracks (e.g. histone modifications, transcription factors etc).

The sample files provided contain the signal values for 200 bp bins and in the following structure:
`fly_chr2L.txt`

These raw files should be stored in `sample_analysis/Unmapped_input_files`.


<div style="display:block;margin-right:auto;text-align:right">
<a href="http://www.ajou.ac.kr/en/index.jsp" rel="nofollow"><img border="0" height="35" src="https://sites.google.com/site/kasohn/_/rsrc/1391695588708/software/ajou3.gif?height=35" /></a>
<a href="https://research.unsw.edu.au/people/dr-joshua-ho" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696102021/software/unsw2.png?height=30" /></a>
<a href="http://www.victorchang.edu.au/home/our-research/faculty-detail/?faculty_name=dr-joshua-w-k-ho" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696055263/software/vccri.png?height=30" /></a>
<a href="http://www.brighamandwomens.org/Departments_and_Services/medicine/services/genetics/default.aspx" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696035753/software/brig1.jpg?height=30" /></a>
<a href="http://compbio.med.harvard.edu" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696008781/software/harvard2.png?height=30" /></a>
<a href="http://www.snubi.org" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696120344/software/seoul3.png?height=30" />
</a></div>
