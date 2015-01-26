# hiHMM

**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [hiHMM](#)
	- [About](#)
	- [An Instruction of how to use hiHMM as an example](#)
		- [Step 1: Create hiHMM input files](#)
		- [Step 2: Remove unmappable regions using hiHMM_Pre1.R](#)
		- [Step 3: Run hiHMM using driv_hihmm.m](#)
		- [Step 4: Annotate states in the emission matrix](#)
		- [Step 5: Plot the emission and transition matrices and recolour the output bed files using](#)
		- [Step 6: Reintroduce unmappable regions as a new state using hiHMM_Post2.r](#)
		- [Results](#)
	- [References](#)
		- [Download the raw data](#)
		- [Fly and worm files](#)
	- [License](#)
	
## About

**hiHMM (hierarchically-linked infinite hidden Markov model)** is a new Bayesian non-parametric method to jointly infer chromatin state maps in multiple genomes (different cell types, developmental stages, even multiple species) using genome-wide histone modification data. 

## An Instruction of how to use hiHMM as an example

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

### Step 1: Create hiHMM input files

There should be one file for each condition and each chromosome, and should be named as `condition_chromosome.txt`. Each file should contain all the ChIP-seq tracks (e.g. histone modifications, transcription factors etc).

The sample files provided contain the signal values for 200 bp bins and in the following structure: `fly_chr2L.txt`.

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig1.png)

These raw files should be stored in `sample_analysis/Unmapped_input_files`.

### Step 2: Remove unmappable regions using hiHMM_Pre1.R

Unmappable regions are removed from the hiHMM input files using the information contained in the `sample_analysis/Mappability_files` folder, which contains the **mappable** regions for each chromosome for each condition. These files are .bed files in the format `chromosome\tstart position\tend position`. For example, for fly: `mappable.dm3.lucy.bed`.

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig2.png)

The R script then removes all regions that are not mappable, for example the input file for **fly chr2L** now
begins at **5100 bp**, since the region **1 bp ? 4991 bp** is unmappable according to the mappability file: `fly_chr2L.txt`.

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig3.png)

To run this script, users need to specify:
* The working directory to be the analysis folder
* `unmapped_file_dir`: the location of the hiHMM unmapped input files
* `mapped_file_dir`: the location where the mapped files should be written to
* `source_dir`: the location of the Mappability folder
* `Mappable`: the mappable files and a species/condition ID, for example **fly** and **worm**

The files with the unmappable regions removed are written to `sample_analysis/Mapped_input_files`.

### Step 3: Run hiHMM using driv_hihmm.m

The analysis folder, input and output folders need to be specified. Parameters such as the conditions,
mapped file names, chromosome labels, <u>model number</u>, <u>bin size</u> (same as the input files e.g. 200 bp), etc, also need to be specified.

In the example given, the results from hiHMM, including the chromatin annotations (bed files e.g.
`hihmm.model2.K7.fly.bed` for fly) and emission/transition matrices, will be written to the
`sample_analysis/output` folder.

### Step 4: Annotate states in the emission matrix

The emission matrix (`sample_analysis/output/train-hihmm-model2.emission.csv`) will contain the **K**
chromatin states along the rows and the ChIP-seq tracks used along the columns.

Note that for Model 1,the output file will contain all the emission matrices for all conditions. For Model 2, there will only be one emission matrix that has been jointly inferred across all conditions.

The states need to be functionally annotated (e.g. promoter, enhancer) in the following format, and
separated by spaces: 

```
State_Number State_Name Species/Condition
```

There are a number of state names and colourings that have been pre-defined in the colourise
function in the `hiHMM_Post1.R` script, including promoter, enhancer, gene, transcription, repressed,
heterochromatin, and low signal. *Note that the last element is optional for Model 2, but should be
included in Model 1 as a sample identifier.* See `hiHMM_Post1.R` for more details on naming states. 

An
example of the named emission matrices are given below:
* `train-hihmm-model1.emission_named.csv`: Note that the first emission matrix is for **Worm (W)** and the
second for **Fly (F)**, as specified by the order in the **condition** parameter in `driv_hihmm.m`.

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig4.png)

* `train-hihmm-model2.emission_named.csv`: A single set of states jointly inferred

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig5.png)

**Important!**: The emission matrix with annotated states needs to be saved with the `_named` suffix, i.e.
`train-hihmm-model2.emission_named.csv`.

### Step 5: Plot the emission and transition matrices and recolour the output bed files using
hiHMM_Post1.R

The emission and transition matrices can be plotted as PDF files using this script. For Model 1, an
emission matrix will be plotted for each condition, while only a single emission matrix will be plotted for
Model 2. A transition matrix will be plotted for each condition regardless of which model is used. These
will be stored in the same `sample_analysis/output` folder.

To run this script, users need to specify:
* The working directory to be the hiHMM output folder (subfolder in the analysis folder)
* **m**: the model number (1 or 2)
* **c**: the number of conditions or samples. In this example the number of samples is 2, for fly and worm
* **fworder**: ensure that the ChIP-seq tracks (e.g. histone modifications) used in the analysis are present in this list and in the desired order
* **colourise**: ensure all annotated states are accounted for in this function with the desired colours

Examples of plotted emission and transition matrices for the two models are provided below:

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig6.png)

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig7.png)

The chromatin state segments in the .bed files will also be recoloured according to the same colouring pattern in the emission matrix. These files will have the **recoloured** suffix, for example, `hihmm.model2.K7.fly.bed` will become `hihmm.model2.K7.fly.recoloured.bed` and stored in the same `sample_analysis/output` folder.

### Step 6: Reintroduce unmappable regions as a new state using hiHMM_Post2.r
The unmappable regions that were removed in Step 2 prior to running hiHMM will now be added back to the `.bed` files as **State 0**.
To run this script, users need to specify:
* the working directory to be the hiHMM output folder (subfolder in the analysis folder)
* **outdir**: the output directory where the remapped files will be written to
* **bin_size**: bin size (same as the hiHMM input files)
* **mappability_dir**: the location of the Mappability folder
* **Mappable**: the mappable files along with a species/condition ID, for example “fly” and “worm”
* **chr_lengths_dir**: the Chromosome Lengths folder
* **chr_lengths**: the chromosome lengths files along with the same species/condition ID i.e. “fly” and “worm”
The output files will have the `ReMapped` suffix, so for example, `hihmm.model2.K7.fly.recoloured.bed` will become `hihmm.model2.K7.fly.recoloured.ReMapped.bed` and written to the `sample_analysis/output/ReMapped` folder.

### Results

The following image shows a screenshot of the IGV Genome Browser for fly showing the different types of `.bed` files that are produced from Steps 3, 5 and 6 respectively.

![](https://dl.dropboxusercontent.com/u/5327300/hihmm/instruction/fig8.png)

## References

### Download the raw data
Here is a link to the EncodeX browser where you can download the normalized ChIP-Seq data that was used in our analysis, and much more!

[EncodeX browser](http://encode-x.med.harvard.edu/data_sets/chromatin)

### Fly and worm files

We also provide the auxillary files you will need to recreate our analysis, including chromosome lengths and unmappable regions.
[Download aux. files](https://sites.google.com/site/kasohn/files/fly_worm_aux.zip?attredirects=0&d=1)

## License

```
The MIT License (MIT)

Copyright (c) 2015 Kyung-Ah Sohn 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

<div style="display:block;margin-right:auto;text-align:right">
<a href="http://www.ajou.ac.kr/en/index.jsp" rel="nofollow"><img border="0" height="35" src="https://sites.google.com/site/kasohn/_/rsrc/1391695588708/software/ajou3.gif?height=35" /></a>
<a href="https://research.unsw.edu.au/people/dr-joshua-ho" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696102021/software/unsw2.png?height=30" /></a>
<a href="http://www.victorchang.edu.au/home/our-research/faculty-detail/?faculty_name=dr-joshua-w-k-ho" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696055263/software/vccri.png?height=30" /></a>
<a href="http://www.brighamandwomens.org/Departments_and_Services/medicine/services/genetics/default.aspx" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696035753/software/brig1.jpg?height=30" /></a>
<a href="http://compbio.med.harvard.edu" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696008781/software/harvard2.png?height=30" /></a>
<a href="http://www.snubi.org" rel="nofollow"><img border="0" height="30" src="https://sites.google.com/site/kasohn/_/rsrc/1391696120344/software/seoul3.png?height=30" />
</a></div>
