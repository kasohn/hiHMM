#####################################################################################################
### Mappability script 1 of 2 - hiHMM: Bayesian non-parametric joint inference of chromatin state maps (Sohn et al. 2014)
### Djordje Djordjevic - d.djordjevic@victorchang.edu.au - 2013 
###
### This file contains a script to remove unmappable regions from your hiHMM formatted input files 
###
### In order to make this script generalisable you must follow a simple naming convention:
### Assumes that the hiHMM input file naming convention is "sample ID"_"Chromosome"
### Assumes each hiHMM input file contains data for ONLY ONE chromosome
### Make sure that "sample ID" is not nested / repeated, as in do not use FLY and FLY2
### 
#####################################################################################################

#####################################################################################################
### START OF SCRIPT TO REMOVE UNMAPPABLE REGIONS FROM HIHMM INPUT FILES

### This script uses the IRanges library.
### If it is not installed, then install it.
if (!require("IRanges")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("IRanges")
  require("IRanges")
}

### Set the working directory to the analysis folder.
setwd("../sample_analysis/")

### Set the directory where your hiHMM formatted input files are.
unmapped_file_dir <- "./Unmapped_input_files/"
### Set the directory where the output files will be written to.
### Note that this will be the location of the input file directory used by hiHMM (specified by the "indir" parameter in driv_hihmm.m)
mapped_file_dir <- "./Mapped_input_files/"
dir.create(mapped_file_dir)

### Set the location of the mappability files.
source_dir <- "./Mappability_files/"

Mappable <- list()

### Load mappability files in bed format (chromosome \t start position \t end position)
### 1 file per sample ID
### Mappability files contains regions that ARE mappable
### Load the mappability file corresponding to each sample ID into the Mappable list
### ie. Mappable[["Sample ID"]] <- the bed file

Mappable[["fly"]] <- read.table(paste(source_dir, "mappable.dm3.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
Mappable[["worm"]] <- read.table(paste(source_dir, "mappable.ce10.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")

### Example for different developmental stages across different organisms
#Mappable[["FlyL3"]] <- read.table(paste(source_dir, "mappable.dm3.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
#Mappable[["FlyEL"]] <- read.table(paste(source_dir, "mappable.dm3.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
#Mappable[["FlyAH"]] <- read.table(paste(source_dir, "mappable.dm3.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
#Mappable[["WormAH"]] <- read.table(paste(source_dir, "mappable.ce10.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
#Mappable[["WormL3"]] <- read.table(paste(source_dir, "mappable.ce10.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")

files <- list()

### Load the hiHMM input files. The files MUST follow the naming convention "sample ID"_"chromosome"
### where "sample ID" is consistent with the "Mappable" list 
for (ID in names(Mappable)) {

	files[[ID]] <- list.files(unmapped_file_dir)[grepl(paste(ID, "_", sep=""), list.files(unmapped_file_dir))]

	for (file in files[[ID]]) {

		rt <- read.table(paste(unmapped_file_dir, file, sep=""), sep="\t", header=TRUE, row.names=1)

		### Determine chromosome based on naming convention ID_CHROMOSOME
		chr <- strsplit(file, "_|\\.")[[1]]
		chr <- chr[chr!=""][2]

		regions <- IRanges(Mappable[[ID]][Mappable[[ID]][,1]==chr,2], Mappable[[ID]][Mappable[[ID]][,1]==chr,3])

		bins <- as.numeric(rownames(rt))
		keep <- unique(queryHits(findOverlaps(IRanges(bins, bins), regions)))
		
		### Write the files with the unmappable regions removed.
		write.table(rt[keep,], paste(mapped_file_dir, ID, "_", chr, ".txt", sep=""), sep="\t", quote=FALSE, col.names=NA)

	}

}

### END OF SCRIPT TO REMOVE UNMAPPABLE REGIONS FROM HIHMM INPUT FILES
#####################################################################################################
