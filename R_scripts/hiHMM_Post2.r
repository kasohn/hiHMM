#####################################################################################################
### Mappability script 2 of 2- hiHMM: Bayesian non-parametric joint inference of chromatin state maps (Sohn et al. 2014)
### Djordje Djordjevic - d.djordjevic@victorchang.edu.au - 2013 
###
### This file contains a script to reintroduce unmappable regions to the hiHMM output bedfiles as the "zero" state 
###
### In order to make this script generalisable you must follow a simple naming convention:
### Assumes that the hiHMM input file naming convention is "sample ID"_"Chromosome"
### Assumes each hiHMM input file contains data for ONLY ONE chromosome
### Make sure that "sample ID" is not nested / repeated, as in do not use FLY and FLY2
### 
#####################################################################################################

#####################################################################################################
### START OF SCRIPT TO REINTRODUCE UNMAPPABLE REGIONS INTO FINAL (RECOLOURED) HIHMM OUTPUT .BED FILES
###
### Note that this script requires the same "Mappable" list to be loaded from the first mappability script in hiHMM_Pre1.r

### This script uses the IRanges library.
### If it is not installed, then install it.
if (!require("IRanges")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("IRanges")
  require("IRanges")
}

options(scipen=999)

### Set the working directory to the hiHMM output folder.
setwd("../sample_analysis/output")

### Set the directory where the output files will be written to.
outdir <- "Remapped/"
dir.create(outdir)

### Set the bin size, which is used when reintroducing unmappable regions into the final chromatin state tracks.
### The bin size is the same as the hiHMM input files and the bin size used by hiHMM.
### For this reason you should only process files of the same bin size when using this script.
bin_size <- 200

### Set the location of the mappability files
mappability_dir <- "../Mappability_files/"

Mappable <- list()

### Load mappability files in bed format (chromosome \t start position \t end position)
### 1 file per sample ID
### Mappability files contains regions that ARE mappable
### Load the mappability file corresponding to each sample ID into the Mappable list
### ie. Mappable[["Sample ID"]] <- the bed file

Mappable[["fly"]] <- read.table(paste(mappability_dir, "mappable.dm3.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")
Mappable[["worm"]] <- read.table(paste(mappability_dir, "mappable.ce10.lucy.bed", sep=""), sep="\t", header=FALSE, quote="")

### Set the location of the chromosome lengths files
chr_lengths_dir <- "../Chromosome_Lengths/"

chr_lengths <- list()

### chr lengths files are in the bed format (chromosome \t start position (probably 0) \t end position)
### sample ID must be consistent.
chr_lengths[["fly"]] <- read.table(paste(chr_lengths_dir, "dm3.txt", sep=""), sep="\t", header=FALSE, quote="")
chr_lengths[["worm"]] <- read.table(paste(chr_lengths_dir, "ce10.txt", sep=""), sep="\t", header=FALSE, quote="")

### Define a suffix that identifies the bed files to be processed. "recoloured.bed" is used here
### because the final bed files have been recoloured and renamed in a previous step (hiHMM_Post1.r)
### based on how the states are named (promoter, gene body, etc.) 
final_BED_Identifier <- "recoloured.bed"

fin_files <- list()

### Load the bed files. "sample ID" must be present in the bed file name and consistent with the Mappable list
for (ID in names(Mappable)) {

	fin_files[[ID]] <- list.files()[grepl(paste(ID, ".*", final_BED_Identifier, sep=""), list.files())]

	if (!isEmpty(fin_files[[ID]])) {

		split_mappabilities <- split(Mappable[[ID]], Mappable[[ID]][,1])
		split_mappabilities <- split_mappabilities[names(split_mappabilities)%in%chr_lengths[[ID]][,1] ]

		regions <- lapply(split_mappabilities, function(X){return(IRanges(start=X[,2],end=X[,3]))})

		unmappable_regions <- lapply(names(regions), function(X,regions){

			RX <- regions[[X]]
			chr_length <- chr_lengths[[ID]][chr_lengths[[ID]][,1]==X,3]

			starts <- end(RX)+1
			ends <- start(RX)-1

			if (!0%in%start(RX)) {
				starts <- c(0,starts)
			} else {
				ends <- ends[2:length(ends)]
			}
			if (!chr_length%in%end(RX)) {
				ends <- c(ends,chr_length)
			} else {
				starts <- starts[1:length(starts)-1]
			}

			if (!isEmpty(starts)) {
				return(IRanges(starts,ends))
			}
		}, regions=regions)
		names(unmappable_regions) <- names(regions)

		for (file in fin_files[[ID]]) {
		### Double check the format of the bed file
		### If the bed file contains a header line you must use the skip=1 parameter

		#bed <- read.table(file, sep="\t", skip=1)
		bed <- read.table(file, sep="\t")

		### Keep track of the colours from the bed file
		t <- table(bed[,4], bed[,9])
		cols <- rbind(c(0,"0,0,0"), cbind(sort(as.numeric(rownames(t))), unlist(lapply(sort(as.numeric(rownames(t))), function(s){
			colnames(t)[t[s,]!=0]
		}))
		))

		split_bed <- split(bed, bed[,1])
		split_bed_regions <- lapply((split_bed), function(X) {return(IRanges(start=X[,2], end=X[,3]))})

		sm_bed <- lapply(names(split_bed_regions), function(chr, SB, BR, UR, ID) {

			chr_length <- chr_lengths[[ID]][chr_lengths[[ID]][,1]==chr,3]
			sm <- cbind(x=seq(bin_size/2, chr_length, by=bin_size), y=0)

			for (state in unique(SB[[chr]][,4])) {
				sm[,2][queryHits(findOverlaps(IRanges(sm[,1], sm[,1]), BR[[chr]][SB[[chr]][,4]==state]))] <- state
			}

			if (!is.null(UR[[chr]])) {

				sm[,2][queryHits(findOverlaps(IRanges(sm[,1], sm[,1]), UR[[chr]]))] <- 0
				return(sm)

			} else {
				return(sm)
			}

		}, SB=split_bed, BR=split_bed_regions, UR=unmappable_regions, ID=ID)

		names(sm_bed) <- names(split_bed_regions)

		final_bed <- lapply(names(sm_bed), function(chr){
			sm <- sm_bed[[chr]]
			chr_length <- chr_lengths[[ID]][chr_lengths[[ID]][,1]==chr,3]

			diff <- sm[2:nrow(sm),2]-sm[1:(nrow(sm)-1),2]
			diff[length(diff)+1] <- 1

			end_bins <- which(diff!=0)
			end_pos <- sm[end_bins,1]+bin_size/2
			start_bins <- c(1,end_bins+1)
			start_bins <- start_bins[1:(length(start_bins)-1)]
			start_pos <- sm[start_bins,1]-bin_size/2+1

			states <- sm[end_bins,2]
			colours <- cols[match(as.character(states),cols[,1]),2]
			return(cbind(chr, as.character(start_pos), as.character(end_pos), states, 1000, "+", as.character(start_pos), as.character(end_pos), colours))
		})

		ff <- do.call(rbind, final_bed)
		fname <- strsplit(file, "\\.")[[1]]

		### Write the final bed file
		fname <- paste(paste(fname[1:(length(fname)-1)], collapse="."), "ReMapped.bed", sep=".")
		write.table(ff, paste(outdir, fname, sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

		}
	}
}

options(scipen=0)

### END OF SCRIPT TO REINTRODUCE UNMAPPABLE REGIONS INTO FINAL (RECOLOURED) HIHMM OUTPUT .BED FILES
#####################################################################################################
