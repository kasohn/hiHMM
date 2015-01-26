#####################################################################################################
### Plot emission matrix script - hiHMM: Bayesian non-parametric joint inference of chromatin state maps (Sohn et al. 2014)
### Djordje Djordjevic - d.djordjevic@victorchang.edu.au - 2013 
### This script will help you generate a heatmap of the chromatin state output from hiHMM
###
### Please note: 
### You need to assign a unique name to each state BEFORE you run this script. You can do this in excel by adding a first column (no column name).
### When naming your states, the format should be	Integer State name	. The states will be ordered by the first integer (up to 2 digits), then only the State name will be displayed. 
### For example you might name your states:
###
### 1 Promoter 1
### 2 Promoter 2
### 3 Enhancer 1
### ....
###
### When using Model 1 with multiple samples, you should also add a sample identifier to the end of your state names. 
### So for example if analysing two samples (fly (F) and worm (W)) with Model 1 and 25 states, you might name your states:
###
### 1 Promoter 1 F
### 2 Promoter 2 F
### ....
### 26 Promoter 1 W
### 27 Promoter 2 W
### ....
###
### Please note that these names are used to colour the states in the final heatmap as well as when recolouring the bedfile.
### This mapping from state name to colour can be seen and changed in the Colourise function below. 
###
#####################################################################################################

#####################################################################################################
### START OF SCRIPT TO PLOT EMISSION MATRIX

### This script uses the gplots and RColorBrewer libraries.
### If these are not installed, then install them.
if (!require("gplots")) {
  require("gplots")
}
if (!require("RColorBrewer")) {
  require("RColorBrewer")
}

### VARIABLES TO SET
### Set the working directory to the hiHMM output folder.
setwd("../sample_analysis/output")

### Set the model number, i.e. 1 or 2.
m <- 2

### Set the number of conditions/samples.
c <- 2

### Read in the PRE NAMED emission matrix from the working directory.
### Note the addition of "_named" to the default output filename.
EM <- read.csv(paste("train-hihmm-model", m, ".emission_named.csv", sep = ""), row.names=1)
TM <- read.csv(paste("train-hihmm-model", m, ".transition.csv", sep = ""))

### This parameter K defines the number of states. Alternatively, can directly set this.
K <- (dim(TM)[1] / c)

### Extract the histone modifications used in the analysis.
marks <- colnames(EM)

### fworder defines the order that the histone modification marks will be displayed in the heatmap. 
### Ensure all the histone modifications used in the analysis are listed here in the desired order.
### Here is an example order we used in our sample analysis for the fly vs. worm comparison.
fworder <- c("H3K4me3","H3K4me2","H3K4me1","H2Bub","H3K27ac","H3K23ac","H3K9ac","H3K9acS10P","H4K16ac","H4K8ac","H3K79me3","H3K79me2","H3K79me1","H3K27me1","H4K20me1","H3K36me3","H3K36me1","H3K27me3","H3K9me3","H3K9me2","H3K9me1")
marks.ord  <- c(intersect(fworder, marks), setdiff(marks, fworder))

## Check whether all marks are present in the ordering
if (length(setdiff(marks, fworder)) != 0) {
    print("Marks not present in the ordering which were added to the end of the emission matrix:")
    print(setdiff(marks, fworder))
}

### These are the colours used to colour different categories of states
reds <- rev(c(rgb(1,0,0),rgb(0.9,0,0),rgb(0.8,0,0),rgb(0.7,0,0),rgb(0.6,0,0)))
yellows <- c(rgb(1,1,0),rgb(1,0.9,0.2),rgb(1,0.8,0.3),rgb(1,0.9,0.4),rgb(0.9,0.9,0))
greens <- rev(c(rgb(0.3,1,0.3),rgb(0,1,0),rgb(0.2,0.9,0.4),rgb(0,0.9,0),rgb(0.3,0.88,0.2),rgb(0,0.8,0),rgb(0,0.7,0),rgb(0.2,0.7,0.3),rgb(0,0.6,0),rgb(0.3,0.6,0.3),rgb(0,0.5,0),rgb(0.3,0.5,0.2)))
blues <- rev(c(rgb(0,0,1),rgb(0,0,0.9),rgb(0,0,0.8),rgb(0,0,0.7),rgb(0,0,0.6)))
purples <- rev(c(rgb(1,0,1),rgb(0.8,0.1,1),rgb(0.9,0.1,0.8),rgb(0.6,0,0.6),rgb(0.5,0.1,0.5),rgb(0.2,0,0.5)))

### These are the colours used to display the ChIP signal enrichment, from white to red.
hmcols <- c(rgb(1,1,1),rgb(1,0.9,0.9),rgb(1,0.8,0.8),rgb(1,0.7,0.7),rgb(1,0.6,0.6),rgb(1,0.5,0.5),rgb(1,0.4,0.4),rgb(1,0.3,0.3),rgb(1,0.2,0.2),rgb(1,0.1,0.1),rgb(1,0,0),rgb(1,0,0),rgb(1,0,0),rgb(1,0,0))

### This function maps state names to corresponding colours, returned as the vector v
colourise <- function(X){
    v <- vector()
    for (i in 1:length(X)) {
        if (grepl("Promoter",X[i], ignore.case=TRUE)) {
            v[i] <- reds[max((i%%5+1),1)]
        } else if (grepl("Enhancer",X[i], ignore.case=TRUE)) {
            v[i] <- yellows[max((i%%5+1),1)]
        } else if (grepl("Gene",X[i], ignore.case=TRUE) | grepl("Transcription",X[i], ignore.case=TRUE)) {
            v[i] <- greens[max((i%%9+1),1)]
        } else if (grepl("Repressed",X[i], ignore.case=TRUE)) {
             v[i] <- purples[max((i%%5+1),1)]
        } else if (grepl("Heterochromatin",X[i], ignore.case=TRUE)) {
             v[i] <- blues[max((i%%5+1),1)]
        } else if (grepl("LowSignal",X[i], ignore.case=TRUE) | grepl("Low Signal",X[i], ignore.case=TRUE) | grepl("Weak Signal",X[i], ignore.case=TRUE)) {
            v[i] <- brewer.pal(9,"Greys")[2+i%%3]
        } else {
            v[i] <- brewer.pal(9,"Oranges")[2+i%%3]
            print("Unknown state; coloured orange")
        }
    }
    return(v)
}

### The following code plots the heatmap for the emission matrices
### The order of samples/conditions is the same as they were defined in the 'condition' variable in the driv_hihmm.m code

### A function to extract the rows belonging to each of the emission matrices, particularly when using Model 1
### on multiple conditions, since the output emission and transition matrix files contain rows for all conditions.
bin <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

### Re-order the columns in the emission matrix according to the order of marks specified above
EM <- EM[,marks.ord]

### Store all emission matrices in a list
EM_list <- list()
if (m == 1) {
	### For Model 1, split matrix into c entries (one for each condition)
	indexes <- bin(1:nrow(EM), c)
	EM_list <- lapply(indexes, function(X){return(EM[X,])})
} else if (m == 2) {
	### For Model 2, there is only one matrix
	EM_list[[1]] <- EM
}

### Plot the heatmap for each emission matrix
for (i in 1:length(EM_list)) {
  
  EM2 <- EM_list[[i]]

	### Here we truncate the displayed emission values for easily interpretable and consistent visualisation
	EM2[EM2<(0)] <- (0)
	EM2[EM2>3.5] <- (3.5)
  
  ### You can use this line to remove the final sample identifier from your state names when analysing multiple samples using Model 1. Otherwise leave it commented out.
  #rownames(EM2) <- unlist(lapply(lapply(strsplit(rownames(EM2)," "),function(X){return(X[1:(length(X)-1)])}),function(X){paste(X,collapse=" ")}))

	### Now we sort the matrix by the first integers in the rownames
	EM3 <- EM2[sort.int(as.numeric(substr(rownames(EM2),1,2)), index.return=T)$ix,]

	### Now we determine the vector of colours that matches our chromatin states
	rowcols <- colourise(rownames(EM3))

	### Plot and save the heatmap. May need to adjust the width and height (in inches for pdf) accordingly.
  if (m == 1) {
    pdf.file <- paste("EmissionMatrix_model", m, "_cond", i, ".pdf", sep="")
    plot.title <- paste("Chromatin State Condition", i, sep=" ")
  } else if (m == 2) {
    pdf.file <- paste("EmissionMatrix_model", m, "_combined.pdf", sep="")
    plot.title <- "Chromatin State Combined"
  }
  pdf(file=pdf.file, width=10, height=10)
  ### May need to adjust the margins accordingly.
	heatmap.2(as.matrix(EM3), col=hmcols, RowSideColors=rowcols, trace="none",margins=c(16,16),notecol="darkgreen",main=plot.title,ylab="State",sepcolor="gray",sepwidth=c(0.05,0.05),colsep=0:50,rowsep=0:50,Colv=FALSE,Rowv=FALSE,cexRow=1.8,cexCol=1.9,density.info="none",dendrogram="none")
	dev.off()

}


### The following code will plot a quick heatmap of the state transition matrix for each condition, though this is less refined.

### Store all transition matrices in a list
### Note that the tranisition matrix will contain the output from all conditions, regardless of which model is used.
TM_list <- list()
indexes <- bin(1:nrow(TM), c)
TM_list <- lapply(indexes, function(X){return(TM[X,])})

### Plot the heatmap for each transition matrix
for (i in 1:length(TM_list)) {
	TM2 <- TM_list[[i]]
	
	if (m == 1) {
		rownames(TM2) <- rownames(EM_list[[i]])
		colnames(TM2) <- rownames(EM_list[[i]])
	} else if (m == 2) {
		rownames(TM2) <- rownames(EM_list[[1]])
		colnames(TM2) <- rownames(EM_list[[1]])
	}
	TM2 <- TM2[sort.int(as.numeric(substr(rownames(TM2),1,2)), index.return=T)$ix,sort.int(as.numeric(substr(colnames(TM2),1,2)), index.return=T)$ix]
	
	### Plot and save the heatmap. May need to adjust the width and height (in inches for pdf) accordingly.
	pdf(file=paste("TransitionMatrix_model", m, "_cond", i, ".pdf", sep=""), width=10, height=10)
  ### May need to adjust the margins accordingly.
	heatmap.2(as.matrix(TM2), col=hmcols, trace="none", Rowv=FALSE, Colv=FALSE, margins=c(15,15),cellnote=round(TM,2),notecol="white",main=paste("Transition Matrix Condition", i, sep=" "),ylab="State", dendrogram="none")
	dev.off()
}

### You can now open the pdf in illustrator or other vector manipulation software and finalise the look of it
###
### END OF SCRIPT TO PLOT EMISSION MATRIX
#####################################################################################################
 


#####################################################################################################
### START OF SCRIPT TO RECOLOUR BED FILES

states <- as.numeric(substr(rownames(EM_list[[1]]),1,2))

### Generate the corresponding state colours using the colourise function and colour definitions as above
rowcols <- colourise(rownames(EM_list[[1]][sort.int(as.numeric(substr(rownames(EM_list[[1]]),1,2)), index.return=T)$ix,]))
rowcols.rgb <- col2rgb(rowcols)
rowcols.rgb <- rowcols.rgb[,states]

### Get the bed files for the corresponding model and state number.
### Ensure these are the bed files that are output by hiHMM only, and don't contain any other suffixes e.g. "recoloured" 
bed.files <- list.files()[grep(paste("^hihmm.model", m, ".K", K, ".*[^.recoloured].bed", sep=""), list.files())] 

for (bed.name in bed.files) {
	
  print(paste("Recolouring", bed.name, sep=" "))
	bed <- read.table(bed.name, skip=1, stringsAsFactors=FALSE)

	### Recolour bed files
	for(i in 1:ncol(rowcols.rgb)) {
		col <- paste(rowcols.rgb[,i],collapse=",")
		bed$V9[bed$V4==i] <- col
	}

	### Renumber the states as per your naming / ordering 
	### This can be easily modified to also rename the states based on the state names
	bed$V4 <- unlist(lapply(bed$V4,function(X){return(states[X])}))

	### Rewrite the bed file adding the "recoloured" suffix
	bed.output.name <- gsub(".bed", ".recoloured.bed", bed.name)
	write.table(bed, file=bed.output.name, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	
}

### The written file is now recoloured and ready for reintroduction of unmapped regions
###
### END OF SCRIPT TO RECOLOUR .BED FILES
#####################################################################################################
