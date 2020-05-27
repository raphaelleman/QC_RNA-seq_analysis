############################
# HTSeq + DESeq2 analysis
############################

# Code extract from https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/

# if library DESeq2 isn't yet installed use:
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# OR
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)

# initialize variables
directory = NULL
sampleNames = NULL
sampleCondition = NULL

argsFull <- commandArgs()
helpMessage=paste("Usage: scriptDESeq2.r\n
    [Mandatory] \n
        -I, --input /path/to/work directory\n\t\tpath with the count files + output directory
    [Option] \n
        -N, --sampleNames name1,name2\n\t\tcommas-separated sample names, by default is the file names\n
        -C, --sampleCondition cond1,cond2\n\t\tcommas-separated sample conditions, by default is the prefix of file names (4 first characters)\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript scriptDESeq2.r -I /path/to/work directory")

# get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<1){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        directory=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-N"|args[i]=="--sampleNames"){
        sampleNames=unlist(strsplit(args[i+1],",",fixed=T));i = i+2
    }else if(args[i]=="-C"|args[i]=="--sampleCondition"){
        sampleCondition=unlist(strsplit(args[i+1],",",fixed=T));i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

DESeq2 <- function(directory=NULL, sampleNames=NULL, sampleCondition=NULL, outputPrefix = "Gmax_DESeq2"){
    if(is.null(directory)){message(helpMessage);stop("You don't define the work directory\n")}

    sampleFiles = list.files(directory,".count$",full.names=F)

    if(is.null(sampleNames)){sampleNames=basename(sampleFiles)}

    if(is.null(sampleCondition)){sampleCondition=substr(sampleNames,1,4)}

    if (length(sampleFiles)!=length(sampleNames) | length(sampleFiles)!=length(sampleCondition)){
        message(helpMessage);stop("No sample names or sample condition for all count files\n")
    }

    setwd(directory)

    message(paste("Work directory:\n\t",directory,"\n"))

    message(paste("sampleFiles:\n\t",sampleFiles,"\n",collapse="\n\t"))

    message(paste("sampleNames:\n\t",sampleNames,"\n",collapse="\n\t"))

    message(paste("sampleCondition:\n\t",sampleCondition,"\n",collapse="\n\t"))

    sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

    treatments = unique(sampleCondition)

    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = directory,
                                           design = ~ condition)
    colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                          levels = treatments)

    #guts
    dds <- DESeq(ddsHTSeq,fitType = "mean")
    res <- results(dds)

    # order results by padj value (most significant to least)
    res <- res[order(res$padj),]

    # save data results and normalized reads to csv
    resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
    names(resdata)[1] <- 'gene'
    write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

    # send normalized counts to tab delimited file for GSEA, etc.
    write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

    # produce DataFrame of results of statistical tests
    mcols(res, use.names = T)
    write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

}

DESeq2(directory,sampleNames,sampleCondition)
