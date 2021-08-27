############################
# DEXSeq analysis
############################

# Code extract from https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/

# if library DESeq2 isn't yet installed use:
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEXSeq")
# OR
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("DEXSeq")

library(DEXSeq)

# initialize variables
inDir = NULL
sampleNames = NULL
sampleCondition = NULL
flattenedFile = NULL
libType = NULL

argsFull <- commandArgs()
helpMessage=paste("Usage: scriptDEXSeq.r\n
    [Mandatory] \n
        -I, --input /path/to/work inDir\n\t\tpath with the count files + output directory
        -F, --flattenedFile /path/to/exonCoordinate.gtf(.gff)\n\t\tpath to GTF or GFF file used for the exon/intron counting
    [Option] \n
        -N, --sampleNames name1,name2\n\t\tcommas-separated sample names, by default is the file names\n
        -C, --sampleCondition cond1,cond2\n\t\tcommas-separated sample conditions, by default is the prefix of file names (4 first characters)\n
        -P, --libType single-end,paired-end\n\t\tcommas-separated sample library type: single-end or paired-end, by default all samples are treated as paired-end
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript scriptDEXSeq.r -I /path/to/work inDir -F /path/to/refgene_exon_BRCA.gtf")

# get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<3){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        inDir=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-F"|args[i]=="--flattenedFile"){
        flattenedFile=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-N"|args[i]=="--sampleNames"){
        sampleNames=unlist(strsplit(args[i+1],",",fixed=T));i = i+2
    }else if(args[i]=="-C"|args[i]=="--sampleCondition"){
        sampleCondition=unlist(strsplit(args[i+1],",",fixed=T));i = i+2
    }else if(args[i]=="-P"|args[i]=="--libType"){
        libType=unlist(strsplit(args[i+1],",",fixed=T));i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

DEXSeq <- function(inDir=NULL, flattenedFile = NULL, sampleNames=NULL, sampleCondition=NULL, libType=NULL, outputPrefix="Dexseq_R_result"){
    if(is.null(inDir)|is.null(flattenedFile)){message(helpMessage);stop("You don't define the work directory and/or the reference file\n")}

    countFiles = list.files(inDir,".count$",full.names=F)

    if(is.null(sampleNames)){sampleNames=basename(countFiles)}

    if(is.null(sampleCondition)){sampleCondition=substr(sampleNames,1,4)}

    if(is.null(libType)){libType=rep("paired-end",length(countFiles))}

    if (length(countFiles)!=length(sampleNames) | length(countFiles)!=length(sampleCondition) | length(countFiles)!=length(libType)){
        message(helpMessage);stop("No sample names or sample condition or library type for all count files\n")
    }

    setwd(inDir)

    message(paste("Work inDir:\n\t",inDir,"\n"))

    message(paste("flattenedFile:\n\t",flattenedFile,"\n"))

    message(paste("countFiles:\n\t",paste(countFiles,collapse="\n\t"),"\n"))

    message(paste("sampleNames:\n\t",paste(sampleNames,collapse="\n\t"),"\n"))

    message(paste("sampleCondition:\n\t",paste(sampleCondition,collapse="\n\t"),"\n"))

    message(paste("libType:\n\t",paste(libType,collapse="\n\t"),"\n"))

    # reformating files
    for (i in countFiles){

    	data=readLines(i)
    	dataCount = data[grep(":",data)]
    	dataGene=data[-grep(":",data)]
    	dataError=dataGene[grep("_",dataGene)]
    	dataGene=dataGene[-grep("_",dataGene)]
    	dataGene=gsub("\"", "", dataGene, i,fixed = T)
    	dataCorr = c(paste(dataGene, dataCount, sep=""),dataError)
    	dataCorr = unlist(strsplit(split="\t",x=dataCorr, fixed=T))
    	regionId = dataCorr[seq(from=1,to=length(dataCorr),by = 2)]
    	dataCorr = c(dataCount,dataError)
    	dataCorr = unlist(strsplit(split="\t",x=dataCorr, fixed=T))
    	count = dataCorr[seq(from=2,to=length(dataCorr),by = 2)]
    	dataMerge = data.frame(regionId,count )
    	dataMergeEx = dataMerge[grep("Exo",dataMerge$regionId),]
        dataMergeInt = dataMerge[grep("Int",dataMerge$regionId),]
    	write.table(dataMergeEx,gsub('.count','.countExon',i),row.names=F,col.names=F,quote=F,sep="\t")
    }

    countFilesExon = paste(countFiles,"Exon",sep="")

    sampleTable = data.frame(row.names = basename(countFilesExon),
			condition = sampleCondition,
			libType = libType)

    lf <- lapply(countFilesExon, function(x) read.table(x, header = FALSE, stringsAsFactors = FALSE))
    if (!all(sapply(lf[-1], function(x) all(x$V1 == lf[1]$V1))))
    	stop("Count files have differing gene ID column.")

    dcounts <- sapply(lf, `[[`, "V2")
    rownames(dcounts) <- lf[[1]][, 1]
    dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ]
    rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
    colnames(dcounts) <- basename(countFilesExon)
    splitted <- strsplit(rownames(dcounts), ":")
    exons <- sapply(splitted, "[[", 2)
    genesrle <- sapply(splitted, "[[", 1)

    #import GTF file

    aggregates <- read.delim(flattenedFile , stringsAsFactors = FALSE, header = FALSE)

    colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")

    aggregates$strand <- gsub("\\\\.", "*", aggregates$strand)
    aggregates <- aggregates[which(aggregates$class == "exonic_part"),]
    aggregates <- aggregates[-grep("Int",aggregates$attr),]
    aggregates$attr <- gsub("\"", "", aggregates$attr)
    SplitAttrib = unlist(strsplit(aggregates$attr,";"))
    aggregates$gene_id = gsub(" gene_id ","",SplitAttrib[grep(" gene_id ",SplitAttrib)])
    transcripts = gsub("transcripts ","",SplitAttrib[grep("transcripts ",SplitAttrib)])
    transcripts <- strsplit(transcripts, "+",fixed=T)
    exonids <- gsub(" exonic_part_number ","",SplitAttrib[grep(" exonic_part_number ",SplitAttrib)])
    exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
    names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
    names(transcripts) <- rownames(exoninfo)
    matching <- match(rownames(dcounts), names(exoninfo))

    stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
    stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))

    treatments = unique(sampleCondition)

    # dxd object creation
    dxd <- DEXSeqDataSet(dcounts, sampleData=sampleTable, design= ~ sample + exon + condition:exon,
    				exons , genesrle , exoninfo[matching], transcripts[matching])

    dxd = estimateSizeFactors( dxd )

    dxd = estimateDispersions( dxd  , fitType='mean')

    dxd = testForDEU( dxd)

    dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

    dxr1 = DEXSeqResults( dxd )

    formulaFullModel = ~ sample + exon + condition:exon
    formulaReducedModel = ~ sample + exon

    dxd = estimateDispersions( dxd, formula = formulaFullModel  , fitType='mean')

    dxd = testForDEU( dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel)

    dxr2 = DEXSeqResults( dxd )
    write.table(dxr2,"dexseq_result.txt",quote=FALSE)
    pdf("BRCA1.pdf",width=20, height=10)
    plotDEXSeq( dxr2, "BRCA1", legend=TRUE, cex.axis=1.2, cex=1.3,lwd=2 )
    dev.off()

    pdf("BRCA2.pdf",width=20, height=10)
    plotDEXSeq( dxr2, "BRCA2", legend=TRUE, cex.axis=1.2, cex=1.3,lwd=2 )
    dev.off()
}

DEXSeq(inDir, flattenedFile, sampleNames, sampleCondition, libType)
