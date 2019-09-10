options(stringsAsFactors=FALSE)

argsFull <- commandArgs()
helpMessage=paste("Usage: compareHtseqCount.r\n
    [Mandatory] \n
        -I, --input /path/to/merged count\n\t\tpath to merged count of HTSeq
        -O, --output /path/to/output pdf/\n\t\tpdf file to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript compareHtseqCount.r -I /path/to/merged count -O ./outputTest.pdf")

#get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<3){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        inputFile=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputDir=args[i+1];i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

data =read.table(inputFile,sep="\t",header=TRUE)

pdf(outputDir, width=10, height=7)
boxplot(BRCA1~run,data=data[data$BRCA1!=0,],
    log="y",ylog=T,
    ylab = "Read count (log scale)",
    main="BRCA1")
boxplot(BRCA2~run,data=data[data$BRCA2!=0,],
    log="y",ylog=T,
    ylab = "Read count (log scale)",
    main="BRCA2")
dev.off()
