options(stringsAsFactors=FALSE)
library(VennDiagram)
library(gplots)

argsFull <- commandArgs()
plot=FALSE
helpMessage=paste("Usage: getVennDiagram.r\n
    [Mandatory] \n
        -I, --input /path/to/pooled matrix\n\t\tRead count matrix (.txt)
        -L, --list /path/to/list of junction file\n\t\tList of junction file (.txt)
        -O, --output /path/to/output file.pdf/\n\t\tFile to save the results\n
        --plot \n\t\tPrint plot (only if less than 5 ensembles)
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript getVennDiagram.r -I path/to/matrices -O ./outputTest/")

#get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
listDir=NULL
i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        inputFile=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputDir=args[i+1];i = i+2
    }else if(args[i]=="-L"|args[i]=="--list"){
        listDir=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="--plot"){
        plot=TRUE;i = i+1
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

######################################
# get VENN diagramm
######################################
dataVenn = read.table(inputFile,sep="\t",header=TRUE)
if(!is.null(listDir)){
    list = read.table(listDir,sep="\t",header=FALSE)
    dataVenn = dataVenn[which(dataVenn$ID%in%list$V1),]
    print(dim(dataVenn))
}
labNames=colnames(dataVenn)[2:ncol(dataVenn)]

eval(parse(text=paste(labNames,"_ID<-dataVenn$ID[dataVenn$",labNames,"!=0]",sep="",collapse=";")))
eval(parse(text=paste("vennList<-list(",paste(labNames,"ID",sep="_",collapse=","),")",sep="")))

if(plot){pdf(outputDir)}
tmp = venn(vennList,names = labNames,show.plot=plot)
summary(attr(tmp,'intersections'))
if(plot){dev.off()}
