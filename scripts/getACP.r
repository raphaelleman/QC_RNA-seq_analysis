options(stringsAsFactors=FALSE)
library(stats)
library(ade4)

argsFull <- commandArgs()
listDir=NULL

helpMessage=paste("Usage: getACP.r\n
    [Mandatory] \n
        -I, --input /path/to/pooled matrix\n\t\tRead count matrix (.txt)
        -L, --list /path/to/list of junction file\n\t\tList of junction file (.txt)
        -O, --output /path/to/output file.pdf/\n\t\tFile to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript getACP.r -I path/to/matrices -O ./outputTest/")

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
    }else if(args[i]=="-L"|args[i]=="--list"){
        listDir=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

dataACP = read.table(inputFile,sep="\t",header=TRUE)
if(!is.null(listDir)){
    list = read.table(listDir,sep="\t",header=FALSE)
    dataACP = dataACP[which(dataACP$ID%in%list$V1),]
    print(dim(dataACP))
}

matrixACP = as.matrix(dataACP[,2:ncol(dataACP)])

res_pca=dudi.pca(as.matrix(matrixACP),scannf=F,nf=6,center=FALSE)
PC1 = res_pca$co[,1]
PC2 = res_pca$co[,2]
print(res_pca$co)
pdf(file=outputDir)
plot(PC2~PC1,pch=2:ncol(dataACP),col=rainbow(ncol(dataACP)-1),xlab="PC1",ylab="PC2")

# If the legend overlaps the drawed points you can adapt its positions by x and y value
legend(x=min(PC1),y = (min(PC2)+0.4*abs(min(PC2)-max(PC2))),legend=colnames(dataACP)[2:ncol(dataACP)],pch=2:ncol(dataACP),col=rainbow(ncol(dataACP)-1))
dev.off()
