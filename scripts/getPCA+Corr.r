options(stringsAsFactors=FALSE)
library(stats)
library(ade4)

argsFull <- commandArgs()
listDir=NULL
legendPos=0

helpMessage=paste("Usage: getPCA.r\n
    [Mandatory] \n
        -I, --input /path/to/pooled matrix\n\t\tRead count matrix (.txt)
        -O, --output /path/to/output file.pdf/\n\t\tFile to save the results\n
        -L, --list /path/to/list of junction file\n\t\tList of junction file (.txt)
        --legend Integer\n\t\tPosition of legend [Default=0], 0: lower left, 1: upper left, 2: upper right, 3: lower right
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript getPCA.r -I path/to/matrices -O ./outputTest/")

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
    }else if(args[i]=="--legend"){
        legendPos=as.numeric(args[i+1]);i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

dataACP = read.table(inputFile,sep="\t",header=TRUE,dec=".")
if(!is.null(listDir)){
    list = read.table(listDir,sep="\t",header=FALSE)
    dataACP = dataACP[which(dataACP$ID%in%list$V1),]
    print(dim(dataACP))
}

matrixACP = as.matrix(dataACP[,2:ncol(dataACP)])
matrixACP = t(apply(matrixACP,1,function(x){x[is.na(x)]=0;return(x)}))

res_pca=dudi.pca(as.matrix(matrixACP),scannf=F,nf=6,center=FALSE)
PC1 = res_pca$co[,1]
PC2 = res_pca$co[,2]
print(res_pca$co)

if(legendPos==0){
    xPos=min(PC1); yPos=min(PC2)+abs((max(PC2)-min(PC2))/2)
}else if(legendPos==1){
    xPos=min(PC1); yPos=max(PC2)
}else if(legendPos==2){
    xPos=min(PC1)+abs((max(PC1)-min(PC1))/2); yPos=max(PC2)
}else if(legendPos==3){
    xPos=min(PC1)+abs((max(PC1)-min(PC1))/2); yPos=min(PC2)+abs((max(PC2)-min(PC2))/2)
}

pdf(file=outputDir)
plot(PC2~PC1,pch=2:ncol(dataACP),col=rainbow(ncol(dataACP)-1),xlab="PC1",ylab="PC2")

# If the legend overlaps the drawed points you can adapt its positions by x and y value
legend(x=xPos,y = yPos,legend=colnames(dataACP)[2:ncol(dataACP)],pch=2:ncol(dataACP),col=rainbow(ncol(dataACP)-1))
dev.off()
