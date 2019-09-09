options(stringsAsFactors=FALSE)
library(stats)

argsFull <- commandArgs()
helpMessage=paste("Usage: getACP.r\n
    [Mandatory] \n
        -I1, --input1 /path/to/untreated data\n\t\tNMD inhibitor untreated data of IVS10-2 (.txt)
        -I2, --input2 /path/to/treated data\n\t\tNMD inhibitor treated data of IVS10-2 (.txt)
        -O, --output /path/to/output folder/\n\t\tFolder to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript getACP.r -I path/to/matrices -O ./outputTest/")

#get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<3){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-I1"|args[i]=="--input1"){
        inputFile1=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-I2"|args[i]=="--input2"){
        inputFile2=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputDir=args[i+1];i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

#####
#import data
#####
#here  the decimal separator is an ","

dataUntreated = read.table(inputFile1,header=TRUE, sep="\t",dec=",")
dataTreated = read.table(inputFile2,header=TRUE, sep="\t",dec=",")

panel.hist = function(x,...)
{
usr =par("usr");on.exit(par(usr))
par(usr=c(usr[1:2],0,1.5))
h = hist(x,plot=FALSE)
breaks = h$breaks; nB = length(breaks)
y=h$counts; y = y/max(y)
rect(breaks[-nB],0, breaks[-1],y,col="cyan",...)
}
panel.cor = function(x,y,digits=2, prefix="",cex.cor,...)
{
usr = par("usr"); on.exit(par(usr))
par(usr = c(0,1,0,1))
r = abs(cor(x[!is.na(x)&!is.na(y)],y[!is.na(x)&!is.na(y)],method = "pearson"))
txt = format(c(r,0.123456789),digits=digits)[1]
txt=paste(prefix, txt, sep="")
if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
text(0.5,0.5,txt,cex=cex.cor * r)
}
setwd(outputDir)

dataCorrUntrPercent = dataUntreated[dataUntreated$nbOccur>=4,which(substr(colnames(dataUntreated),1,2)=="P_")]
dataCorrTrPercent = dataTreated[dataTreated$nbOccur>=4,which(substr(colnames(dataTreated),1,2)=="P_")]
pdf("CorrExprUntreatedIVS10-2.pdf")
pairs(dataCorrUntrPercent,lower.panel=panel.smooth, diag.panel=panel.hist,upper.panel=panel.cor)
dev.off()
pdf("CorrExprTreatedIVS10-2.pdf")
pairs(dataCorrTrPercent,lower.panel=panel.smooth, diag.panel=panel.hist,upper.panel=panel.cor)
dev.off()

dataCorrUntrRead = dataUntreated[dataUntreated$nbOccur>=4,which(substr(colnames(dataUntreated),1,2)=="IV")]
dataCorrTrRead = dataTreated[dataTreated$nbOccur>=4,which(substr(colnames(dataTreated),1,2)=="IV")]
logRead <- function(x){x[x>0&!is.na(x)]=log(x[x>0&!is.na(x)]);return(x)}

dataCorrUntrRead = apply(dataCorrUntrRead,2,logRead)
dataCorrTrRead = apply(dataCorrTrRead,2,logRead)

pdf("CorrReadUntreatedIVS10-2.pdf")
pairs(dataCorrUntrRead,lower.panel=panel.smooth, diag.panel=panel.hist,upper.panel=panel.cor)
dev.off()
pdf("CorrReadTreatedIVS10-2.pdf")
pairs(dataCorrTrRead,lower.panel=panel.smooth, diag.panel=panel.hist,upper.panel=panel.cor)
dev.off()
