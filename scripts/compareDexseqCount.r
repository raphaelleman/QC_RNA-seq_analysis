argsFull <- commandArgs()
helpMessage=paste("Usage: compareDexseqCount.r\n
    [Mandatory] \n
        -I, --input /path/to/merged count\n\t\tpath to merged count of DEXSeq
        -O, --output /path/to/output pdf/\n\t\tpdf file to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript compareDexseqCount.r -I /path/to/merged count -O ./outputTest.pdf")

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

data = read.table(inputFile,sep="\t",header=TRUE)
laboratories = levels(as.factor(data$run))

pdf(outputDir, width=20, height=7)
par(mar=c(10,4,4,1))
for(l in laboratories){
    tmp=data[data$run==l,]
    for (g in c("BRCA1","BRCA2")){
        tmp2 = tmp[,grep(g,colnames(tmp))]
        v= unlist(tmp2)
        v[v==0]=1
        t=sub(paste(g,".",sep=""),"",rep(colnames(tmp2),each=nrow(tmp2)))
        n=factor(x=t,levels=unique(t)[order(as.numeric(substr(colnames(tmp2),nchar(colnames(tmp2))-1,nchar(colnames(tmp2)))))])
        boxplot(v~n,y.log="TRUE",log="y",main=paste(l,g,sep="_"),las=2,col=rep(c("red","blue"),ncol(tmp2)/2))
    }
}
dev.off()
pdf(sub(".pdf","2.pdf",outputDir), width=15, height=7)
n=data[,"run"]
for(e in colnames(data)[grep("BRCA",colnames(data))]){
    v= data[,e]
    v[v==0]=1
    boxplot(v~n,y.log="TRUE",log="y",main=e)
}
dev.off()
