options(stringsAsFactors=FALSE)
library(gplots)

argsFull <- commandArgs()
plot=FALSE
helpMessage=paste("Usage: getHeatmap.r\n
    [Mandatory] \n
        -I, --input /path/to/pooled matrix\n\t\tRead count matrix (.txt)
        -O, --output /path/to/output file.pdf/\n\t\tFile to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript getHeatmap.r -I path/to/matrices -O ./outputTest/")

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
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

######################################
# Convert count matrix
######################################
dataMatrix = read.table(inputFile,sep="\t",header=TRUE)
newMatrix = as.matrix(dataMatrix[,2:ncol(dataMatrix)])
row.names(newMatrix) = dataMatrix[,1]
newMatrix = t(apply(newMatrix,1,function(x){x[x!=0]=1;return(x)}))
dataToHeatmap = t(newMatrix)%*%newMatrix
print(dataToHeatmap)

######################################
# Heatmap
######################################
col<- colorRampPalette(c("blue", "red"), bias = 0.9)(nrow(dataToHeatmap)*ncol(dataToHeatmap))

pdf(outputDir, width = 12, height = 12)

heatmap.2(dataToHeatmap, Rowv=T, Colv="Rowv", scale='row',
    symm =T, col=col, cellnote = dataToHeatmap, notecex=1.5,
    notecol="white", trace="none", key = TRUE,
    margins = c(10,10), density.info="density")

dev.off()
