options(stringsAsFactors=FALSE)

argsFull <- commandArgs()

helpMessage=paste("Usage: mergeMatrixCount.r\n
    [Mandatory] \n
        -I, --input /path/to/matrices directory\n\t\tRead count matrix (.txt)
        -O, --output /path/to/output file.txt/\n\t\tFile to save the results\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript mergeMatrixCount.r -I path/to/matrices -O ./outputTest/")

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

matrices = list.files(path = inputFile, pattern = ".txt", all.files = TRUE,
           full.names = TRUE, recursive = TRUE)

print(matrices)

#Compute geometric mean for each junction
gm_mean = function(x){prod(x[x > 0])^(1/length(x))}
contigToChr <- function(text){
    text1=unlist(strsplit(text, ".", fixed = TRUE))[1]
    chr='.'
    if(substr(text1,1,2)=="NC"){
        chr=paste('chr',as.numeric(substr(text1,4,nchar(text1))),sep="")
    }
    return(chr)
}
poolMatrix = NULL
transcriptOfInterest=c("NM_000059","NM_007294","NM_007300","NR_027676","NM_007298","NM_007297","NM_007299")

for (matrix in matrices){
    print(paste("read: ",matrix,sep=""))
    data = read.table(file=matrix,sep="\t",header=TRUE)
    data = data[which(data$NM%in%transcriptOfInterest),]
    if(nchar(data$chr[1])>5){
        data$chr=unlist(mapply(contigToChr,data$chr))
    }
    print("Calculate geometric mean")
    data$gm = apply(data[c(7:ncol(data))],1,gm_mean)
    colnames(data)[ncol(data)]<-basename(matrix)
    print("Merge data")
    data$ID=paste(data[,1],data[,2],data[,3],data[,4],sep="_")
    if(is.null(poolMatrix)){
        poolMatrix=data[,c("ID",basename(matrix))]
    }else{
        poolMatrix <- merge(poolMatrix, data[,c("ID",basename(matrix))], by.x="ID", by.y="ID", all.x=TRUE, all.y=TRUE)
    }
}

write.table(poolMatrix,file = outputDir, quote = FALSE, sep = "\t",na = "0", dec = ".", row.names = FALSE, col.names = TRUE)
