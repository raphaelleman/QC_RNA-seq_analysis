options(stringsAsFactors=FALSE)

countType=NULL
argsFull <- commandArgs()
helpMessage=paste("Usage: mergeCountFile.r\n
    [Mandatory] \n
        -D, --directories \"/path/to/directory1;/path/to/directory2\"\n\t\tList of directories containing the count file \";\"-separated
        -n, --names \"laboratory1;laboratory1\"\n\t\tList of laboratories names \";\"-separated
        --countType htseq\n\t\tThe count type (\"htseq\" or \"dexseq\")
        -O, --output /path/to/output file.txt/\n\t\tPath to the output file\n
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript mergeCountFile.r -D \"/path/to/directory1;/path/to/directory2\" -O ./outputTest.txt --names \"laboratory1;laboratory2\" --countType htseq")

#get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<3){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-D"|args[i]=="--directories"){
        directories=args[i+1];i = i+2
    }else if(args[i]=="-n"|args[i]=="--names"){
        names=args[i+1];i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputDir=args[i+1];i = i+2
    }else if(args[i]=="--countType"){
        countType=args[i+1];i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

directories=unlist(strsplit(directories,";",fixed=TRUE))
names=unlist(strsplit(names,";",fixed=TRUE))
if (length(directories)!=length(names)){
    stop("lengths of directories and laboratories names are differents")
}

if(is.null(countType)|(countType!="htseq" & countType!="dexseq")){
    stop("You didn't define correctly the count type (\"htseq\" or \"dexseq\")")
}

readHtseqFile <- function(pathFile){
    data=read.table(pathFile,sep="\t",header=F)
    tmp = as.data.frame(t(data$V2))
    colnames(tmp)=data$V1
    return(tmp)
}
readDexseqFile <- function(pathFile){
    data=readLines(pathFile)
	dataCount = data[grep(":",data)]
	dataGene=data[-grep(":",data)]
	dataError=dataGene[grep("_",dataGene)]
	dataGene=dataGene[-grep("_",dataGene)]
	dataGene=gsub("\"", "", dataGene,fixed = T)
	if (length(dataGene) == length(dataCount)){
		dataCorr = c(paste(dataGene, dataCount, sep=""),dataError)
	}else{
		dataCorr = c(paste(dataGene, dataCount[0:length(dataGene)], sep=""),
			dataCount[(length(dataGene)+1):length(dataCount)],dataError)
	}
    dataSplit = unlist(strsplit(split="\t",x=dataCorr, fixed=T))
	regionId = dataSplit[seq(from=1,to=length(dataSplit),by = 2)]
    count = as.numeric(dataSplit[seq(from=2,to=length(dataSplit),by = 2)])
    dataMerge = data.frame(V1 = count )
    tmp = as.data.frame(t(dataMerge))
    colnames(tmp)= regionId
    return(tmp)
}
mergeCount=NULL

for (i in 1:length(directories)) {
    countFiles = list.files(path = directories[i], pattern = ".count", full.names=TRUE)
    countLab = NULL
    for(j in countFiles){
        if(countType=="htseq"){data=readHtseqFile(j)}else{data=readDexseqFile(j)}
        countLab=rbind(countLab,data)
    }
    countLab = countLab[,grep("BRCA",colnames(countLab))]
    print(names[i])
    print(dim(countLab))
    countLab$run= rep(names[i],nrow(countLab))
    mergeCount=rbind(mergeCount,countLab)
}

write.table(mergeCount,outputDir,row.names=FALSE,sep="\t",quote=FALSE)
