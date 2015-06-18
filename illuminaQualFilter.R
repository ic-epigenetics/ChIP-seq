ofile <- as.character(read.table("tofilter")[[1]])
library(ShortRead)
for(i in 26:length(ofile)){
cat("Reading",ofile[i],as.character(Sys.time()),sep=" ",fill=T)
reads <- readFastq(ofile[i])
ideal.readlength <- as.numeric(names(tail(sort(table(width(reads[1:100]))),n=1)))
quals1 <- as(FastqQuality(quality(quality(reads[1:min(c(length(reads),20000000))]))),"matrix")
if(length(reads)>20000000){
        quals2 <- as(FastqQuality(quality(quality(reads[20000001:min(c(length(reads),40000000))]))),"matrix")
}
if(length(reads)>40000000){
        quals3 <- as(FastqQuality(quality(quality(reads[40000001:min(c(length(reads),60000000))]))),"matrix")    
}
if(length(reads)>60000000){
        quals4 <- as(FastqQuality(quality(quality(reads[60000001:min(c(length(reads),80000000))]))),"matrix")    
}
if(length(reads)>80000000){
        quals5 <- as(FastqQuality(quality(quality(reads[80000001:min(c(length(reads),100000000))]))),"matrix")
}
if(length(reads)>100000000){
        quals6 <- as(FastqQuality(quality(quality(reads[100000001:min(c(length(reads),120000000))]))),"matrix")    
}
trim <- apply(quals1,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
if(length(reads)>20000000){
        trim2 <- apply(quals2,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
		trim <- c(trim,trim2)
}
if(length(reads)>40000000){
        trim3 <- apply(quals3,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
        trim <- c(trim,trim3)
}
if(length(reads)>60000000){
        trim4 <- apply(quals4,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
        trim <- c(trim,trim4)
}
if(length(reads)>80000000){
        trim5 <- apply(quals5,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
        trim <- c(trim,trim5)
}
if(length(reads)>100000000){
        trim6 <- apply(quals6,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
        trim <- c(trim,trim6)
}
bad.q <- which(!is.infinite(trim) & trim<(0.8*ideal.readlength))
good.reads <- setdiff(c(1:length(reads)),bad.q)
reads <- reads[good.reads]
trim <- trim[good.reads]
trim[is.infinite(trim)] <- ideal.readlength
reads <- ShortReadQ(sread=subseq(sread(reads),start=rep(1,length(trim)),end=trim),quality=new(Class=class(quality(reads)),quality=subseq(quality(quality(reads)),start=rep(1,length(trim)),end=trim)),id=id(reads))
#no_n <- which(alphabetFrequency(sread(reads[1:100]))[,'N']==0)
newfile <- strsplit(ofile[i],split=".fastq")[[1]][1]
newfile <- paste(newfile,"_filtered.fq",sep="")
#writeFastq(reads[no_n],newfile)
writeFastq(reads,newfile,compress=F)
}
