source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(groupName="all")

source("chapter1/Sc1.r")
#setwd("/Desktop/R_Scratch/PAGDUR")
save.image("chapter1/myworkspace")
load("chapter1/myworkspace")
ls()
rm(aa,bb)
ls()
rm(list=ls())

readLines("chapter1/Sc1.r", n=3)

mydata = read.table("chapter1/snps.txt", T, sep="\t")
mydata[1:5,1:3]
indices=c(1,7,10,2,4)
mydata$name[indices]

pheno = read.table("chapter1/animals.txt", header=T, sep="\t")
pheno

pheno2 = read.table("chapter1//animals2.txt", header=T, na.strings = "*")
pheno2

alldata=cbind(mydata,pheno)
alldata=alldata[,-1]
alldata=alldata[,c(3,4,1,2)]
alldata

alldata=alldata[order(alldata$weight,decreasing=T),]

summary(alldata)

barplot(c(summary(alldata$allele1), summary(alldata$allele2)), main="Allele counts",col=c(1,2,3,1,2,3))

pooled = c(as.character(alldata$allele1),as.character(alldata$allele2))
pooled=(summary(factor(pooled)))
barplot(pooled, col = c(1,2,3))
legend("topleft", c("-","A","B"))

head(c(as.character(alldata$allele1),as.character(alldata$allele2)))

boxplot(alldata$weight)

plot(density(alldata$weight,na.rm=T))

plot(sort(alldata$weight),col="blue",main="Sorted weights",xlab="animal",ylab="weight")
lines(sort(alldata$weight),col="red")

boxplot(alldata$weight~alldata$allele1)

write.table(alldata,"chapter1/alldata.txt", quote=F,row.names=F,sep="\t")
     
plot(density(alldata$weight,na.rm=T), main ="Density of weights", col ="blue")
dev.print("chapter1/density.pdf", device = pdf)

plot(sort(alldata$weight))
lines(sort(alldata$weight), col="red")
