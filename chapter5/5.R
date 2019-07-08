
# Chapter 5 - Gene Expression Analysis ------------------------------------


# Intoduction to expression profiling studies and DGE analysis in general is extremely broad but a solid introduction if the field is unfamiliar

library(affy) #From bioconductor
library(ShortRead)
library(affyPLM)
filenames = c(paste("ctrl", 1:5, ".CEL", sep = ""),
              paste("treat", 1:5, ".CEL", sep = ""))

Names = c(paste("C", 1:5, sep = ""),
          paste("T", 1:5, sep = ""))
slides = ReadAffy(filenames = paste("chapter5/", filenames, sep = ""),
                  sampleNames = Names)



seq=readFastq("chapter5/RNAseq.fastq")
summary(seq)
head(sread(seq))
head(quality(seq))

PLM = fitPLM(slides)
par(mfrow=c(2,2))
image(slides[,1],main="Log Intensities")
image(PLM,type="weight",which=1,xlab=XLabel,main="Weights")
image(PLM,type="resids",which=1,xlab=XLabel,main="Residuals")
image(PLM,type="sign.resids",which=1,xlab=XLabel,main="Sign of residuals")
# dev.print(file="images/BadSlide.pdf", device=pdf,height=6,width=6)


