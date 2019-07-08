library(cluster)
library(RSQLite)
library(cluster)
library(HSAUR)
library(gplots)
library(affyPLM)

setwd("~/Documents/GitHub/MISC/AnalysisofGenomes/")

con = dbConnect(dbDriver("SQLite"))
dbname = "chapter3/SNPsmall"

dbWriteTable(
  con,
  "snpmap",
  "chapter3/SNPmap.txt",
  header = TRUE,
  append = T,
  sep = "\t"
)
dbWriteTable(
  con,
  "SNP",
  "chapter3/SNPsample.txt",
  append = T,
  header = TRUE,
  skip = 9,
  sep = "\t"
)

con = dbConnect(dbDriver("SQLite"), dbname = "chapter3/SNPsmall")
dbListTables(con)

dbListFields(con, "snpmap")
dbListFields(con, "SNP")

dbGetQuery(con, "select count (*) from snpmap")

animalids = dbGetQuery(con, "select distinct animal from SNP")
animalids = as.vector(animalids$animal)

hold = dbGetQuery(con,
                  paste("select * from SNP where animal='", animalids[1], "'", sep = ""))
dim(hold)

snpids = as.vector(dbGetQuery(con, "select distinct name from snpmap")[, 1])
snp = dbGetQuery(con, paste("select * from SNP where snp='", snpids[1], "'", sep =
                              ""))

snp$allele1 = factor(snp$allele1)
snp$allele2 = factor(snp$allele2)

summary(snp$allele1)
summary(snp$allele2)

snp = data.frame(snp, genotype = factor(
  paste(snp$allele1, snp$allele2, sep = ""),
  levels = c("AA", "AB", "BB")
))


kmData = data.frame(x = snp$x, y = snp$y)
kmm    = kmeans(kmData, 3)

#par(mfrow = c(1, 2))
plot(kmData, col = kmm$cluster)

plot(
  snp$x,
  snp$y,
  col = snp$genotype,
  pch = as.numeric(snp$genotype),
  xlab = "x",
  ylab = "y",
  main = snp$snp[1],
  cex.main = 0.9
)
legend(
  "bottomleft",
  paste(levels(snp$genotype), " (", summary(snp$genotype), ")", sep = ""),
  col = 1:length(levels(snp$genotype)),
  pch = 1:length(levels(snp$genotype))
)

length(which(snp$gcscore < 0.7))
summary(snp$gc)

alleles = factor(c(as.character(snp$allele1), as.character(snp$allele2)), levels =
                   c("A", "B"))
summary(alleles) / (sum(summary(alleles))) * 100

obs = summary(factor(snp$genotype, levels = c("AA", "AB", "BB")))
obs

hwal = summary(factor(c(
  as.character(snp$allele1), as.character(snp$allele2)
), levels = c("A", "B")))
hwal = hwal / sum(hwal)
print(hwal)

exp = c(hwal[1] ^ 2, 2 * hwal[1] * hwal[2], hwal[2] ^ 2) * sum(obs)
names(exp) = c("AA", "AB", "BB")

# Note yates correction
xtot = sum((abs(obs - exp) - c(0.5, 1, 0.5)) ^ 2 / exp)
pval = 1 - pchisq(xtot, 1)
pval

sumslides = matrix(NA, length(animalids), 4)
rownames(sumslides) = animalids
colnames(sumslides) = c("-/-", "A/A", "A/B", "B/B")

#Numeric matrix
geno = matrix(9, length(snpids), length(animalids))
for (i in 1:length(animalids)) {
  hold = dbGetQuery(con,
                    paste("select * from SNP where animal='", animalids[i], "'", sep = ""))
  hold = data.frame(hold, genotype = factor(
    paste(hold$allele1, hold$allele2, sep = ""),
    levels = c("--", "AA", "AB", "BB")
  ))
  hold = hold[order(hold$snp), ]
  
  sumslides[i, ] = summary(hold$genotype)
  temp = hold$genotype
  levels(temp) = c(9, 0, 1, 2)
  geno[, i] = as.numeric(as.character(temp))
  # change to 9 genotypes under GC score cutoff
  geno[which(hold$gcscore < 0.6), i] = 9
}

rownames(geno) = hold$snp
colnames(geno) = animalids
dim(sumslides)

#uniqSNPs = dbGetQuery(con, "select distinct snp from SNP") uniqSNPs =
#uniqSNPs[,1]
#
#temp=matrix(F,length(snpids),1) for (i in 1:1000) { hold = dbGetQuery(con,
#paste("select gcscore from SNP where snp='", uniqSNPs[i], "'", sep = ""))
#if(median(hold$gcscore)<0.5){ temp[i]=T } }



head(sumslides)

hetero = sumslides[, 3] / rowSums(sumslides[, -1])


upper = mean(hetero) + 3 * sd(hetero)
lower = mean(hetero) - 3 * sd(hetero)

out = length(which(hetero > upper))
out = out + length(which(hetero < lower))
plot(
  sort(hetero),
  1:83,
  col = "blue",
  cex.main = 0.9,
  cex.axis = 0.8,
  cex.lab = 0.8,
  ylab = "sample",
  xlab = "heterozygosity",
  main = paste(
    "Sample heterozygosity\nmean:",
    round(mean(hetero), 3),
    " sd:",
    round(sd(hetero), 3)
  ),
  sub = paste("mean: black line ", 3,
              "SD: red line number of outliers:", out),
  cex.sub = 0.8
)

#Mean and cutoffs
abline(v = mean(hetero))
abline(v = mean(hetero) - 3 * sd(hetero), col = "red")
abline(v = mean(hetero) + 3 * sd(hetero), col = "red")
dev.print(
  file = "images/samplehetero.pdf",
  device = pdf,
  width = 6,
  height = 6
)

animcor = cor(geno)
hmcol = greenred(256)
heatmap(
  animcor,
  col = hmcol,
  symm = T,
  labRow = " ",
  labCol = " "
)


# Single SNP Analysis -----------------------------------------------------

genotypes = read.table(
  "chapter3/SNPxSample.txt",
  header = T,
  sep = "\t",
  na.strings = "9",
  colClasses = "factor"
)

dim(genotypes)
for (i in 1:length(genotypes[1,])) {
  levels(genotypes[, i]) = c("AA", "AB", "BB", NA)
}

indexsnp = apply(genotypes, 1, function(x)
  length(which(is.na(x) == T)))
indexsnp = which(indexsnp == length(genotypes[1,]))
indexsample = apply(genotypes, 2, function(x)
  length(which(is.na(x) == T)))
indexsample = which(indexsample == length(genotypes[, 1]))
length(indexsample)
length(indexsnp)

genotypes = genotypes[-indexsnp, ]

# Note: Avoid using loops if at all possible. Apply is prefered.

weight = rnorm(83, 50, 10)
summary(weight)
sd(weight)
plot(density(weight), col = "blue", main = "Density Plot of Weights")
abline(v = mean(weight), col = "red")
lines(density(rnorm(83000, 50, 10)), col = "green", lty = 2)

singlesnp = function(trait, snp)
{
  if (length(levels(snp)) > 1)
    lm(trait ~ snp)
  else
    NA
}

pvalfunc = function(model) {
  if (class(model) == "lm")
    anova(model)[[5]][1]
  else
    NA
}

results = apply(genotypes, 1, function(x)
  singlesnp(weight, factor(t(x))))

pvals = lapply(results, function(x)
  pvalfunc(x))
names(results) = row.names(genotypes)
pvals = data.frame(snp = row.names(genotypes),
                   pvalue = unlist(pvals))
length(which(pvals$pvalue < 0.01))
index = sort(pvals$pvalue, index.return = T)[[2]][1:5]


estimates = NULL
for (i in 1:5)
  estimates = rbind(estimates, coefficients(summary(results[[index[i]]])))
estimates = cbind(rep(c("AA mean", "AB dev", "BB dev"), 5),
                  estimates, rep(names(results)[index], each = 3))
estimates = data.frame(estimates)
names(estimates) = c("genotype", "effect", "stderror", "t-value", "p-value", "snp")
for (i in 2:5)
  estimates[, i] = signif(as.numeric(as.character(estimates[, i])), 2)

#Note variance in answers based on random generation of weight data above
print(estimates)

map = dbGetQuery(con, "select * from snpmap")
merged = merge(pvals, map, by.x = 1, by.y = 1)

head(merged)


plot(
  merged$position[which(merged$chromosome == 1)],-log(merged$pvalue[which(merged$chromosome ==
                                                                            1)]),
  xlab = "Map Position",
  ylab = "Log Odds",
  col = "blue",
  pch = 20,
  main = "Chromosome 1"
)
abline(h = -log(0.01), col = "red")

dev.print(
  file = "images/chr1.pdf",
  device = pdf,
  width = 6,
  height = 4
)
dev.off()


# Multiple Testing --------------------------------------------------------

# Corrects for error in multiple testing
# Note that it is overly conservative and can cause extensive false negatives

bfcor = 0.01/dim(pvals)[1]

length(which(pvals$pvalue < bfcor))

sort(pvals$pvalue)[1:5]
