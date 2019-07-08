library(wesanderson)
library(RColorBrewer)
library(made4)
library(lme4)
library(nlme)
library(car)


sires = read.table("chapter2/siredata.txt", header = T)
prog = read.table("chapter2/progdata.txt", header = T)

colnames(prog)
head(prog)
dim(sires)

summary(sires$weight)

# could just use na.strings="-"
prog[which(prog$weight == "-"),]
prog = prog[-which(prog$weight == "-"),]
prog$weight = as.numeric(as.character(prog$weight))
summary(prog$weight)
summary(sires$weight)

# Filter outliers

plot(prog$weight)
prog = prog[-which(prog$weight > 400),]
summary(prog$weight)

index = grep("m", names(prog))
index


missing = numeric()
for (i in 1:length(index)) {
  missing = c(missing, which(prog[, index[i]] == "-"))
}
missing

# Remove with missing values
missing = unique(missing) 
prog = prog[-missing,]
summary(prog$m11)

for (i in index) {
  prog[, i] = factor(prog[, i])
}
summary(prog$m11)

boxplot(prog$weight ~ prog$sire, col = 1:length(levels(prog$sire)))
boxplot(prog$weight ~ prog$sex, col = 1:2, main = "Weigth x Sex")

plot(density((prog$weight)))

plot(
  prog$weight,
  col = prog$sex,
  pch = as.numeric (prog$sex),
  main = "XY plot of weight by animal",
  xlab = "animal",
  ylab = "weight"
)
legend("topleft",
       levels(prog$sex),
       col = 1:2,
       pch = 1:2)


indexms = grep("m", names(sires))
indexms = matrix(indexms, length(indexms) / 2, 2, byrow = T)
print(indexms)

indexm = grep("m", names(prog))
indexm = matrix(indexm, length(index) / 2, 2, byrow = T)
print(indexm)

compatible = matrix(NA, length(sires$id), length(indexms[, 1]))

# Check alleles of offspring vs that of sires -- most common should be shared with sire
for (j in 1:length(indexms[, 1])) {
  for (i in 1:length(sires$id)) {
    indexs = which(prog$sire == sires$id[i])
    sirealleles = sires[i, indexms[j,]]
    sirealleles = c(as.character(sirealleles[, 1]), as.character(sirealleles[, 2]))
    
    hold = prog[indexs, indexm[j,]]
    hold = factor(c(as.character(hold[, 1]), as.character(hold[, 2])))
    hold = sort(summary(hold), decr = T)
    topalleles = names(hold)[1:2]
    compatible[i, j] = length(comparelists(sirealleles, topalleles)$Set.Dif)
    if (i == 1 & j == 1) {
      cat("allele counts in offspring\n")
      print(hold)
      cat("most common alleles in offspring\n")
      print(topalleles)
      cat("sire alleles\n")
      print(sirealleles)
    }
  }
}
print(compatible)
compatible = matrix(NA, length(prog$id), length(indexms[, 1]))
# Check for mendelian inconsistancies
for (j in 1:length(indexms[, 1]))
{
  for (i in 1:length(sires$id))
  {
    indexs = which(prog$sire == sires$id[i])
    
    sirealleles = sires[i, indexms[j,]]
    sirealleles = c(as.character(sirealleles[, 1]), as.character(sirealleles[, 2]))
    
    for (k in 1:length(indexs)) {
      hold = prog[indexs[k], indexm[j,]]
      topalleles = c(as.character(hold[, 1]),
                     as.character(hold[, 2]))
      compatible[indexs[k], j] = length(comparelists(sirealleles, topalleles)$intersect)
    }
  }
}
compatible = data.frame(compatible)
for (i in 1:length(compatible[1,]))
  compatible[, i] = factor(as.character(compatible[, i]))
cat("\nSummary of alleles in common between + sires and offspring\n")
colnames(compatible) = c("M1", "M2", "M3", "M4", "M5")


m = as.matrix(compatible)
mat = mapply(m, FUN = as.numeric)
mat = matrix(data = mat,
             nrow = dim(m)[1],
             ncol = dim(m)[2])
# Heatmap visualizing the fact that each prog must contain at least 1 allele from its associated sire
# Note  single problematic entry with 0
heatmap.2(mat, col = wes_palette("Zissou1", 4))

index = which(compatible$M1 == 0)
print(prog[index,])

# Remove read that is erroneous
prog = prog[-index,]

alleles = summary(factor(c(
  as.character(prog$m11), as.character(prog$m12)
)))
alleles = alleles / sum(alleles)
alleles

hold = data.frame(m11 = as.character(prog$m11), m12 = as.character(prog$m12))
hold[, 1] = as.character(hold[, 1])
hold[, 2] = as.character(hold[, 2])
sorted = character()

for (i in 1:length(hold[, 1])) {
  sorted = rbind(sorted, sort(as.character(hold[i,]))) #Sort character vector to correct for ordering
  genotypes = paste(sorted[, 1], sorted[, 2], sep = "_")
}
genotypes = summary(factor(genotypes))
genotypes = genotypes / sum(genotypes)
# See visualization of allele frequencies

# Alternative to split.screen
par(mfrow = c(2, 1))
barplot(sort(genotypes, decr = T), col = brewer.pal(11, "Set1"))
barplot(sort(alleles, decr = T), col = wes_palette("IsleofDogs1", length(alleles)))
dev.print(
  "chapter2/GtypeAndAllelicFreq.pdf",
  device = pdf,
  width = 8.5,
  height = 11
)
dev.off()

# Note that sires have allelic frequencies shown below so we expected the progeny, given they are randomly selected from a similar population
barplot(summary(factor(c(
  as.character(sires$m11), as.character(sires$m12)
))), col = wes_palette("BottleRocket2", 4))

allgeno = NULL
for (i in 1:length(indexm[, 1]))
{
  hold = data.frame(prog[indexm[i,]])
  hold[, 1] = as.character(hold[, 1])
  hold[, 2] = as.character(hold[, 2])
  hold[, 2] = as.character(hold[, 2])
  sorted = character()
  # Note that this is nested with the same index as above
  for (i in 1:length(hold[, 1]))
    sorted = rbind(sorted, sort(as.character(hold[i,])))
  genotypes = paste(as.character(sorted[, 1]),
                    as.character(sorted[, 2]), sep = "_")
  allgeno = cbind(allgeno, genotypes)
}
colnames(allgeno) = c("M1", "M2", "M3", "M4", "M5")
markers = data.frame(prog[, 1:4], allgeno)
head(markers)

write.table(
  markers,
  "chapter2/cleandata.txt",
  quote = F,
  sep = "\t",
  row.names = F
)

attach(sires)
attach(prog)

results = lm(weight ~ M1, data = markers)
summary(results)

par(mfrow = c(2, 2))
plot(results)
dev.off()

qqnorm(predict(results), col = "blue")

fligner.test(weight ~ M1, markers)
shapiro.test(prog$weight)

anova(results)

results = lm(weight ~ . - id - sire, data = markers)

model1 = lm(weight ~ sex + sire + M1, data = markers)
model2 = lm(weight ~ sex + M1, data = markers)
model3 = lm(weight ~ M5, data = markers)
model4 = lm(weight ~ 1, data = markers)

# Note that with normal anova ordering is important -- we dont want the effect to get sucked up by the first terms
# before hitting our marker
anova(model1, model2)
anova(model2, model3)
anova(model3, model4)
anova(model2, model4)

model5 = lm(weight ~ sex * sire * M1, data = markers)
anova(model5, model1)

summary(lm(weight ~ sex + M1, data = markers[which(prog$sire == "sire1"),]))

plot.design(weight ~ sex + sire + M5, data = markers, col = "blue")
anova(lm(weight ~ sex + sire + M1, data = markers))
anova(lm(weight ~ sex + sire + M5, data = markers))



# Check for type III(marginal) errors
Anova(lm(weight ~ sex + sire + M1 + M2 + M3 + M4 + M5, data = markers), type = 3)


# Consider using linear mixed models and generalized linear models


linear = coefficients(lm(weight ~ sex + sire + M5, data = markers))
random = fixef(lme(weight ~ sex + M5, random =  ~ 1 |
                     sire, data = markers))

linear = data.frame(effect = names(linear), fixed = linear)
random = data.frame(effect = names(random), random = random)

# See http://www.bodowinter.com/tutorial/bw_LME_tutorial1.pdf
# and
# http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf

# Note merge function for combining data
comparison = merge(linear, random, by = "effect")
comparison = data.frame(comparison, difference = comparison$fixed - comparison$random)
comparison
