library(made4)
library(rrBLUP)
# Beyond GWAS


# Notes on R as a language ------------------------------------------------

# Note that R doesn't like many traditional programming practices and
# is fairly picky for computationally intensive manipulation
# Memory reallocation provides a simplistic example of where proper programming practices can help.
#
# View apply as "loop hiding" and proper use of compiled functions as ideal
geno = readRDS("chapter4/genotypes.rds")

# Proper use of vectorization and compiled function over loops
freqB = rowSums(geno) / (dim(geno)[2] * 2)
freqA = 1 - freqB

pheno = rnorm(2000, 100, 10)
pvmat = numeric(10000)
pvloopfunc = function() {
  for (i in 1:10000) {
    pvmat[i] = coef(summary(lm(pheno ~ geno[i,])))[2, 4]
  }
}

pvapplyfunc = function() {
  apply(geno, 1, function(x)
    coef(summary(lm(pheno ~ x)))[2, 4])
}

system.time(pvloopfunc())
system.time(pvapplyfunc())



# Note that M*M2 does elementwise multiplication and %*% for proper matrix multiplication
# We also have t(M)


# Linear Unbiased Prediction ----------------------------------------------

gwas = readRDS("chapter4/gwasData.rds")
pheno = read.table("chapter4/phenotypes.txt", header = T, sep = "\t")$Pheno
map = read.table("chapter4/map.txt", header = T, sep = "\t")

# Single regression onto individual SNPs

effect = numeric(10000)
pvalues = numeric(10000)
for (i in 1:10000) {
  res = coef(summary(lm(pheno ~ gwas[i,])))[2, c(1, 4)]
  effect[i] = res[1]
  pvalues[i] = res[2]
}
head(effect)
head(pvalues)

# Emphasis on additive effect of genotypes and the major faults of assuming only additive effects

head(effect)
head(sort(pvalues))
tail(sort(pvalues))

par(mfrow = c(2, 1))

plot(
  gwas[which(pvalues == min(pvalues)),],
  pheno,
  xlab = "Genotypes",
  ylab = "Phenotypes",
  main = paste0("Effect Size:", round(effect[which(pvalues == min(pvalues))], 2))
)

#Corresponding model
abline(lm(pheno ~ gwas[which(pvalues == min(pvalues)),]),
       col = "red", lwd = 2.5)

plot(
  gwas[which(pvalues == max(pvalues)),],
  pheno,
  xlab = "Genotypes",
  ylab = "Phenotypes",
  main = paste0("Effect Size:", round(effect[which(pvalues == max(pvalues))], 2))
)
#Corresponding model
abline(lm(pheno ~ gwas[which(pvalues == max(pvalues)),]),
       col = "red", lwd = 2.5)
dev.print(
  file = "images/MinXMaxPValueModel.pdf",
  device = pdf,
  width = 6,
  height = 6
)


y = pheno
X = matrix(0, 2000, 2)
X[, 1] = 1
interceptM = numeric(10000) # beta0 or a (intercept)effectM=numeric(10000) # beta1 or b (slope)
effectM = numeric(10000)

for (i in 1:10000) {
  X[, 2] = gwas[i, ]
  XtX = t(X) %*% X
  lhs = solve(XtX)
  rhs = t(X) %*% y
  sol = lhs %*% rhs
  interceptM[i] = sol[1, 1]
  effectM[i] = sol[2, 1]
}

yhat = interceptM[1] + effectM[1] * gwas[1,]

#snpBLUP
h2 = 0.5
y = pheno
p = rowMeans(gwas) / 2
# From eq 4.7 on pg 118
d = 2 * sum(p * (1 - p))
ve = var(y) * (1 - h2) # residual variance
va = var(y) * h2 # additive variance
lambda = d * (ve / va)

# From 4.6 on pg 118
X = t(gwas - (p * 2)) # freq. adjusted X matrix
#Computationally intensive -- dont worry about the time this takes
XtX = t(X) %*% X
diag(XtX) = diag(XtX) + lambda

ones = rep(1, length(y))
oto = t(ones) %*% ones # 1n'1n
otX = t(ones) %*% X # 1n'X
Xto = t(X) %*% ones # X'1n

# Full matrix
lhs = rbind(cbind(oto, otX), cbind(Xto, XtX))
oty = t(ones) %*% y
Xty = t(X) %*% y
rhs = rbind(oty, Xty)
effectBLUP = solve(lhs) %*% rhs
head(effectBLUP)

# Skip mean fitting
XtX = t(X) %*% X
diag(XtX) = diag(XtX) + lambda
Xty = t(X) %*% pheno
effectBLUP2 = solve(XtX) %*% Xty # SNP solutions

#Mean is the first value in the solutions matrix
# we want to skip that.
plot(effect, effectBLUP[2:10001])

trueQTLs = read.table("chapter4/trueQTL.txt", header = T, sep = "\t")

plot(
  effect,
  pch = 20,
  cex = 2.5,
  xlab = "SNP",
  ylab = "Effect",
  main = "SNP Effects"
)
abline(v = trueQTLs$indexQTL, col = "gray")
points(effectBLUP,
       col = "red",
       pch = 17,
       cex = 1.5)
points(
  trueQTLs$indexQTL,
  trueQTLs$QTLval,
  col = "green",
  pch = 17,
  cex = 1.5
)

legend(
  "topright",
  c("single SNP", "snpBLUP", "true QTL"),
  fill = c("black", "red", "green"),
  cex = 0.8
)
effectBLUP = effectBLUP[-1] # Remove mean estimate from vector
pvalBlup = 2 * pt(-abs(effectBLUP / sd(effectBLUP)), df = length(effectBLUP) - 1)
sigSNP = which(pvalBlup < 0.01 / length(pvalBlup))
length(sigSNP)

length(intersect(sigSNP, trueQTLs$indexQTL))[1]


sigSNP = which(pvalues < 0.01 / length(pvalues))
length(sigSNP)

length(intersect(sigSNP, trueQTLs$indexQTL))[1]

mbcol = as.factor(map$chrom)
chrs = unique(map$chrom)
mapcol = getcol(length(chrs))
chrseps = numeric(length(chrs))
xdist = length(map[, 3])
cum = 0
chrpos = numeric(length(chrs))
for (i in 1:length(chrs)) {
  index = which(map[, 2] == chrs[i])
  chrseps[i] = index[length(index)]
  xdist[index] = map[, 3][index] + cum
  cum = cum + map[, 3][index[length(index)]]
  chrpos[i] = cum
}

chrpos[2:length(chrs)] = chrpos[2:length(chrs)] - ((chrpos[2:length(chrs)] -
                                                      chrpos[1:(length(chrs) - 1)]) / 2)
chrpos[1] = chrpos[1] / 2

par(mfrow = c(2, 1))
plot(
  xdist,
  -log10(pvalues),
  col = (mapcol[mbcol]),
  pch = 20,
  xlab = "Chr",
  ylab = expression(paste("-", log[10], "p-value", sep = "")),
  axes = F,
  main = "single SNP regressions"
)

abline(h = -log10(0.05 / nrow(map)), lty = 2)
axis(1, at = chrpos, labels = chrs, las = 1)
axis(2, )

plot(
  xdist,-log10(pvalBlup),
  col = (mapcol[mbcol]),
  pch = 20,
  xlab = "chromosome",
  ylab = expression(paste("-", log[10], "p-value", sep = "")),
  axes = F,
  main = "snpBLUP"
)
abline(h = -log10(0.05 / nrow(map)), lty = 2)
axis(1, at = chrpos, labels = chrs, las = 1)
axis(2,)

dev.print(
  file = "images/sigSNP.pdf",
  device = pdf,
  width = 6,
  height = 6
)

sol = mixed.solve(y, X)
effPack = sol$u
cor(effPack, effectBLUP)

out = data.frame(single = effSingle,
                 snpblup = effBlup,
                 rrblup = effPack)

write.table(
  out,
  "chapter4/effects.txt",
  quote = F,
  sep = "\t",
  row.names = F
)

# Genomic Prediction ------------------------------------------------------

# Notes on partitioning phenotypic variance


# Prediction with snpBLUP -------------------------------------------------
gwas = readRDS("chapter4/gwasData.rds")
pheno = read.table("chapter4/phenotypes.txt", header = T, sep = "\t")$Pheno
effect = read.table("chapter4/effects.txt", header = T, sep = "\t")

tgv = read.table("chapter4/trueGeneticValue.txt",
                 header = T,
                 sep = "\t")$TGV

p = rowMeans(gwas) / 2
X = t(gwas - (p * 2))
pred = X %*% effect$snpblup
head(pred)

par(mfrow = c(1, 2))
plot(pheno,
     pred,
     ylab = "Predicted",
     xlab = "Observed",
     pch = 20)
model1 = lm(pred ~ pheno)
abline(model1, lwd = 2, col = "blue")

# correlation
cor(pheno, pred)

plot(tgv,
     pred,
     ylab = "Predicted",
     xlab = "True",
     pch = 20)
model2 = lm(pred ~ tgv)
abline(model2, col = "blue")
summary(model2)
# correlation with true additive value
cor(pred, tgv)

dev.print(
  file = "images/ModelComparison.pdf",
  device = pdf,
  width = 8,
  height = 6
)

# Validation

validG = readRDS("chapter4/validGeno.rds")
validP = read.table("chapter4/validPheno.txt", header = T, sep = "\t")$Pheno
validT = read.table("chapter4/validTGV.txt", header = T, sep = "\t")$TGV

validX = t(validG - (p * 2))
validPred = validX %*% effect$snpblup

plot(validT,
     validPred,
     ylab = "Predicted",
     xlab = "Observed",
     pch = 20)
model = lm(validPred ~ validT)
abline(model, lwd = 2, col = "blue")
summary(model)

dev.print(
  file = "images/Validation.pdf",
  device = pdf,
  width = 8,
  height = 6
)

predS = X %*% effect$single
plot(tgv,
     predS,
     ylab = "Predicted(via single SNP regressions)",
     xlab = "observed phenotypes",
     pch = 20)
model = lm(predS ~ tgv)
abline(model, lwd = 1.5, col = "green")

dev.print(
  file = "images/ValidationSingleSNP.pdf",
  device = pdf,
  width = 6,
  height = 6
)

cor(validT, validPred)


# gBLUP relationship matrix of similarities

freqAvg = rowMeans(gwas, na.rm = T)
p = freqAvg / 2
M = gwas - 1
P = 2 * (p - 0.5)
W = M - P
WtW = t(W) %*% W
d = 2 * sum(p * (1 - p))
G = WtW / d



# Run through gBLUP
h2 = 0.5 # Note that this was set as lambda in snpBLUP -- 50% of variation to additive genetic value
y = pheno
ve = var(y) * (1 - h2)
va = var(y) * h2
lambda = ve / va
Z = matrix(0, ncol(gwas), ncol(gwas))
diag(Z) = 1
ZtZ = t(Z) %*% Z
Gil = solve(G) * lambda
lhs = solve(ZtZ + Gil)
rhs = t(Z) %*% y
sol = lhs %*% rhs


#Compare gBLUP and snpBLUP -- Note that R rounds this to 1.
cor(pred, sol)

# gBLUP does not calculate SNP effects and instead relies on genomic relationships

#We can solve for SNP effects
backSolve = 1 / d * W %*% solve(G) %*% sol


#Finally, predict given unknown phenotypes
All = cbind(gwas, validG)
freqAvg = rowMeans(gwas, na.rm = T)
p = freqAvg / 2
M = All - 1
P = 2 * (p - 0.5)
W = M - P
WtW = t(W) %*% W
d = 2 * sum(p * (1 - p))
Gall = WtW / d

missindex = 2001:4000
Ginv = Gall[-missindex, -missindex]


diag(Ginv) = diag(Ginv) + lambda
Ginv = solve(Ginv)
solAll = Gall[, -missindex] %*% Ginv %*% pheno
cor(validPred, solAll[2001:4000])


# Pop Gen -----------------------------------------------------------------

sos = readRDS("chapter4/sosData.rds")

M = matrix(NA, nrow(sos), 3)
colnames(M) = c("hanwoo", "angus", "brahman")

M[, 1] = apply (sos[, which(colnames(sos) == "hanwoo")], 1, function(x)
  sum(x) / (length(x) * 2))

M[,2]=apply (sos[,which(colnames(sos)=="angus")],1,function(x) sum(x)/(length(x)*2))

M[,3]=apply (sos[,which(colnames(sos)=="brahman")],1,function(x) sum(x)/(length(x)*2))

rMeans = rowMeans(M)
alleleVar = rMeans * (1 - rMeans) #This is simply p*q (freqency variance)

meanDev = M - rMeans
FST = meanDev ^ 2 / alleleVar
head(FST)


plot(FST[,1],type="l",xlab="SNP", ylab="Fst",col="gray")
smoothed=runmed(Fst[,1],k=101,endrule="constant") #Run median instead of normall rollavg
lines(smoothed,type="l",col="red",lwd=2)

FSTm=mean(smoothed)
FSTsd=sd(smoothed)
sig=FSTm+3*Fstsd #Outliers
abline(h=sig)

dev.print(
  file = "images/Fst_1.pdf",
  device = pdf,
  width = 6,
  height = 6
)


#Note some compatibility issues with pegas

# See Book for Other Tools/Examples of Common Procedures ------------------


