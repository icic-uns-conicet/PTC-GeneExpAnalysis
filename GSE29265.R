# Title: Sporadic vs. Post-Chernobyl Papillary vs. Anaplastic Thyroid Cancers
# Organism:	Homo sapiens
# Experiment type:	Expression profiling by array
# Summary:	We profiled the gene expression of anaplastic thyroid cancers of Belgian patients. We compared these with the expression profile of a cohort of papillary thyroid tumors both from the Chernobyl Tissues Bank (CTB) and French patients with no history of exposure to radiations, along with their patient-matched healthy adjacent thyroid.
# Link in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29265
# Overall design:	A group of nine anaplastic thyroid carcinomas (ATCs), together with twenty papillary thyroid carcinomas (PTCs) paired with their respective adjacent tissues (N), were hybridized onto Affymetrix U133 Plus 2.0 arrays. PTCs were obtained from the Ambroise Paré Hospital (Boulogne, France; sporadic PTC, n=10) and from the Chernobyl Tissue Bank (CTB, www.chernobyltissuebank.com) (post-Chernobyl PTC, n=10). ATCs were obtained from the Jules Bordet Institute (Brussels) and the Ambroise Paré Hospital (France).

# Platforms (1)
# GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# RMA signals

## Load packages for GEO dataset management
 
if(!require(BiocManager)) library(BiocManager)
if(!require(GEOquery)) library(GEOquery)
if(!require(Biobase)) library(Biobase)

## Load dataset - Gene Set Expression object
gse2 = getGEO("GSE29265") #, AnnotGPL = T, getGPL = T) 
es2=gse2$GSE29265_series_matrix.txt.gz
e2 = exprs(es2) 
dim(e2) 
# [1] 54675   49 

feno2=pData(es2) ## Samples

geno2=fData(es2) ## Genes

## Filter out ATC samples; leave in CONTROL and PTC
tejido2 = as.factor(feno2$description)
summary(tejido2)

iATC2 = which(grepl("Anaplastic",tejido2, ignore.case = TRUE))
iPTC2 = which(grepl("Papillary",tejido2, ignore.case = TRUE))
iC2 = which(grepl("control",tejido2, ignore.case = TRUE))
namesATC2 = rownames(feno2[iATC2,])
namesPTC2 = rownames(feno2[iPTC2,])
namesControl2 = rownames(feno2[iC2,])
write.csv(namesATC2, "GSE29265_samplesOUT_ATC.csv")
write.csv(namesPTC2, "GSE29265_samplesIN_PTC.csv")
write.csv(namesControl2, "GSE29265_samplesIN_control.csv")

e2 = e2[,-iATC2] # Remove ATC samples from expression matrix
feno2 = feno2[-iATC2,] # Remove ATC samples from samples data.frame
dim(e2) # [1] 54675    40 (20 PTC, 20 control)
dim(feno2) # [1] 40 59

## Filter genes
## FACTOR // Differential expression filtering with LIMMA

fc2 = rep(0, ncol(e2)) 
fc2[iPTC2] = 1
fc2 = as.factor(fc2)
summary(fc2)

library(limma)
head(e2)
dim(e2)
length(fc2)

design2 = model.matrix(~0 + fc2)
colnames(design2) <- c("CONTROL", "PTC")
fit2 = lmFit(e2, design2) # Linear model fitting
contrast.matrix2 = makeContrasts("CONTROL-PTC", levels = design2)
PTC_fits2 <- contrasts.fit(fit2, contrast.matrix2)

# Given a microarray linear model fit, compute moderated t-statistics, 
# moderated F-statistic, and log-odds of differential expression 
# by empirical Bayes moderation of the standard errors towards a common value

ebFit2 <- eBayes(PTC_fits2)
tT2 = topTable(ebFit2, number=nrow(e2), coef=1, adjust.method = "BH", p.value=0.01) # Contraste CONTROL-PTC
dim(tT2)
## EJEMPLOS uno FC positivo, otro FC negativo
head(tT2) # todos negativos
tT2[1:20,]
boxplot(e2[rownames(tT2)[15],]~fc2) # logFC=3.228284; >0, sub expr en ptc
boxplot(e2[rownames(tT2)[1],]~fc2) # logFC= -3.853652; <0, sobre expr en ptc

probes2=rownames(tT2) # 4130 probes
conGeneSymbol2 = filter(geno2, ID %in% probes2)
conGeneSymbol2 = subset(conGeneSymbol2, select=c(ID,`Gene Symbol`))
dim(conGeneSymbol2) #4130 2
length(unique(conGeneSymbol2$`Gene Symbol`)) # 2692
SYMBS2 = unique(conGeneSymbol2$`Gene Symbol`)
GEO2file = filter(GENES, GENSYMBOL %in% SYMBS2) # 2450 2
dim(GEO2file) # 2450 2
write.csv(GEO2file,paste(DIR, "GEO2.csv",sep=""))

## Cuando FC<0, sub expresados en PTC
## Cuando FC>0, sobre expresados en PTC
## Separarlos

iSOBREgeo2 = which(tT2$logFC>0) # sobre-expr en PTC: 2224 genes
probesGEO2sobre = rownames(e2[iSOBREgeo2,])

genSOBRE2 = e2[probesGEO2sobre[1],]
boxplot(genSOBRE2~fc2,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('GEO1: ',probesGEO2sobre[1] ), cex=0.5)

iSUBgeo2 = which(tT2$logFC<0) # under-expr en PTC: 1906 genes
probesGEO2sub = rownames(e2[iSUBgeo2,])

genSUB2 = e2[probesGEO2sub[1],]
boxplot(genSUB2~fc2,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('GEO1: ',probesGEO2sub[1] ), cex=0.5)

A2sobre = filter(geno2, ID %in% probesGEO2sobre)
A2sobre = subset(A2sobre, select=c(ID,`Gene Symbol`))
dim(A2sobre) # 2224 2
length(unique(A2sobre$`Gene Symbol`)) # 1742
B2sobre = unique(A2sobre$`Gene Symbol`)
B2sobre = filter(GENES, GENSYMBOL %in% B2sobre) # 
dim(B2sobre) # 1501 2
B2sobre=cbind(B2sobre, UpDown="UpPTC")

A2sub = filter(geno2, ID %in% probesGEO2sub)
A2sub = subset(A2sub, select=c(ID,`Gene Symbol`))
dim(A2sub) # 7422 2
length(unique(A2sub$`Gene Symbol`)) # 1535
B2sub = unique(A2sub$`Gene Symbol`)
B2sub = filter(GENES, GENSYMBOL %in% B2sub) # 
dim(B2sub) # 4957 2
B2sub = cbind(B2sub, UpDown="DownPTC")

GEO2completo=rbind(B2sobre,B2sub)
write.csv(GEO2completo,paste(DIR, "GEO2b.csv",sep=""))

write.csv(e2, paste(DIR, "e2.csv"))
save(e2, file = "e2.RData")

#### REVISION ####
#### genes en comun entre los assesed de los tres experimentos
set1 = rownames(CONTEOS) # matriz expr de TCGA-THCA
set2 = GENES$ENSID
length(intersect(set1, set2))
dim(GENES)
# en GENES estan todos los de TCGA, luego buscar los de GSE en GENES
View(geno2)
# en geno2 estan todos los genes de GSE
length(unique(geno2$`Gene Symbol`)) # 23521
# son mucho menos los q realmente se analizan porq el numero inicial corresponde a probes
length(intersect(geno2$`Gene Symbol`,GENES$GENSYMBOL)) # 19516
