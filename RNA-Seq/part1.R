#starting with matrix count

library(limma)
library(edgeR)
library(dplyr)

#file path
files <- c("GSM1545535_10_6_5_11.txt", 
           "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", 
           "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", 
           "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", 
           "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
path=paste0("./data/",files)
#readDGE - edgeR
#read and merge set of files containing count data
x <- readDGE(path, columns=c(1,3))
#View(x$counts)
#View(x$samples)

#if already in a single file, import in R, convert to DGEList-object using DGEList function

#keep part of col names for simplicity. substring(string, start,end)
colnames(x) <- substring(colnames(x), 19, nchar(colnames(x)))

#sample info (group)
x$samples$group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP")) 

#batch info
x$samples$lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))

library(Mus.musculus)
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
                keytype="ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes

#convert raw counts to CPM and log-CPM

cpm <- cpm(x)
#perform on x$counts df 
#lib.size in samples df calculated by taking sum of all counts in a sample
#cpm=count per genes*e6/lib.size -> normalize count for each gene according to that sample's total count

lcpm <- cpm(x,log=TRUE)
#lcpm = log2(cpm + 2/ave.lib.size.in.million)

#filter lowly expressed genes
keep.exprs <- filterByExpr(x) #output table with boolean val for eahc gene, threshold: CPM = 10/ave.lib.size
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

#normalization factors of each sample
x <- calcNormFactors(x, method = "TMM")

#limma package
#MDS plot
plotMDS(lcpm)

#creating a design matrix and contrast (from x$samples df info)
design <- model.matrix(~0+x$samples$group+x$samples$lane)
colnames(design) <- gsub("x\\$samples\\$group", "", colnames(design))
colnames(design) <- gsub("x\\$samples\\$lane", "", colnames(design))
#design

#set up contrast w limma
contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design))
contr.matrix

#voom method
v <- voom(x, design, plot=TRUE)
View(v$targets)
#v$genes = x$genes, v$targets=x$samples

#fit a model to expression value of each gene
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 

#take info across all genes -> more precise estimates of gene-wise variablity
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit))

#treat takes into account min log-FC
tfit <- treat(vfit, lfc=1) 
dt <- decideTests(tfit) 
#summary(dt)
#View(dt@.Data)
#0=notDE,1=upregulated,-1=down-regulated
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

#examine top DE genes
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf) 
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf) 
head(basal.vs.lp)
