#starting with matrix count

library(limma)
library(edgeR)

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

