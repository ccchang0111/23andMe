# http://www.vincebuffalo.com/blog/2012/03/12/using-bioconductor-to-analyze-your-23andme-data.html
# https://plus.google.com/+AndrewCowanMD/posts/TThUZKUMYxV
# installing Bioconductor - gwascat
source("https://bioconductor.org/biocLite.R")
biocLite("gwascat")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene") # my 23andMe uses hg19 or "reference human assembly build 37"
biocLite("org.Hs.eg.db")
biocLite("IRanges")
# Check out 'GenomicFeatures' functions:
# browseVignettes("GenomicFeatures")
# browseVignettes("GenomicRanges")

library(gwascat)
library(ggplot2)
filepath = '/user/xxx'
d <- read.table(filepaht,
                sep="\t", header=FALSE,
                colClasses=c("character", "character", "numeric", "character"),
                col.names=c("rsid", "chrom", "position", "genotype"))

tmp <- d$chrom
d$chrom = ordered(d$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))

ggplot(d) + geom_bar(aes(chrom))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # txdb is an object now

#####################
## analyzing risks ##
#####################
data("gwrngs19")
data(ebicat37)
gwrngs.emd <- as.data.frame(elementMetadata(gwrngs19)) #elementMetadata pull out the full info, whereas gwrngs19 is just GRanges object printing limited info
dm <- merge(d, gwrngs.emd, by.x="rsid", by.y="SNPs")

risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele) # find the allele genotype (single character)
i.have.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))
}, risk.alleles, dm$genotype)
dm$i.have.risk <- i.have.risk

my.risk <- dm[dm$i.have.risk, ]
rel.cols <- c(colnames(d), "Disease.Trait", "Risk.Allele.Frequency", "p.Value", "i.have.risk", "X95..CI..text.")
head(my.risk[order(my.risk$Risk.Allele.Frequency), rel.cols], 1)
dm[which(dm$rsid == "rs10790268"), "Initial.Sample.Size"]

####################################
# search keywords in Disease.Trait #
####################################
grep("rheum", dm$Disease.Trait, ignore.case=TRUE, value=TRUE) # quick search of described disease
#rows_disease = grepl("Celiac disease", dm$Disease.Trait, ignore.case=TRUE) # search Disease.Trait column for the word "Rheuma" and isn't case sensitive
rows_disease = grepl("rheum", dm$Disease.Trait, ignore.case=TRUE) # search Disease.Trait column for the word "Rheuma" and isn't case sensitive
rows_asian = grepl("Asian", dm$Initial.Sample.Size, ignore.case = TRUE) # search just the Asian studies

my_asian_disease_data = dm[rows_disease & rows_asian,] # in the matched data "dm", look for diseases that match above two selections
# check col names colnames(dm)
# select rwos that also contain Asian sample
my_asian_disease_data[c('rsid','genotype','Strongest.SNP.Risk.Allele','num.Risk.Allele.Frequency','OR.or.beta','PUBMEDID')]