#Visualize the results

#In local computer, using `R`:
  
  #### Barplot for SFS
#Barplot for SFS
sfs <- scan('cache/out.sfs')
sfs <- sfs[-c(1,length(sfs))]
sfs <- sfs/sum(sfs)

barplot(sfs,names=1:length(sfs), main='SFS')
barplot(sfs, col="#cdc0b0", xlab="No. of segregating sites", 
        ylab="Minor allele frequency", 
        names=1:length(sfs))

library("data.table")
#Histgram distribution of the theta values
#theta <- read.table("cache/theta.txt", header=TRUE)
theta <- fread("cache/theta.txt", data.table =FALSE)
hist(theta$Pairwise) 
hist(theta$Pairwise, col="green", xlab="Theta pairwise values (log10)") 

#Scatter plot of the Fst values
fst <- read.table("cache/fst_win.txt", skip=1, header=FALSE)
names(fst)[c(3,5)] <- c("midp", "fst")
plot(fst$midp, fst$fst, xlab="Physical position", ylab="Fst", col="#5f9ea0", pch=16)

# install.package("data.table")
library("data.table")
## simply read in wouldn't work
gff <- fread("largedata/Zea_mays.B73_RefGen_v4.46.chromosome.Mt.gff3", skip="#", header=FALSE, data.table=FALSE)
## grep -v means select lines that not matching any of the specified patterns
gff <- fread(cmd='grep -v "#" largedata/Zea_mays.B73_RefGen_v4.46.chromosome.Mt.gff3', header=FALSE, data.table=FALSE)

names(gff) <- c("seq", "source", "feature", "start", "end", "score", "strand", "phase", "att")
table(gff$feature)

g <- subset(gff, feature %in% "gene")
g$geneid <- gsub(".*gene:|;biotype.*", "", g$att)

### + strand
gp <- subset(g, strand %in% "+") 
# nrow(gp) 75

### get the 5k upstream of the + strand gene model
gp_up <- gp
gp_up$end <- gp_up$start - 1
gp_up$start <- gp_up$end - 5000 

### get the 5k downstream of the + strand gene model
gp_down <- gp
gp_down$start <- gp_down$end + 1
gp_down$end <- gp_down$start + 5000 

g <- subset(gff, feature %in% "gene")
g$geneid <- gsub(".*gene:|;biotype.*", "", g$att)

### - strand
gm <- subset(g, strand %in% "-") 
# nrow(gm) 82

### get the 5k upstream of the + strand gene model
gm_up <- gm
gm_up$start <- gm_up$end + 1
gm_up$end <- gm_up$start + 5000 

### get the 5k downstream of the + strand gene model
gm_down <- gm
gm_down$end <- gm_down$start - 1
gm_down$start <- gm_down$end - 5000 

gup <- rbind(gp_up, gm_up)
fwrite(gup, "cache/mt_gene_up5k.txt", sep="\t", row.names = FALSE, quote=FALSE)

gdown <- rbind(gp_down, gm_down)
fwrite(gdown, "cache/mt_gene_down5k.txt", sep="\t", row.names = FALSE, quote=FALSE)

### - strand
gm <- subset(g, strand %in% "-") 
dim(gm) # 82

fwrite(g, "cache/mt_gene.txt", sep="\t", row.names = FALSE, quote=FALSE)

## Intepret the theta results
library("data.table")
library("GenomicRanges")
library("plyr")


theta <- fread("cache/theta.txt", data.table=FALSE)
names(theta)[1] <- "seq"

up5k <- read.table("cache/mt_gene_up5k.txt", header=TRUE)

### define the subject file for theta values
grc <- with(theta, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

### define the query file for genomic feature
grf <- with(up5k, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))

### find overlaps between the two
tb <- findOverlaps(query=grf, subject=grc)
tb <- as.matrix(tb)

out1 <- as.data.frame(grf[tb[,1]])
out2 <- as.data.frame(grc[tb[,2]])
### for each genomic feature, find the sites with non-missing data
out <- cbind(out1, out2[, "start"]) 
names(out)[ncol(out)] <- "pos"


#define unique identifier and merge with the thetas
out$uid <- paste(out$seqnames, out$pos, sep="_")
theta$uid <- paste(theta$seq, theta$Pos, sep="_")

df <- merge(out, theta[, c(-1, -2)], by="uid")
# for each upstream 5k region, how many theta values

mx <- ddply(df, .(geneid), summarise,
            Pairwise = mean(Pairwise, na.rm=TRUE),
            thetaH = mean(thetaH, na.rm=TRUE),
            nsites = length(uid))

get_mean_theta <- function(gf_file="cache/mt_gene_up5k.txt"){
  # gf_file: gene feature file [chr, ="cache/mt_gene_up5k.txt"]
  
  theta <- fread("cache/theta.txt", data.table=FALSE)
  names(theta)[1] <- "seq"
  
  up5k <- read.table(gf_file, header=TRUE)
  
  ### define the subject file for theta values
  grc <- with(theta, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))
  
  ### define the query file for genomic feature
  grf <- with(up5k, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))
  
  ### find overlaps between the two
  tb <- findOverlaps(query=grf, subject=grc)
  tb <- as.matrix(tb)
  
  out1 <- as.data.frame(grf[tb[,1]])
  out2 <- as.data.frame(grc[tb[,2]])
  ### for each genomic feature, find the sites with non-missing data
  out <- cbind(out1, out2[, "start"]) 
  names(out)[ncol(out)] <- "pos"
  
  #define unique identifier and merge with the thetas
  out$uid <- paste(out$seqnames, out$pos, sep="_")
  theta$uid <- paste(theta$seq, theta$Pos, sep="_")
  
  df <- merge(out, theta[, c(-1, -2)], by="uid")
  # for each upstream 5k region, how many theta values
  
  mx <- ddply(df, .(geneid), summarise,
              Pairwise = mean(Pairwise, na.rm=TRUE),
              thetaH = mean(thetaH, na.rm=TRUE),
              nsites = length(uid))
  return(mx)
}


## Plot the results
### apply the function
up5k <- get_mean_theta(gf_file="cache/mt_gene_up5k.txt")
down5k <- get_mean_theta(gf_file="cache/mt_gene_down5k.txt")
gene <- get_mean_theta(gf_file="cache/mt_gene.txt")

library("ggplot2")

up5k$feature <- "up 5k"
down5k$feature <- "down 5k"
gene$feature <- "genic"
res <- rbind(up5k, down5k)
res_nongenic <- res
res_nongenic$feature <- "intergenic"
res_t <- rbind(res, gene, res_nongenic)
ggplot(res_t, aes(x=feature, y=Pairwise, fill=feature)) + 
  geom_violin(trim=FALSE)+
  labs(title="Theta value", x="", y = "Log10 (theta)")+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_brewer(palette="Blues") + 
  theme_classic()


## Intepret the Fst results
get_mean_Fst <- function(gf_file="cache/mt_gene_up5k.txt"){
  # gf_file: gene feature file [chr, ="cache/mt_gene_up5k.txt"]
  
  fst <- fread("cache/fst_win.txt", data.table=FALSE)
  names(fst)[2] <- "seq"
  up5k <- read.table(gf_file, header=TRUE)
  ### define the subject file for Fst values
  grc_f <- with(fst, GRanges(seqnames=seq, IRanges(start=chr, end=chr)))
  ### define the query file for genomic feature
  grf <- with(up5k, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))
  
  ### find overlaps between the two
  tb_f <- findOverlaps(query=grf, subject=grc_f)
  tb_f <- as.matrix(tb_f)
  
  out1_f <- as.data.frame(grf[tb_f[,1]])
  out2_f <- as.data.frame(grc[tb_f[,2]])
  ### for each genomic feature, find the sites with non-missing data
  out_f <- cbind(out1_f, out2_f[, "start"]) 
  names(out_f)[ncol(out_f)] <- "pos"
  #define unique identifier and merge with the Fst
  out_f$uid <- paste(out_f$seqnames, out_f$pos, sep="_")
  fst$uid <- paste(fst$seq, fst$chr, sep="_")
  df_f <- merge(out_f, fst[, c(-1, -2, -3, -4)], by="uid")
  # for each upstream 5k region, how many Fst values
  names(df_f)[ncol(df_f)] <- "fst"
  
  mx_f <- ddply(df_f, .(geneid), summarise,
                Fst = mean(fst, na.rm=TRUE),
                nsites = length(uid))
  return(mx_f)
}

### apply the function
up5k_f <- get_mean_Fst(gf_file="cache/mt_gene_up5k.txt")
down5k_f <- get_mean_Fst(gf_file="cache/mt_gene_down5k.txt")
gene_f <- get_mean_Fst(gf_file="cache/mt_gene.txt")

library("ggplot2")
up5k_f$feature <- "up 5k"
down5k_f$feature <- "down 5k"
gene_f$feature <- "genic"
res_f <- rbind(up5k_f, down5k_f)
res_nongenic_f <- res_f
res_nongenic_f$feature <- "intergenic"
res_t_f <- rbind(res_f, gene_f, res_nongenic_f)
ggplot(res_t_f, aes(x=feature, y=Fst, fill=feature)) + 
  geom_violin(trim=FALSE)+
  labs(title="Fst distribution", x="", y = "Fst")+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_brewer(palette="greens") + 
  theme_classic()