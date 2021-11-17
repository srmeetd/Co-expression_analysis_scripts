# import data
library (edgeR)
library(limma)
library (gplots)
library(goseq)
library(WGCNA)

library(biomaRt)

setwd("D:/dataset/")


Grouping = read.table("D:/coexpressed/Grouing.txt", sep = "\t", header =F)

#file_listeria = list.files(pattern = "L2_mRNA.genes.results.txt$")
file_wt = list.files(pattern = "NI2_mRNA.genes.results.txt$")
file_salmonella = list.files (pattern = "S2_mRNA.genes.results.txt$")
files = c(file_wt,file_salmonella)
file = files[grep("EU",files)]
#file = file_salmonella
names(file) = file

import = function (file)
{
  name = basename(file)
  data = read.table (file, as.is = T, header  = T, sep = "\t")
  column = data [,"expected_count"]
  names(column) = data [,"gene_id"]
  return(column)
}

data.import=lapply(file, import)
dfs = do.call(cbind, data.import)
df1 = df[,grep ("EU",colnames(df))] 
colnames(df1) = gsub("_mRNA.genes.results.txt|_S2","",colnames(df1))
colnames(df1) = gsub("^.*\\_","\\1",colnames(df1))

homo_alte_ld = read.csv ("z:/Shared/Thesis_data/SER.csv", header = FALSE)
#homo_alte = read.csv ("z:/Shared/Allelic_imbalance/hom_alt.csv", header =FALSE)
homo_alte$V1 = gsub("^.*\\-","",homo_alte$V1)
homo_alte$V1 = gsub("\\_.*","",homo_alte$V1)
homo_alte$V1 = gsub("\\s+","",homo_alte$V1)
homo_alte_ld$V1 = gsub ("trimmed-","",homo_alte_ld$V1)
ENSG00000049239 = df1[grep("ENSG00000149131",rownames(df1)),]
ENSG00000049239_hetro = ENSG00000049239[names(ENSG00000049239)%in%homo_alte_ld$V1]
#names(ENSG00000049239_hetro) = rep ("sampels",length(names(ENSG00000049239_hetro)))
ENSG00000049239_hetro = as.data.frame(ENSG00000049239_hetro)
colnames(ENSG00000049239_hetro) = "Samples"
ENSG00000049239_hetro$SNP = "Htero"

ENSG00000049239_homo = ENSG00000049239[!names(ENSG00000049239)%in%homo_alte_ld$V1]
#names(ENSG00000049239_homo) = rep ("sampels",length(names(ENSG00000049239_homo)))
ENSG00000049239_homo = as.data.frame(ENSG00000049239_homo)
colnames(ENSG00000049239_homo) = "Samples"
ENSG00000049239_homo$SNP = "Homo"
#homo_alte = intersect (homo_alte_ld$V1,homo_alte_RS7755$V1)
#homo_alte = read.csv ("z:/Shared/Alleic_imabalnce_super_Regs_miRNAs/LD_CD36/Homo_alter_rs1049673.txt", header =FALSE)
ENSG00000049239_all = rbind (ENSG00000049239_homo,ENSG00000049239_hetro)
ENSG00000049239_all$log = log2(ENSG00000049239_all$Samples)
ENSG00000049239_all$Samples =NULL
#heterozygous_ld = read.csv ("z:/Shared/Alleic_imabalnce_super_Regs_miRNAs/Hetrozygous_rs12706950.txt", header = FALSE)
htero = read.csv ("z:/Shared/Allelic_imbalance/rs7755_hetro1.csv",header=FALSE)
htero$V1  = gsub("^.*\\-|\\_.*|\\s+","",htero$V1)
#heterozygous = intersect(heterozygous_ld$V1,htero_rs7755$V1)
#heterozygous = read.csv("z:/Shared/Alleic_imabalnce_super_Regs_miRNAs/LD_CD36/hetro_rs1049673.txt",header = FALSE)

homo_ref = read.csv ("z:/Shared/Allelic_imbalance/homo_ref.csv", header = FALSE)
homo_ref$V1  = gsub("^.*\\-|\\_.*|\\s+","",homo_ref$V1)

#tmp = rbind (homo_alte,heterozygous)


#homo_ref = df1[,colnames(df1)%in%homo_ref$V1]
Grouping$V1 = gsub ("_mRNA.genes.results.txt","",Grouping$V1)
seeker=function(x)
{
  return(as.character(Grouping[Grouping[,"V1"]==x,"V2"]))
}
col_name = sapply (colnames(df1), seeker)
Group = col_name

Group = col_name
GroupColors=as.character(Group)
GroupColors=gsub("L2","Blue",GroupColors)
GroupColors=gsub("NI2","Red",GroupColors)
GroupColors=gsub("S2","Black",GroupColors)


remove_zero = df1[rowSums(df1)>0, ] 
remove_zero = remove_zero[rowMeans(remove_zero)>10,]
count = DGEList (count = remove_zero) ##
TMM_count = calcNormFactors(count, method = "TMM")

CPM_Count = cpm(TMM_count, normalized.lib.sizes=TRUE, log=T, prior.count=1)
CPM_Count_1 = CPM_Count[grep ("ENSG00000173334", rownames(CPM_Count)),]
rownames(CPM_Count_1) = NULL
names(CPM_Count_1) = gsub ("_mRNA.genes.results.txt","",names(CPM_Count_1))

bamnes = read.csv("z:/Shared/Thesis_data/Nmaes_unstimulated.csv")
bamnes$ExperimentTitle = gsub (":","_",bamnes$ExperimentTitle)

seeker=function(x)
{
  return(as.character(bamnes[bamnes[,"ExperimentTitle"]==x,"given_names"]))
}
col_name = sapply (names(CPM_Count_1), seeker)
names(CPM_Count_1) =  col_name

CPM_Count_1_htro = CPM_Count_1[names(CPM_Count_1)%in%names(IL17RA_genotype_hetrro)]
CPM_Count_1_homo_alt =  CPM_Count_1[,colnames(CPM_Count_1)%in%rownames(hetro_ref)]

CPM_Count_1_homo_ref = CPM_Count_1[names(CPM_Count_1)%in%names(IL17RA_genotype_homo)]

group <- c(rep("homozygous_reference", length(CPM_Count_1_homo_ref)), 
             rep("heterozygous", length(CPM_Count_1_htro)),
           rep("homozygous_alternate",length(CPM_Count_1_homo_alt)))
Salmonella <- data.frame(infection = c(CPM_Count_1_homo_ref,CPM_Count_1_htro,CPM_Count_1_homo_alt)
                         , group = group)

Salmonella %>%
  dplyr::mutate (group = factor(group, levels = c("homozygous_reference", 
                         "heterozygous", "homozygous_alternate"))) -> sorted_column
colnames(sorted_column) = c("Log_CPM_values","Genotypes")
sorted_column$SNP = "rs1372389110"


#tiff("z:/Shared/Thesis_data/H6_PD_rs72641816.NOCHANGE.PNG" ,
 #   width = 5.5, height = 4.5, units = 'in', res = 600)

p = ggplot(ENSG00000049239_all, aes(x=SNP, y=log,fill =SNP)) +
  geom_violin(trim=FALSE, fill = "white", color="black")+
  geom_boxplot(width=0.1,outlier.size = 0.1) + theme(panel.background = element_blank()) +
  theme(axis.text=element_text(size=8),axis.text.x = element_text(size = 8,color = "black"), 
        axis.text.y = element_text(size = 8, color = "black")) + 
  geom_smooth(method = "lm", aes(group=0.8),color = "black", lwd = 0.5)+
  theme(legend.position="right") + theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank())
  #
#p = p +  guides(fill=TRUE)
 
p 

ggsave(filename="z:/Shared/Thesis_data/H6_PD__3rd_noCHANGE.PNG", 
       plot=p,
       width = 4, height = 2, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 600)

#p=  p+ stat_compare_means(comparisons = my_comparisons,paired = FALSE)

#dev.off()

stripchart(infection~group, Salmonella, method = "overplot", offset=1/3,vertical = TRUE, pch=19,col = c("BLUE","GREEN","Red"),ylab = "Normalized")
mean = tapply(Salmonella$infection, Salmonella$group, mean)
points(c(1, 2,3), mean, pch="-", cex=3, col="black")
legend("bottom", legend=("Mean of normalized"),col=("red"), lty=1:3, cex=0.95)



pdf(file = "MDS_Plot_all_samples.pdf", title = 'MDS PLOT', paper = "a4")
plotMDS(CPM_Count_1,col=GroupColors)
dev.off()

Group=as.factor(Group)
Group=relevel(Group,ref="NI2")
#design <- model.matrix(~ 0+Group + ethnicity)
design <- model.matrix(~ 0+Group)
colnames(design) <- gsub("Group","",colnames(design))
dge=calcNormFactors(count)
v=voom(dge,design,plot=TRUE,normalize.method="none")
fit= lmFit(v, design)
contrast.matrix <- makeContrasts(L2-NI2,S2-NI2, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
comparisons=list()
comparisons[["L2-NI2"]]=1
#comparisons[["S2-NI2"]]=2
#comparisons[["L2-S2"]]=3


volcano.plot=function(x)
{
  cat("Working with ",x,".\n",sep="")
  tt=topTable(fit2, coef=comparisons[[x]], adjust="BH",number=nrow(v))
  basename=x
  write.csv(tt,file=paste(basename,".diff_exp.csv",sep=""))
  pdf(width=10,height=10,file=paste(basename,"_volcano.pdf",sep=""))
  plot(tt[,"logFC"],-log10(tt[,"adj.P.Val"]),ylab="-log10(BH QVal)",xlab="logFC",pch=19,cex=0.30,main=basename)
  abline(h=-log10(0.05),lty=3,col="red")	
  abline(v=c(-1,1),lty=3,col="red")
  dev.off()
  signif.upreg=tt[tt[,"logFC"]>=1&tt[,"adj.P.Val"]<=0.05,]
  signif.downreg=tt[tt[,"logFC"]<=-1&tt[,"adj.P.Val"]<=0.05,]
  write.csv(signif.upreg,file=paste(basename,".signif.upreg.csv",sep=""))
  write.csv(signif.downreg,file=paste(basename,".signif.downreg.csv",sep=""))
  signif.change=c(rownames(signif.downreg),rownames(signif.upreg))
  return(signif.change)
}
signif.genes=unique(unlist(lapply(names(comparisons),FUN=volcano.plot)))

signif.data=v[rownames(v)%in%signif.genes,]
listeria = signif.data[,-grep ("NI2_mRNA.genes.results.txt", colnames(signif.data))]
listeria = listeria$E
#listeria_only = v[rownames(v)%in%lister_coexpressed_trib1$V1,] 


pdf(file = "Heatmap.pdf", title = 'HEATMAP',paper = "a4r")
heatmap.2 (as.matrix(listeria),trace="none",cexCol = 0.5,ColSideColors=GroupColors)
dev.off()
getwd()

###########################
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = "85")
ensembl_wt = as.character(unlist(rownames(signif.data)))
res_wt = getBM(attributes= c('ensembl_gene_id','hgnc_id','hgnc_symbol'), 
               filters = 'ensembl_gene_id', 
               values = ensembl_wt, 
               mart = ensembl)

anno = merge(signif.data$E,res_wt,by.x = "row.names",by.y = "ensembl_gene_id")
rownames(anno) = make.names(anno$hgnc_symbol,unique = T)

trib1 = fread ("z:/Shared/Thesis_data/Macrophages_RNA_expression/modules_listeria/modulered.txt",header = FALSE)
trib1_anno = merge(trib1,res_wt,by.x = "V1",by.y = "ensembl_gene_id")


pa_background_mz =  anno
analyzed_mz =  trib1_anno



library(goseq)
library(KEGGREST)
library(plyr)

macrophages_targted_genes = unique(pa_background_mz$Row.names)
downregs_mRNA_lps_upregs_miRNAs = unique(subset(analyzed_mz, select =  V1))
gene_name = as.data.frame(macrophages_targted_genes)
colnames(gene_name) = "ID"
colnames(downregs_mRNA_lps_upregs_miRNAs) = "ensembl_transcript_id"
downregs_mRNA_lps_upregs_miRNAs$value = 1
tmp = merge (gene_name,downregs_mRNA_lps_upregs_miRNAs, by.x = "ID", by.y = "ensembl_transcript_id", all = TRUE)
tmp[is.na(tmp)] <- 0
gene = as.integer(tmp[,2])
names(gene) = tmp$ID
gene = gene[!duplicated(names(gene))]
bias_dataa = pa_background_mz[pa_background_mz$Row.names%in% names(gene),]
bias_data = rownames(bias_dataa)
bias_data = subset(bias_data, select = c("ensembl_gene_id","size"))
bias_data = bias_data[!duplicated(bias_data$ensembl_gene_id),]
rownames(bias_data) = bias_data$ensembl_gene_id
bias_data$ensembl_gene_id = NULL
colnames(bias_data) = NULL
bias_data1 = as.vector(t(bias_data))
names(bias_data1) = rownames(bias_data) 

genes = gene[names(gene)%in%names(bias_data1)]

```

```{r}
library(goseq)
pwf=nullp(genes,bias.data = bias_data1)

GO.counts=goseq(pwf,"hg19","ensGene")

getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    out[[term]] <- tmp
  }
  out
}

goterms <- GO.counts$category 
goList <- getGeneLists(pwf, goterms, "hg38", "ensGene")
GO.counts$EnsemblID <- sapply(GO.counts$category, function(x) goList[[x]])
GO.counts <- GO.counts[GO.counts$numInCat > 10,]
GO.counts$Enrich <- (GO.counts$numDEInCat/sum(gene))/(GO.counts$numInCat/length(gene))
GO.counts$FDR <- p.adjust(GO.counts$over_represented_pvalue, method="BH")
downregsmRNA_upregs_miRNAs_GO = GO.counts