library("GSEABase")

library(DESeq2)
library(gplots)

### Load count matrix of RNAseq data and identify differentially expressed genes using DESeq2
a549<-read.table("./RNAseq_data/A549_Rawcount.txt",header=T,sep="\t")

colData<-cbind(c(1,1,1,1,2,2,2,2),c(rep("AMG",4),rep("DMSO",4)),"paired-end")
rownames(colData)<-colnames(a549)
colnames(colData)<-c("source","condition","type")

dds<-DESeqDataSetFromMatrix(countData=a549,colData=colData,design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds<-DESeq(dds)
exp_mat<-rlog(dds,fitType="local")
expmat<-assay(exp_mat)
res<-results(dds)
sort_res<-res[order(res$pvalue),]
a549_res<-subset(sort_res,padj<0.01)
a549_res<-a549_res[sort.list(a549_res[,4],decreasing=T),]

a549_up<-rownames(a549_res)[which(a549_res$log2FoldChange<0)]
a549_dn<-rownames(a549_res)[which(a549_res$log2FoldChange>0)]

DEG_list<-list(a549_up,a549_dn)

### Load signatures into a list and compare with DEGs
gene_univ<-rownames(expmat)

signature<-read.table("./Signatures/Invasive_vs_noninvasive.txt",header=T,sep="\t")
pro<-as.character(signature[which(signature$Difference>0),1])
pro<-pro[pro%in%gene_univ]
ind<-as.character(signature[which(signature$Difference<0),1])
ind<-ind[ind%in%gene_univ]

tpx2<-as.character(read.csv("./Signatures/TPX2_AURKB_node.csv",header=T)$name)
tpx2<-tpx2[tpx2%in%gene_univ]
col<-as.character(read.csv("./Signatures/COL1A2_node.csv",header=T)$name)
col<-col[col%in%gene_univ]

sig_list<-list(pro,ind,tpx2,col)

univ<-length(gene_univ)
over_tab<-matrix(NA,nr=4,nc=2)
rownames(over_tab)<-paste(c("Pro-invasivev","Indolent","TPX/AURKB""),lengths(sig_list),sep=":")
colnames(over_tab)<-paste(c("A549_AMG_up","A549_AMG_dn"),lengths(DEG_list),sep=":")
for(i in 1:4){
        gs1<-sig_list[[i]]
        len1<-length(gs1)
        for(j in 1:2){
                gs2<-DEG_list[[j]]
                len2<-length(gs2)
                overlap<-intersect(gs1,gs2)
                over<-length(overlap)
                sel_mat<-matrix(c(univ-len1-len2+over,len1-over,len2-over,over),nr=2)
                OR<-sprintf("%.2f",fisher.test(sel_mat,alternative="greater")$estimate)
                pva<-fisher.test(sel_mat,alternative="greater")$p.value
                pva<-ifelse(pva>0.01,sprintf("%.2f",pva),pva)
                print(sort(overlap))
                over_tab[i,j]<-paste(over,OR,pva,sep="/")
        }
}

### Functional enrichment of the DEGs against HALLMARK genesets
source("./0_msigdb_FET.R")

up_h<-msigdb_fish(a549_up,gene_univ,"h")
up_h<-up_h[sort.list(as.numeric(up_h[,5])),]
dn_h<-msigdb_fish(a549_dn,gene_univ,"h")
dn_h<-dn_h[sort.list(as.numeric(dn_h[,5])),]

par(mfrow=c(1,2))
up_h_sel<-up_h[which(as.numeric(up_h[,6])<0.01),]
up_h_sel<-up_h_sel[1:10,]
up_h_nam<-up_h_sel[,1]
up_h_pva<- -log10(as.numeric(up_h_sel[,6]))
dn_h_sel<-dn_h[which(as.numeric(dn_h[,6])<0.01),]
dn_h_sel<-dn_h_sel[1:10,]
dn_h_nam<-dn_h_sel[,1]
dn_h_pva<- log10(as.numeric(dn_h_sel[,6]))
bp<-barplot(rev(dn_h_pva),col="lightblue",horiz=T,xlab="-log10(FDR)")
text(-0.3,bp,rev(dn_h_nam),pos=2)
bp<-barplot(rev(up_h_pva),col="pink",horiz=T,xlab="-log10(FDR)")
text(0.3,bp,rev(up_h_nam),pos=4)


