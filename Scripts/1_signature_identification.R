hg19<-read.table("./Information/refgene_hg19.txt",header=T,sep="\t")
gene_tab<-read.table("./Information/Human_entrezid_symbol_all.txt",header=T,sep="\t")


### Load RNaseq matrix
sam_info<-read.csv("../../NIHADKEY_UPDATED.csv")
sam_gen<-as.character(sam_info$Gender)
new_mat<-read.table("./RNAseq_data/NIHAD_FPKM_normalized.txt",header=T,sep="\t")


### Confirm sex of samples using a Y-chromosome gene, RPS4Y1
rps4y1<-as.numeric(new_mat["RPS4Y1",match(sam_info[,1],colnames(new_mat))])
pdf("NIHAD_gender.pdf")
hist(as.numeric(rps4y1),br=20,xlab="FPKM - RPS4Y1 (log2)",main=NULL)
abline(v=4,col="red",lty=2)
dev.off()
rps4y1_gen<-ifelse(rps4y1>4,"M","F")
gen_tab<-cbind(as.character(sam_info[,1]),sam_gen,rps4y1_gen)
colnames(gen_tab)<-c("Sample","Annotated","Predicted")
write.table(gen_tab,"NIHAD_sample_genders.txt",row.names=F,col.names=T,quote=F,sep="\t")


### Load sample information
ref_his<-c("AC","MP","PAP","SOL","LPA","AIS","MIA","")
col_map_histo<-c("orange","red","brown","red4","yellow","lightblue","blue","gray")
histo_col<-col_map_histo[match(histo,ref_his)]


### Unsupervised clustering using most varying genes having larger variance (>1)
vari<-apply(new_mat,1,var)
cut<-vari[sort.list(vari,decreasing=T)][1000]
sel_mat<-as.matrix(new_mat[vari>=cut,])
sel_chr<-as.character(hg19[match(rownames(sel_mat),hg19$name2),3])
sel_gen<-as.character(gen_tab[match(colnames(sel_mat),gen_tab[,1]),3])
gen_col<-ifelse(sel_gen=="M","blue","red")
gene_col<-ifelse(sel_chr%in%c("chrX","chrY"),"yellow","green")

pdf("NIHAD_FPKM_variance.pdf")
hist(vari,br=50,xlab="Variance (log2(FPKM))",xlim=c(0,5),cex.lab=1.5,cex.axis=1.2)
abline(v=1,lwd=2,lty=2,col="red")
dev.off()


### Invasive and non-invasive tumors are determined by histological information of samples
pdf("NIHAD_FPKM_1000_histo.pdf")
hm<-heatmap.2(sel_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=histo_col)
legend(0.92,1.01,cex=0.6,fill=col_map_histo,ref_his,bty="n")
dev.off()

col_den<-hm$colDendrogram
all_col<-as.hclust(col_den)
all_col_tree<-cutree(all_col,k=2)
sort_sample<-rownames(hm$carpet)

inv_sam<-names(all_col_tree)[all_col_tree==1]
non_sam<-names(all_col_tree)[all_col_tree==2]

### Identify DEGs between two clusters
inv_mat<-new_mat[,inv_sam]
non_mat<-new_mat[,non_sam]

inv_mean<-apply(inv_mat,1,mean)
non_mean<-apply(non_mat,1,mean)
diff<-inv_mean-non_mean

inv_pva<-sapply(seq.int(dim(inv_mat)[[1]]), function(i) t.test(as.numeric(inv_mat[i,]),as.numeric(non_mat[i,]))$p.value)
inv_qva<-p.adjust(inv_pva,method="BH",n=length(inv_pva))

unsuper_comp<-cbind(rownames(inv_mat),gene_tab[match(rownames(inv_mat),gene_tab[,2]),1],inv_mean,non_mean,diff,inv_pva,inv_qva)
unsuper_comp<-unsuper_comp[sort.list(inv_qva),]
colnames(unsuper_comp)<-c("GeneName","EntrezID","Inv_mean","Non_mean","Difference","p-value","q-value")

### DEGs are determined based on fold change >log2(1.5) FDR<0.01
sig_comp<-unsuper_comp[which(abs(as.numeric(unsuper_comp[,5]))>0.585&as.numeric(unsuper_comp[,7])<0.01),]
sig_mat<-as.matrix(new_mat[sig_comp[,1],])
sig_mat<-sig_mat[,match(sort_sample,colnames(sig_mat))]
histo<-sam_info[match(colnames(sig_mat),sam_info$CASE),5]
ref_his<-c("AC","MP","PAP","SOL","LPA","AIS","MIA","")
col_map_histo<-c("orange","red","brown","red4","yellow","lightblue","blue","gray")
histo_col<-col_map_histo[match(histo,ref_his)]
sel_chr<-as.character(hg19[match(rownames(sig_mat),hg19$name2),3])

group<-ifelse(colnames(sig_mat)%in%inv_sam,"red","blue")

### Identify DEGs between two clusters defined by DEGs
pdf("NIHAD_FPKM_unsupervised_supervised_histo1.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=histo_col)
legend(0.92,1.01,cex=0.6,fill=c("red","blue"),c("Invasive","Non-invasive"),bty="n")
dev.off()

pdf("NIHAD_FPKM_unsupervised_supervised_group1.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=group)
legend(0.92,1.01,cex=0.6,fill=c("red","blue"),c("Invasive","Non-invasive"),bty="n")
dev.off()

col_den<-hm$colDendrogram
all_col<-as.hclust(col_den)
all_col_tree<-cutree(all_col,k=2)
group_tab<-cbind(names(all_col_tree),all_col_tree)

sort_sample<-rownames(hm$carpet)
inv_sam0<-names(all_col_tree)[all_col_tree==1]
non_sam0<-names(all_col_tree)[all_col_tree==2]

inv_mat<-new_mat[,inv_sam0]
non_mat<-new_mat[,non_sam0]

inv_mean<-apply(inv_mat,1,mean)
non_mean<-apply(non_mat,1,mean)
diff<-inv_mean-non_mean

inv_pva<-sapply(seq.int(dim(inv_mat)[[1]]), function(i) t.test(as.numeric(inv_mat[i,]),as.numeric(non_mat[i,]))$p.value)
inv_qva<-p.adjust(inv_pva,method="BH",n=length(inv_pva))

super_comp<-cbind(rownames(inv_mat),gene_tab[match(rownames(inv_mat),gene_tab[,2]),1],inv_mean,non_mean,diff,inv_pva,inv_qva)
super_comp<-super_comp[sort.list(inv_qva),]
colnames(super_comp)<-c("GeneName","EntrezID","Inv_mean","Non_mean","Difference","p-value","q-value")

### New DEGs are deteremind between invasive and non-invasive tumors
sig_super<-super_comp[which(abs(as.numeric(super_comp[,5]))>0.585&as.numeric(super_comp[,7])<0.01),]

sig_mat<-as.matrix(new_mat[sig_super[,1],])
histo<-sam_info[match(colnames(sig_mat),sam_info$CASE),5]
ref_his<-c("AC","MP","PAP","SOL","LPA","AIS","MIA","")
col_map_histo<-c("orange","red","brown","red4","yellow","lightblue","blue","gray")
histo_col<-col_map_histo[match(histo,ref_his)]

pdf("NIHAD_FPKM_unsupervised_supervised_histo2.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=histo_col)
legend(0.92,1.01,cex=0.6,fill=col_map_histo,ref_his,bty="n")
dev.off()

ref_nod<-c("","N0","N0","N1","N2","N0")
nstage<-as.character(sam_info$N.STAGE[match(colnames(sig_mat),sam_info[,1])])
nodal<-ref_nod[match(nstage,names(table(sam_info$N.STAGE)))]
ref_col<-c("gray","greenyellow","green1","darkgreen")
sam_col<-ref_col[match(nodal,names(table(nodal)))]

pdf("NIHAD_FPKM_unsupervised_supervised_nodal2.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=sam_col)
legend(0.92,1.01,cex=0.6,fill=ref_col,names(table(nodal)),bty="n")
dev.off()

ref_tum<-names(table(sam_info$T.STAGE))
tsg<-as.character(sam_info$T.STAGE[match(colnames(sig_mat),sam_info[,1])])
ref_col<-c("gray","cyan","cyan3","darkcyan")
sam_col<-ref_col[match(tsg,ref_tum)]

pdf("NIHAD_FPKM_unsupervised_supervised_tstage2.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=sam_col)
legend(0.92,1.01,cex=0.6,fill=ref_col,ref_tum,bty="n")
dev.off()

group<-ifelse(colnames(sig_mat)%in%inv_sam0,"red","blue")
pdf("NIHAD_FPKM_unsupervised_supervised_group2.pdf")
hm<-heatmap.2(sig_mat,col="bluered",scale="row",trace="none",cexRow=0.2,cexCol=0.9,hclust=function(x) hclust(x,method="complete"),distfun=function(x) as.dist(1-cor(t(x))),ColSideColors=group)
legend(0.92,1.01,cex=0.6,fill=c("red","blue"),c("Invasive","Non-invasive"),bty="n")
dev.off()

col_den<-hm$colDendrogram
all_col<-as.hclust(col_den)
all_col_tree<-cutree(all_col,k=2)
sort_sample<-rownames(hm$carpet)
inv_sam1<-names(all_col_tree)[all_col_tree==1]
non_sam1<-names(all_col_tree)[all_col_tree==2]
group_tab<-cbind(names(all_col_tree),all_col_tree)
colnames(group_tab)<-c("Sample","Group")
write.table(group_tab,"./Information/NIHAD_FPKM_invasive_noninvasive_group.txt",row.names=F,col.names=T,quote=F,sep="\t")

inv_mat<-new_mat[,inv_sam1]
non_mat<-new_mat[,non_sam1]

inv_mean<-apply(inv_mat,1,mean)
non_mean<-apply(non_mat,1,mean)
diff<-inv_mean-non_mean

inv_pva<-sapply(seq.int(dim(inv_mat)[[1]]), function(i) t.test(as.numeric(inv_mat[i,]),as.numeric(non_mat[i,]))$p.value)
inv_qva<-p.adjust(inv_pva,method="BH",n=length(inv_pva))

super_comp<-cbind(rownames(inv_mat),gene_tab[match(rownames(inv_mat),gene_tab[,2]),1],inv_mean,non_mean,diff,inv_pva,inv_qva)
super_comp<-super_comp[sort.list(inv_qva),]
colnames(super_comp)<-c("GeneName","EntrezID","Inv_mean","Non_mean","Difference","p-value","q-value")

sig_super<-super_comp[which(abs(as.numeric(super_comp[,5]))>0.585&as.numeric(super_comp[,7])<0.01),]
write.table(sig_super,"./Signatures/Invasaive_vs_noninvasive.txt",row.names=F,col.names=T,quote=F,sep="\t")
