library("GSEABase")
library("gplots")

### Load FPKM normalizaed values and separate them into two subgroups
fpkm<-read.table("./RNAseq_data/NIHAD_FPKM_normalized.txt",row.names=1,header=T,sep="\t")
univ_gene<-rownames(fpkm)
univ<-length(univ_gene)


### Pro-invasive and indolence signature genes
signature<-read.table("./Signatures/Invasive_vs_noninvasive.txt",header=T,sep="\t")
inv<-as.character(signature[which(signature$Difference>0),1])
non<-as.character(signature[which(signature$Difference<0),1])


### Fisher's exact test against HALLMARK genesets
hall<-getGmt(con="./Msigdb_database/h.all.v6.0.symbols.gmt",geneIdType=SymbolIdentifier(),collectionType=BroadCollection(category="h"))
up_len<-length(inv)
dn_len<-length(non)
upfinal<-matrix(NA,nr=length(hall),nc=6)
colnames(upfinal)<-c("Geneset","Targets","Overlap","Oddratio","p-value","Overlap_gene")
dnfinal<-upfinal
pva_cut<-0.05/length(hall)
for(i in 1:length(hall)){
        gs<-hall[i]
        gsname<-names(gs)
        gs_list<-unique(unlist(geneIds(gs)))
        gs_list<-gs_list[gs_list%in%univ_gene]
        len_gs<-length(gs_list)
        overlap<-intersect(gs_list,inv)
        over<-length(overlap)
        tab<-matrix(c(univ-up_len-len_gs+over,up_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        upfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
        overlap<-intersect(gs_list,non)
        over<-length(overlap)
        tab<-matrix(c(univ-dn_len-len_gs+over,dn_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        dnfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
}
inv_final<-upfinal[sort.list(as.numeric(upfinal[,5])),]
inv_final<-inv_final[which(as.numeric(inv_final[,5])<pva_cut),]
write.table(inv_final,"Pro_invasive_enrichment_hallmark.txt",row.names=F,col.names=T,quote=F,sep="\t")
non_final<-dnfinal[sort.list(as.numeric(dnfinal[,5])),]
non_final<-t(as.matrix(non_final[which(as.numeric(non_final[,5])<pva_cut),]))
write.table(non_final,"Indolent_enrichment_hallmark.txt",row.names=F,col.names=T,quote=F,sep="\t")


### Fisher's exact test against C2
c2<-getGmt(con="./Msigdb_database/c2.all.v6.0.symbols.gmt",geneIdType=SymbolIdentifier(),collectionType=BroadCollection(category="c2"))
up_len<-length(inv)
dn_len<-length(non)
upfinal<-matrix(NA,nr=length(c2),nc=6)
colnames(upfinal)<-c("Geneset","Targets","Overlap","Oddratio","p-value","Overlap_gene")
dnfinal<-upfinal
pva_cut<-0.05/length(c2)
for(i in 1:length(c2)){
        gs<-c2[i]
        gsname<-names(gs)
        gs_list<-unique(unlist(geneIds(gs)))
        gs_list<-gs_list[gs_list%in%univ_gene]
        len_gs<-length(gs_list)
        overlap<-intersect(gs_list,inv)
        over<-length(overlap)
        tab<-matrix(c(univ-up_len-len_gs+over,up_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        upfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
        overlap<-intersect(gs_list,non)
        over<-length(overlap)
        tab<-matrix(c(univ-dn_len-len_gs+over,dn_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        dnfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
}
inv_final<-upfinal[sort.list(as.numeric(upfinal[,5])),]
inv_final<-inv_final[which(as.numeric(inv_final[,5])<pva_cut),]
non_final<-dnfinal[sort.list(as.numeric(dnfinal[,5])),]
non_final<-non_final[which(as.numeric(non_final[,5])<pva_cut),]

index<-c("LUNG","INVAS","MIGRA","EMT","EPITHELIAL","METAST","ANGIO")
sel_ind<-c()
sel<-unlist(sapply(seq.int(length(index)), function(i) c(sel_ind,grep(index[i],inv_final[,1]))))
inv_sel<-inv_final[sel,]
inv_sel<-inv_sel[sort.list(as.numeric(inv_sel[,5])),]
write.table(inv_sel,"Pro_invasive_enrichment_C2_sel.txt",row.names=F,col.names=T,quote=F,sep="\t")
sel_ind<-c()
sel<-unlist(sapply(seq.int(length(index)), function(i) c(sel_ind,grep(index[i],non_final[,1]))))
non_sel<-non_final[sel,]
non_sel<-non_sel[sort.list(as.numeric(non_sel[,5])),]
write.table(non_sel,"Indolent_enrichment_C2_sel.txt",row.names=F,col.names=T,quote=F,sep="\t")


### Fisher's exact test against GO database
c5<-getGmt(con="./Msigdb_database/c5.all.v6.0.symbols.gmt",geneIdType=SymbolIdentifier(),collectionType=BroadCollection(category="c5"))
up_len<-length(inv)
dn_len<-length(non)
upfinal<-matrix(NA,nr=length(c5),nc=6)
colnames(upfinal)<-c("Geneset","Targets","Overlap","Oddratio","p-value","Overlap_gene")
dnfinal<-upfinal
pva_cut<-0.05/length(c5)
for(i in 1:length(c5)){
        gs<-c5[i]
        gsname<-names(gs)
        gs_list<-unique(unlist(geneIds(gs)))
        gs_list<-gs_list[gs_list%in%univ_gene]
        len_gs<-length(gs_list)
        overlap<-intersect(gs_list,inv)
        over<-length(overlap)
        tab<-matrix(c(univ-up_len-len_gs+over,up_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        upfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
        overlap<-intersect(gs_list,non)
        over<-length(overlap)
        tab<-matrix(c(univ-dn_len-len_gs+over,dn_len-over,len_gs-over,over),nr=2)
        fish_test<-fisher.test(tab,alternative="greater")
        pva<-fish_test$p.value
        odd<-fish_test$estimate
        genelist<-paste(overlap,collapse=",")
        dnfinal[i,]<-c(gsname,len_gs,over,odd,pva,genelist)
}
inv_final<-upfinal[sort.list(as.numeric(upfinal[,5])),]
inv_final<-inv_final[which(as.numeric(inv_final[,5])<pva_cut),]
write.table(inv_final,"Pro_invasive_enrichment_C5.txt",row.names=F,col.names=T,quote=F,sep="\t")
non_final<-dnfinal[sort.list(as.numeric(dnfinal[,5])),]
non_final<-non_final[which(as.numeric(non_final[,5])<pva_cut),]
write.table(non_final,"Indolent_enrichment_C5.txt",row.names=F,col.names=T,quote=F,sep="\t")
com_tab<-cbind(upfinal[,1:5],dnfinal[,3:5])


## Tumor suppressor genes
tsg<-read.table("./Signatures/Human_TSGs.txt",header=T,sep="\t")
tsg<-tsg[tsg$GeneSymbol%in%univ_gene,]
tsg_name<-as.character(tsg$GeneSymbol)
tsg_len<-length(tsg_name)

overlap<-intersect(tsg_name,inv)
over<-length(overlap)
tab<-matrix(c(univ-up_len-tsg_len+over,up_len-over,tsg_len-over,over),nr=2)
fish_test<-fisher.test(tab,alternative="greater")
pva<-fish_test$p.value
odd<-fish_test$estimate
print(paste(paste(sort(overlap),collapse=","),odd,pva))

final<-c("Invasive",over,pva,odd,paste(sort(overlap),collapse=","))
overlap<-intersect(tsg_name,non)
over<-length(overlap)
tab<-matrix(c(univ-dn_len-tsg_len+over,dn_len-over,tsg_len-over,over),nr=2)
fish_test<-fisher.test(tab,alternative="greater")
pva<-fish_test$p.value
odd<-fish_test$estimate
print(paste(paste(sort(overlap),collapse=","),odd,pva))

final<-rbind(final,c("Non-invasive",over,pva,odd,paste(sort(overlap),collapse=",")))
colnames(final)<-c("Geneset","# of overlap","FET pva","Odds ratio","genelist")
write.table(final,"Indolent_TSGs.txt",row.names=F,col.names=T,quote=F,sep="\t")

