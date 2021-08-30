msigdb_fish <- function(query, gene_univ, dataset)
{
    set<-getGmt(con=paste("./Msigdb_database/",dataset,".all.v6.0.symbols.gmt",sep=""),geneIdType=SymbolIdentifier(),collectionType=BroadCollection(category=dataset))
    univ<-length(gene_univ)
    gs1<-query
    len1<-length(unique(query))
    res<-matrix(NA,nr=length(set),nc=7)
    colnames(res)<-c("Geneset","Size","Overlap","OddRatio","p-value","q-value","Genelist")
    for(i in 1:length(set)){
        gs<-set[i]
        gsname<-names(gs)
        gs2<-unique(unlist(geneIds(gs)))
        gs2<-gs2[gs2%in%gene_univ]
        len2<-length(gs2)
        overlap<-sort(intersect(gs1,gs2))
	over<-length(overlap)
        tab<-matrix(c(univ-len1-len2+over,len1-over,len2-over,over),nr=2)
	fish_test<-fisher.test(tab,alternative="greater")
	pva<-fish_test$p.value
        odd<-fish_test$estimate
        res[i,c(1:5,7)]<-c(gsname,len2,over,odd,pva,paste(overlap,collapse=", "))
    }
    res[,6]<-p.adjust(as.numeric(res[,5]),method="BH",n=length(set))
    return(res)	
}

