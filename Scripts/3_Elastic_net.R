### Elastic net is used to estimate invasiveness score

library(caret)
library(glmnet)
library(gplots)
library(ggpubr)

### Our primary dataset is used for training elastic net
group<-read.table("./Information/NIHAD_FPKM_invasive_noninvasive_group.txt",header=T,sep="\t")
sam_inv<-as.character(group[which(group$Group==1),1])
sam_non<-as.character(group[which(group$Group==2),1])

sam_anno<-read.csv("./Information/NIHADKEY_UPDATED.csv")

fpkm<-read.table("./RNAseq_data/NIHAD_FPKM_normalized.txt",header=T,sep="\t")

pro_ind<-read.table("./Signatures/Invasive_vs_noninvasive.txt",header=T,sep="\t")
inv_sig<-as.character(pro_ind[which(pro_ind$Difference>0),1])
non_sig<-as.character(pro_ind[which(pro_ind$Difference<0),1])
signature<-as.character(read.table("./Signatures/Original_signature_ordered.txt")[,1])

sel<-fpkm[signature,]

sel_z<-t(apply(sel,1,scale))
colnames(sel_z)<-colnames(sel)
rownames(sel_z)<-rownames(sel)

load("./RNAseq_data/Shedden_example.RData")

sel_ori<-sel_z[match(rownames(tumr),rownames(sel_z)),]
sig_col<-ifelse(rownames(sel_ori)%in%inv_sig,"yellow","green")

### Binary invasiveness scores
inv_score<-ifelse(colnames(sel_ori)%in%sam_inv,1,0)
y.train<-as.factor(inv_score)
x.train<-t(sel_ori)

# Tune parameteres with caret and glmneti
eGrid <- expand.grid(.alpha=(0:100)*0.01, .lambda=(0:10)*0.01)
control <- trainControl(method="repeatedCV", number=5 , repeats=5)
fitM <- train(x.train, y.train,trControl=control, method="glmnet", tuneGrid=eGrid, family="binomial")

### fitM$bestTune
tumr_z<-t(apply(tumr,1,scale))
colnames(tumr_z)<-colnames(tumr)
rownames(tumr_z)<-rownames(tumr)
sig_col<-ifelse(rownames(tumr_z)%in%inv_sig,"yellow","green")

x.test<-t(tumr_z)
fitF<-glmnet(x.train, y.train, family="binomial",alpha=0.02, lambda=0.1)
yhat<-predict(fitF,s=fitF$lambda,newx=x.test,type="response")

new_yhat<-predict(fitF,s=fitF$lambda,newx=x.train,type="response")

sort_yhat<-yhat[sort.list(yhat[,1]),]

tumr_z_inv<-tumr_z[rownames(tumr_z)%in%inv_sig,]
tumr_z_non<-tumr_z[rownames(tumr_z)%in%non_sig,]
tumr_rat<-length(tumr_z_non[,1])/length(tumr_z_inv[,1])
new_score<-tumr_rat*apply(tumr_z_inv,2,sum)-apply(tumr_z_non,2,sum)

hs<-hist(sort_yhat,br=40,xlab="Relative invasiveness score",main="Shedden",cex.axis=1.2,cex.lab=1.4,cex.main=2)
hs_cnt<-hs$counts
hs_mid<-hs$mids

### Identify local minimas based on histogram
l_tok = 0
r_tok = 0
for(i in 1:50){
	left<-hs_cnt[i]
	right<-hs_cnt[51-i]
	if(i == 1){
		min<-left
		max<-right
	} else {
		if(l_tok==0){
			if(min>=left){
				min<-left
			} else {
				low<-hs_mid[i-1]
				l_tok = 1
			}
		}
		if(r_tok==0){
			if(max>=right){
				max<-right
			} else {
				print(i)
				high<-hs_mid[50-i+2]
				r_tok = 1
			}
		}
	}
}
			
# Clinical
library(survival)
sgroup<-ifelse(sort_yhat<low,1,ifelse(sort_yhat<high,2,3))
sam_sta<-cli[match(names(sgroup),cli[,1]),6]
sam_day<-cli[match(names(sgroup),cli[,1]),12]
sam_gen<-cli[match(names(sgroup),cli[,1]),3]
sam_gen<-ifelse(sam_gen=="Female",1,2)
sam_age<-cli[match(names(sgroup),cli[,1]),4]
sam_stg<-as.character(cli[match(names(sgroup),cli[,1]),9])
ref_stage<-c(1,1,2,2,2)
sam_stg<-ref_stage[match(sam_stg,names(table(sam_stg)))]

lung<-data.frame(names(sgroup),sam_day,sam_sta,sgroup,sam_gen,sam_age,sam_stg,sort_yhat)
colnames(lung)<-c("inst","time","status","group","gender","age","stage","score")
write.table(lung,"Shedden_invasiveness_score.txt",row.names=F,col.names=T,quote=F,sep="\t")

lung$time<-as.numeric(as.character(lung$time))
lung$status<-ifelse(lung$status=="Alive",1,0)
lung$SurvObj<-with(lung,Surv(time,status==0))

# Multivariate
library(survminer)
coxmul<-coxph(SurvObj ~ age + gender + stage + score, data=lung)
pdf("Shedden_score_multicox_stg12.pdf")
ggforest(coxmul,data=lung,fontsize=1)
dev.off()

lung[which(lung$time>60),3]<-1
lung[which(lung$time>60),2]<-60
lung$SurvObj<-with(lung,Surv(time,status==0))

pva1<-coxph(SurvObj ~ group, data=lung)
all<-summary(pva1)$logtest[3]
all<-ifelse(all>0.01,sprintf("%.2f",all),ifelse(all>0.001,sprintf("%.3f",all),ifelse(all>0.0001,sprintf("%.4f",all),all)))

sel_lung<-lung[which(lung$group!=2),]
pva2<-coxph(SurvObj ~ group, data=sel_lung)
two<-summary(pva2)$logtest[3]
two<-ifelse(two>0.01,sprintf("%.2f",two),ifelse(two>0.001,sprintf("%.3f",two),ifelse(two>0.0001,sprintf("%.4f",two),two)))

pdf("Shedden_score_survival.pdf")
km.by.grp<-survfit(SurvObj ~ group, data=lung, type="kaplan-meier", conf.type="none")
plot(km.by.grp,xlab="Time (month)",ylab="5-year Overall Survival",main="Shedden",lwd=1,col=c("green","gray","red"),mark.time=T,cex.axis=1.2,cex.lab=1.4,cex.main=1.8)
legend("bottomleft",col=c("green","gray","red"),lwd=2,lty=1,legend=c("Low (n=145)","Middle (n=168)","High (n=58)"),cex=1.3,bty="n")
legend("left",c("LRT p-values=",paste("All groups:",all),paste("Indolent vs Invasive:",two)),bty="n",cex=1.3)
dev.off()

