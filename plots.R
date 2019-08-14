library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(grid)
library(reshape2)

best_lambda<-sapply(res_std,function(x) which.max(rowMeans(x$auc)))
		max<-sapply(res_std,function(x) max(rowMeans(x$auc)))

best_lambda<-sapply(res,function(x) which.max(rowMeans(x$auc)))
		max<-sapply(res,function(x) max(rowMeans(x$auc)))


##
best_model<-res[[5]]


df<-melt(best_model$auc)
df<-cbind(df,melt(best_model$sensitivity)[,3])
df<-cbind(df,melt(best_model$f1)[,3])
df<-cbind(df,melt(best_model$specificity)[,3])
df<-cbind(df,melt(best_model$f1*best_model$auc)[,3])

df<-round(df,3)

names(df)<-c('lambdaRank','fold','AUC','sensitivity','f1','specificity','mix')


df_plot<-data.frame(lambdaRank=df$lambdaRank[1:100],AUC=rep(0,100),AUC_05=rep(0,100),AUC_95=rep(0,100))

df_plot$AUC<-sapply(unique(df$lambdaRank),function(x) mean(df$AUC[df$lambdaRank==x])  )
df_plot$AUC_05<-sapply(unique(df$lambdaRank),function(x) quantile(df$AUC[df$lambdaRank==x],0.05)  )
df_plot$AUC_95<-sapply(unique(df$lambdaRank),function(x) quantile(df$AUC[df$lambdaRank==x],0.95)  )

df_plot$sensitivity<-sapply(unique(df$lambdaRank),function(x) mean(df$sensitivity[df$lambdaRank==x])  )
df_plot$sensitivity_05<-sapply(unique(df$lambdaRank),function(x) quantile(df$sensitivity[df$lambdaRank==x],0.05)  )
df_plot$sensitivity_95<-sapply(unique(df$lambdaRank),function(x) quantile(df$sensitivity[df$lambdaRank==x],0.95)  )

df_plot$specificity<-sapply(unique(df$lambdaRank),function(x) mean(df$specificity[df$lambdaRank==x])  )
df_plot$specificity_05<-sapply(unique(df$lambdaRank),function(x) quantile(df$specificity[df$lambdaRank==x],0.05)  )
df_plot$specificity_95<-sapply(unique(df$lambdaRank),function(x) quantile(df$specificity[df$lambdaRank==x],0.95)  )

df_plot$f1<-sapply(unique(df$lambdaRank),function(x) mean(df$f1[df$lambdaRank==x])  )
df_plot$f1_05<-sapply(unique(df$lambdaRank),function(x) quantile(df$f1[df$lambdaRank==x],0.05)  )
df_plot$f1_95<-sapply(unique(df$lambdaRank),function(x) quantile(df$f1[df$lambdaRank==x],0.95)  )

df_plot$mix<-sapply(unique(df$lambdaRank),function(x) mean(df$mix[df$lambdaRank==x])  )
df_plot$mix_05<-sapply(unique(df$lambdaRank),function(x) quantile(df$mix[df$lambdaRank==x],0.05)  )
df_plot$mix_95<-sapply(unique(df$lambdaRank),function(x) quantile(df$mix[df$lambdaRank==x],0.95)  )
df_plot$nfeatures<-apply(best_model$nfeatures,1,max)




p <- ggplot(df_plot, aes(x = lambdaRank))
p <- p + geom_line(aes(y = AUC, colour = "AUC"),size=1)
p<-p+geom_ribbon(aes(ymin=AUC_05,ymax=AUC_95,colour='AUC'),alpha=0.08,size=0.15)
p<-p+geom_line(aes(y=sensitivity,col='sensitivity'),size=1)
p<-p+geom_ribbon(aes(ymin=sensitivity_05,ymax=sensitivity_95,col='sensitivity'),alpha=0.08,size=0.15)
p<-p+geom_line(aes(y=specificity,col='specificity'),size=1)
p<-p+geom_ribbon(aes(ymin=specificity_05,ymax=specificity_95,col='specificity'),alpha=0.08,size=0.15)
p<-p+geom_line(aes(y=f1,col='F1'),size=1)
p<-p+geom_ribbon(aes(ymin=f1_05,ymax=f1_95,col='F1'),alpha=0.08,size=0.15)
p<-p+geom_line(aes(y=mix,col='AUC*F1'),size=1)
p<-p+geom_ribbon(aes(ymin=mix_05,ymax=mix_95,col='AUC*F1'),alpha=0.08,size=0.15)
p <- p + scale_fill_brewer(palette="Dark2")
p <- p + labs(y = "Test set fit (10-fold CV)",
                x = "Lambda rank",
                colour = "Measure")
p <- p +theme_bw()+ theme(legend.position = c(0.12, 0.87))

y_nfeatures<-runif(100,-0.12,-0.04)
p<-p+ geom_text(aes(label=nfeatures,y= y_nfeatures),size=5,check_overlap = TRUE)
p<-p+ geom_text(aes(label='max n features (over CV)',x=50,y= 0.0),size=5)

ggsave('incMDD/plots/enet_nobt_fit.png', plot = p, device = 'png', dpi = 450)






bart_res<-readRDS('incMDD/test_results/results_bart_test.rds')
rf_res<-readRDS('incMDD/test_results/results_rf_test.rds')
knn_res<-readRDS('incMDD/test_results/results_knn_test.rds')
svm_res<-readRDS('incMDD/test_results/results_svm_test.rds')
logreg_res<-readRDS('incMDD/test_results/results_logreg_cov_test.rds')


hyperpar_bart<-hyperpar_list()$hyperpar_bart
hyperpar_rf<-hyperpar_list()$hyperpar_rf
hyperpar_knn<-hyperpar_list()$hyperpar_knn
hyperpar_svm<-hyperpar_list()$hyperpar_svm

hyperpar_svm$gamma<-round(hyperpar_svm$gamma,2)

models_best<-data.frame(sel=rep(c('best_f1', 'best_pr_auc', 'best_f1*pr_auc', 'union_f1', 'union_pr_auc', 'union_f1*pr_auc','VARBVS'),each=40), class=rep(c('BART','RF','KNN','SVM'),each=10), titer=1:10, PR_AUC=0,ROC_AUC=0, hyperpar_id=0, hyperpar=0)

bart_pr_mean<-apply(bart_res$roc_auc,1:2,mean)
rf_pr_mean<-apply(rf_res$roc_auc,1:2,mean)
knn_pr_mean<-apply(knn_res$roc_auc,1:2,mean)
svm_pr_mean<-apply(svm_res$roc_auc,1:2,mean)

models_best$hyperpar_id[models_best$class=='BART']<-rep(apply(bart_pr_mean,1,which.max),each=10)
models_best$hyperpar_id[models_best$class=='RF']<-rep(apply(rf_pr_mean,1,which.max),each=10)
models_best$hyperpar_id[models_best$class=='KNN']<-rep(apply(knn_pr_mean,1,which.max),each=10)
models_best$hyperpar_id[models_best$class=='SVM']<-rep(apply(svm_pr_mean,1,which.max),each=10)

for(t in 1:10){
	id<-which(models_best$titer==t & models_best$class=='BART')
	models_best$PR_AUC[id]<-sapply(1:7,function(x) bart_res$roc_auc[x,models_best$hyperpar_id[id[x]],t])

	id<-which(models_best$titer==t & models_best$class=='RF')
	models_best$PR_AUC[id]<-sapply(1:7,function(x) rf_res$roc_auc[x,models_best$hyperpar_id[id[x]],t])

	id<-which(models_best$titer==t & models_best$class=='KNN')
	models_best$PR_AUC[id]<-sapply(1:7,function(x) knn_res$roc_auc[x,models_best$hyperpar_id[id[x]],t])

	id<-which(models_best$titer==t & models_best$class=='SVM')
	models_best$PR_AUC[id]<-sapply(1:7,function(x) svm_res$roc_auc[x,models_best$hyperpar_id[id[x]],t])
}

id<-which(models_best$class=='BART')
for(i in id){
	models_best$hyperpar[i]<-paste0(paste0(names(hyperpar_bart[models_best$hyperpar_id[i],]), ':', hyperpar_bart[models_best$hyperpar_id[i],]),collapse=' , ')
}

id<-which(models_best$class=='RF')
for(i in id){
	models_best$hyperpar[i]<-paste0(paste0(names(hyperpar_rf[models_best$hyperpar_id[i],]), ':', hyperpar_rf[models_best$hyperpar_id[i],]),collapse=' , ')
}

id<-which(models_best$class=='KNN')
for(i in id){
	models_best$hyperpar[i]<-paste0('K:', hyperpar_knn[models_best$hyperpar_id[i]])
}

id<-which(models_best$class=='SVM')
for(i in id){
	models_best$hyperpar[i]<-paste0(paste0(names(hyperpar_svm[models_best$hyperpar_id[i],]), ':', hyperpar_svm[models_best$hyperpar_id[i],]),collapse=' , ')
}

logreg<-readRDS('incMDD/test_results/results_logreg_cov_test.rds')
logreg_pr_auc<-sapply(logreg,function(x) x$roc_auc)

data<-models_best[models_best$PR_AUC>0,]
data_logreg<-data.frame(sel=rep('Covariates only',10),class='Log. Regression',titer=1:10,PR_AUC=logreg_pr_auc,ROC_AUC=0,hyperpar_id=NA,hyperpar=NA)
data<-rbind(data,data_logreg)

cbp1 <- c( "#D55E00",  "#009E73",
          "#F0E442", "#0072B2","#E69F00",  "#CC79A7")

g <- ggplot(data, aes(x=sel, y=PR_AUC))
g<- g + geom_boxplot(aes(fill=factor(class))) +theme_bw()
g<- g + theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  		labs(subtitle="Box plot", 
       		title="Precision-Recall AUC Grouped by Variable Selection Method",
       		x="Variable Selection Method",
       		y="Precision-Recall AUC")

g<-g + labs(fill = "Model")+ scale_fill_manual(values = cbp1)

pr_change<-readRDS('incMDD/train_valid_test/results_varimp_pr.rds')
pr_change<-apply(pr_change,1:2,mean)
selected<-which((apply(pr_change,1,mean)<0 )  &   (rowSums(pr_change==0)<6)    )
selected_mat<-matrix(0,10,length(selected))

for(i in 1:length(selected)){
	selected_mat[,i]<-pr_change[selected[i],]
	selected_mat[,i][order(selected_mat[,i],decreasing=T)[1:2]]<-0
}

selected_mat[selected_mat==0]<-NA
colnames(selected_mat)<-names(selected)
significance<-apply(selected_mat,2,function(x) quantile(x,0.90,na.rm=T))
significant<-names(which(significance<0))
selected_mat<-selected_mat[,significant]
selected_df<-melt(selected_mat)
names(selected_df)<-c('iter','Var','PR_AUC_change')
selected_df[selected_df==0]<-NA

s <- ggplot(selected_df, aes(x=Var, y=PR_AUC_change))
s<- s + geom_boxplot() +theme_bw()
s<- s + theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  		labs(subtitle="Box plot - Decrease of Precision-Recall AUC", 
       		title="Variable Importance Plot",
       		x="Variable",
       		y="Change Precision-Recall AUC")

ggsave('incMDD/plots/var_imp_pr.png', plot = s, device = 'png', dpi = 450)

cov_cor<-numeric(length(selected))

cor_data<-cbind(annotation[,c('sex','age','smoke1','smoke2','smoke3','NK')],x_varbvs[,(significant_int)])
cormat<-cor(cor_data)

 get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

heatmap<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
	 geom_tile(color = "white")+theme_classic()+
	 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
	   midpoint = 0, limit = c(-1,1), space = "Lab", 
	   name="Pearson\nCorrelation") +
	 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
		size = 12, hjust = 1))+
	 coord_fixed()+xlab(NULL)+ylab(NULL)


ggsave('incMDD/plots/heatmap.png', plot = heatmap, device = 'png', dpi = 450)

cpg_selected<-data.frame(name=significant_int)

for(i in 1:9){
	cpg_selected$Regulatory_Feature_Group[i]<-probes$Regulatory_Feature_Group[probes$Name==cpg_selected$name[i]]
	cpg_selected$Gene_Name[i]<-probes$UCSC_RefGene_Name[probes$Name==cpg_selected$name[i]]
	cpg_selected$OpenChromatin[i]<-probes$OpenChromatin_Evidence_Count[probes$Name==cpg_selected$name[i]]
	cpg_selected$DNase_hyper[i]<-probes$DNase_Hypersensitivity_Evidence_Count[probes$Name==cpg_selected$name[i]]
	cpg_selected$TFBS[i]<-probes$TFBS_Evidence_Count[probes$Name==cpg_selected$name[i]]
}






g<-g + labs(fill = "Model")+ scale_fill_manual(values = cbp1)



data$sel_hyperpar<-paste0(data$sel,'__',data$class,':',data$hyperpar)

g2 <- ggplot(data, aes(x=sel_hyperpar, y=AUC))+geom_boxplot(aes(fill=factor(class))) + theme_linedraw() +theme(axis.text.x = element_text(angle=90, vjust=0.6))

g2<-g2 + coord_flip()

ggsave('incMDD/plots/models.png', plot = g, device = 'png', dpi = 450)
ggsave('incMDD/plots/hyperpar.png', plot = g2, device = 'png', dpi = 450)

l<-list()
for( i in 2:5){
	tmp<-readRDS(paste0('incMDD/titer/sel_titer_',i,'.rds'))
	l[[i-1]]<-tmp[[4]]
}

tmp<-intersect(l[[1]],l[[2]])
tmp<-intersect(tmp,l[[3]])
tmp<-intersect(tmp,l[[4]])

tmp

data_bn<-cbind(annotation[,c('sex','smoke12','bmi','age')],meth[,tmp[-1]])

data_bn$sex<-as.factor(data_bn$sex)
data_bn$smoke12<-as.factor(data_bn$smoke12)


structure <- tabu(data_bn, score = "bic-cg",maxp=1)

png('incMDD/plots/bn.png', dpi = 450)
plot(structure)
dev.off()

res<-bn.fit(net, data_bn)

