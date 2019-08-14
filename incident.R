
library(gridExtra)
library(ggplot2)
library(pROC)
library(biglasso)
library(BART)
source('functions.R')
source('classifiers.R')

library(grid)



meth<-readRDS('meth.rds')


standardise<-function(data){
	return((data-mean(data))/sd(data))
}


annotation<-readRDS('incMDD/annotationF.rds')


#annotation$smoke[is.na(annotation$smoke)]<-4
#annotation$smoke12<-1*(annotation$smoke < 3)
#annotation<-annotation[annotation$unrelated==1,]
#incidence<-readRDS('incMDD/both_waves_dnam_incidence')

#annotation$MDDinc<-incidence$depression[match(annotation$Sample_Sentrix_ID,incidence$Sample_Sentrix_ID)]
#annotation$bmi[which(is.na(annotation$bmi),arr.ind=T)]<-mean(annotation$bmi,na.rm=T)
#saveRDS(annotation,'incMDD/annotationF.rds')

annotation$smoke1<-0
annotation$smoke2<-0
annotation$smoke3<-0

annotation$smoke[is.na(annotation$smoke)]<-4

annotation$smoke1[(annotation$smoke)==1]<-1
annotation$smoke2[(annotation$smoke)==2]<-1
annotation$smoke3[(annotation$smoke)==3]<-1

idx<-which(!is.na(annotation$MDDinc) &  annotation$unrelated==1)

#meth<-meth[idx,]
#meth_std<-apply(meth,2,standardise)
#meth<-NA
#saveRDS(meth_std,'incMDD/meth_std.rds')


annotation<-annotation[idx,]

y<-annotation$MDDinc

y<-as.factor(y)


covariates<-c('sex','age','Bcell','CD4T','CD8T','Gran','NK','Mono','bmi','smoke1','smoke2','smoke3')

meth_std_cov<-as.data.frame(cbind(annotation[,covariates],readRDS('incMDD/meth_std.rds')))

meth_std_cov_bm<-as.big.matrix(meth_std_cov)


##### TAKE UNION OF SELVARS on DIFFERENT ALPHA LEVELS

### The bigLASSO package offers the user a cross-validation hyperparameter autotunning function. However, the metric choice is only limited to misclassification error (MSE and MAPE for regression). We prefer to evaluate the goodness of fit with area under reciever operator curve and F1 score - harmonic mean of precision and recall. Therefore, we construct the cross-validation step manually, measuring the AUC and F1 score.


### For given level of alpha, too low lambda values results in unstable convergence and costly computations. The bigLASSO package can precompute a feasible reagion of lambdas at given alpha level when it is given number of lambda (nlambda). At one hand, this is very useful for us as users because we don't have to calibrate the min/max lambda, which can be very sensitive with high dimensional data, but when we want to constructy the k-fold crossvalidation on our own the package would use non-consistent lambda values at each k. Therefore we use the inner package function to find the feasible range of lambdas at different alpha levels


penalty.factor<-numeric(dim(meth_std_cov)[2])
penalty.factor<-penalty.factor+1
penalty.factor[c(1,2,7,10)]<-0.001



covariates<-c('sex','age','Bcell','CD4T','CD8T','Gran','NK','Mono','bmi','smoke12')
meth_cov<-cbind(annotation[,covariates],meth_std)
meth_cov_bm<-as.big.matrix(meth_cov)

bart_test_list<-list()
rf_test_list<-list()
knn_test_list<-list()
svm_test_list<-list()




results_knn_test<-list(pr_auc=array(0,dim=c(7,50,10)),roc_auc=array(0,dim=c(7,50,10)))
results_bart_test<-list(pr_auc=array(0,dim=c(7,36,10)),roc_auc=array(0,dim=c(7,36,10)))
results_rf_test<-list(pr_auc=array(0,dim=c(7,180,10)),roc_auc=array(0,dim=c(7,180,10)))
results_svm_test<-list(pr_auc=array(0,dim=c(7,130,10)),roc_auc=array(0,dim=c(7,130,10)))


results_logreg_cov_test<-list()


for(titer in 3:10){
	test<-sample(1:length(y),round(length(y)/10))
	train_valid<-(1:length(y))[-test]
	model_enet_CV_bootstrap<-list()
	for( i in 1:9){
		print(' ')
		print(paste0('ALPHA=' ,i/10))
		print(' ')
		res<-elastic_net_CV_bootstrap(x=meth_std_cov_bm, y=y, kfolds=10, train_valid=train_valid, nlambda=100, alpha=(i/10), ncores=15)
		model_enet_CV_bootstrap[[i]]<-res
		saveRDS(model_enet_CV_bootstrap,'incMDD/titer/model_enet_CV_bootstrap.rds')
	}
	sel_1<-var_selection_best(model_enet_CV_bootstrap,method='f1',titer)
	sel_2<-var_selection_best(model_enet_CV_bootstrap,method='pr_auc',titer)
	sel_3<-var_selection_best(model_enet_CV_bootstrap,method='pr_auc_f1',titer) 
	sel_4<-var_selection_union(model_enet_CV_bootstrap,method='f1')
	sel_5<-var_selection_union(model_enet_CV_bootstrap,method='pr_auc')
	sel_6<-var_selection_union(model_enet_CV_bootstrap,method='pr_auc_f1') 
	

	x_varbvs<-readRDS('incMDD/x_varbvs.rds')
	fit.varbvs<-varbvs(X=x_varbvs[train_valid,],y=as.numeric(as.character(y[train_valid])),Z=annotation[train_valid,c('sex','age','smoke12','NK')],family='binomial')
	sel_varbvs<-names(fit.varbvs$pip[fit.varbvs$pip>0.95])
	sel_varbvs<-c('sex','age','smoke12','NK',sel_varbvs)
	sel<-list(sel_1,sel_2,sel_3,sel_4,sel_5,sel_6,sel_varbvs)
	
	saveRDS(sel,paste0('incMDD/titer/sel_titer_',titer,'.rds') )
	bart_test_list[[titer]]<-list()
	rf_test_list[[titer]]<-list()
	knn_test_list[[titer]]<-list()
	svm_test_list[[titer]]<-list()


	results_logreg_cov_test[[titer]]<-logreg_cov(annotation,train_valid,test)
	for(l in 1:length(sel) ){
		x<-meth_std_cov[,sel[[l]]]
		x$age<-standardise(x$age)
		#bart_hyperpar_sel<-bart_CV(x,y,train_valid) 
		bart_fit<-bart_hyperpar_test(x,y,train_valid,test,l)
		#bart_test_best<-bart_test(x,y,train_valid,test,bart_fit$hyperpar)
		#bart_test_list[[titer]][[l]]<-bart_test_best$pred

		#rf_sel<-rf_CV(x,y,train_valid) 
		rf_fit<-rf_hyperpar_test(x,y,train_valid,test,l)
		#rf_test_best<-rf_test(x,y,train_valid,test,rf_fit$hyperpar)
		#rf_test_list[[titer]][[l]]<-rf_test_best$pred

		#knn_CV<-knn_CV(x,y,train_valid)
		knn_fit<-knn_hyperpar_test(x,y,train_valid,test)
		#knn_test_best<-knn_test(x,y,train_valid,test,knn_fit$hyperpar)	
		#knn_test_list[[titer]][[l]]<-knn_test_best$pred

		#svm_CV<-svm_CV(x,y,train_valid)
		svm_fit<-svm_hyperpar_test(x,y,train_valid,test,l)
		#svm_test_best<-svm_test(x,y,train_valid,test,svm_fit$hyperpar)
		#svm_test_list[[titer]][[l]]<-svm_test_best$pred

		results_bart_test$roc_auc[l,,titer]<-bart_fit$roc_auc
		results_rf_test$roc_auc[l,,titer]<-rf_fit$roc_auc
		results_knn_test$roc_auc[l,,titer]<-knn_fit$roc_auc
		results_svm_test$roc_auc[l,,titer]<-svm_fit$roc_auc

		results_bart_test$pr_auc[l,,titer]<-bart_fit$pr_auc
		results_rf_test$pr_auc[l,,titer]<-rf_fit$pr_auc
		results_knn_test$pr_auc[l,,titer]<-knn_fit$pr_auc
		results_svm_test$pr_auc[l,,titer]<-svm_fit$pr_auc



		saveRDS(results_bart_test,'incMDD/test_results/results_bart_test.rds')
		saveRDS(results_rf_test,'incMDD/test_results/results_rf_test.rds')
		saveRDS(results_knn_test,'incMDD/test_results/results_knn_test.rds')
		saveRDS(results_svm_test,'incMDD/test_results/results_svm_test.rds')
		saveRDS(results_logreg_cov_test,'incMDD/test_results/results_logreg_cov_test.rds')

		saveRDS(bart_test_list,'incMDD/test_results/results_bart_test_prob.rds')
		saveRDS(svm_test_list,'incMDD/test_results/results_svm_test_prob.rds')
		saveRDS(knn_test_list,'incMDD/test_results/results_knn_test_prob.rds')
		saveRDS(rf_test_list,'incMDD/test_results/results_rf_test_prob.rds')
	}

}

	






#################################



################################


###10-1114

set.seed(2222)

x_cand<-readRDS('incMDD/x_varbvs.rds')

annotation<-readRDS('incMDD/annotationF.rds')
idx<-which(!is.na(annotation$MDDinc) &  annotation$unrelated==1)
annotation<-annotation[idx,]
annotation$smoke1<-0
annotation$smoke2<-0
annotation$smoke3<-0

annotation$smoke[is.na(annotation$smoke)]<-4
annotation$smoke1[(annotation$smoke)==1]<-1
annotation$smoke2[(annotation$smoke)==2]<-1
annotation$smoke3[(annotation$smoke)==3]<-1

y<-annotation$MDDinc
y<-as.factor(y)
covariates<-c('sex','age','NK','bmi','smoke1','smoke3')

x_cov<-cbind(annotation[,covariates],x_cand)
x_cov<-as.matrix(x_cov)



pr_auc_valid<-array(0,dim=c(36,10))
roc_auc_valid<-array(0,dim=c(36,10))

pr_auc_logreg<-numeric(10)
roc_auc_logreg<-numeric(10)

pr_auc_test<-array(0,dim=c(10,10))
roc_auc_test<-array(0,dim=c(10,10))

pr_change_array<-array(0,dim=c(dim(x_cov)[2],10,10))
roc_change_array<-array(0,dim=c(dim(x_cov)[2],10,10))
rownames(pr_change_array)<-colnames(x_cov)
rownames(roc_change_array)<-colnames(x_cov)

varbvs_selection<-list()
	hyperpar<-hyperpar_list()$hyperpar_bart
	folds<-stratified_crossvalidation(y,1:length(y),10)
	iter_vec<-c(1:10,1)
	for(iter in 1:10){
		train<-which((folds!=iter)&(folds!=(iter_vec[iter+1])))
		fit.varbvs<-varbvs(X=x_cov[train,], y=as.numeric(as.character(y[train])), Z=NULL, family='binomial')
		sel_varbvs<-names(fit.varbvs$pip[fit.varbvs$pip>0.95])
		varbvs_selection[[iter]]<-sel_varbvs
	}
	for(iter in 1:10){
		valid<-which(folds==iter_vec[iter+1])
		x<-x_cov[,varbvs_selection[[iter]] ]
		train<-which((folds!=iter)&(folds!=(iter_vec[iter+1])))
		x<-cbind(annotation[,covariates],x)	
		res_valid<-foreach(h=1:36) %dopar% {

			print('----------------------')
			print(paste0('VALIDATION:',iter))
			print('----------------------')

			model<-lbart(x[train,], as.numeric(as.character(y[train])), ndpost=1000, power=hyperpar[h,2], ntree=hyperpar[h,1], k=hyperpar[h,3], nskip=1000,printevery=1000L)
			pred<-predict(model,x[valid,])$prob.test.mean
			roc_auc_valid<-(measure_fit(pred,as.numeric(as.character(y[valid]))))$roc_auc
			pr_auc_valid<-(measure_fit(pred,as.numeric(as.character(y[valid]))))$pr_auc
			list(roc_auc_valid=roc_auc_valid, pr_auc_valid=pr_auc_valid )
		}
		for(h in 1:36){
			pr_auc_valid[h,iter]<-(res_valid[[h]])$pr_auc
			roc_auc_valid[h,iter]<-(res_valid[[h]])$roc_auc
		}
	}
	best_hyperpar<-which.max(rowMeans(pr_auc_valid))

	for(iter in 1:10){
		test<-which(folds==(iter))
		train<-which(folds!=(iter))
		x<-x_cov[,varbvs_selection[[iter]] ]

	
		glm.fit <- glm(y ~.,  data = cbind(y,as.data.frame(x)), family = binomial, subset = train)
		glm.probs <- predict(glm.fit, newdata = as.data.frame(x[test,]), type = "response")

		pr_auc_logreg[iter]<-measure_fit(glm.probs,y[test])$pr_auc
		roc_auc_logreg[iter]<-measure_fit(glm.probs,y[test])$roc_auc
		
		x<-cbind(annotation[,covariates],x)
		res_test<-foreach(K=1:10) %dopar% {
			print('----------------------')
			print(paste0('TEST:',iter,'     K:',K))
			print('----------------------')

			res<-bart_permute(x[train,], y[train], x[test,], y[test], hyperpar[best_hyperpar,])
			list(pr_change=res$permute_pr_auc, roc_change=res$permute_roc_auc, pr_auc=res$pr_auc_fit, roc_auc=res$roc_auc_fit )
		}
		for(K in 1:10){
			pr_change_array[match(colnames(x),colnames(x_cov)),iter,K]<-(res_test[[K]])$pr_change
			roc_change_array[match(colnames(x),colnames(x_cov)),iter,K]<-(res_test[[K]])$roc_change
			pr_auc_test[iter,K]<-(res_test[[K]])$pr_auc
			roc_auc_test[iter,K]<-(res_test[[K]])$roc_auc

		}

		saveRDS(pr_change_array,'incMDD/train_valid_test/results_11CV_varimp_pr.rds')
		saveRDS(roc_change_array,'incMDD/train_valid_test/results_11CV_varimp_roc.rds')
		saveRDS(pr_auc_test,'incMDD/train_valid_test/results_11CV_pr_auc_test.rds')
		saveRDS(roc_auc_test,'incMDD/train_valid_test/results_11CV_roc_auc_test.rds')
		saveRDS(pr_auc_logreg,'incMDD/train_valid_test/results_11CV_pr_auc_logreg.rds')
		saveRDS(roc_auc_logreg,'incMDD/train_valid_test/results_11CV_roc_auc_logreg.rds')

	}


################

################	Goodness of a fit results

################

bart_logreg_diff<-rowMeans(pr_auc_test)-pr_auc_logreg
bart_pr_test<-rowMeans(pr_auc_test)
logreg_pr_test<-pr_auc_logreg
diff_df<-data.frame(matrix(c(bart_pr_test,logreg_pr_test),10,2))
names(diff_df)<-c('BART','LogRegression')
diff_df<-melt(diff_df,value.name='PR_AUC',variable.name='Model')

cbp1 <- c(  "#E69F00","#D55E00")
diff_plot<-ggplot(diff_df)+geom_jitter(aes(x=Model,col=Model,y=PR_AUC),width=0.1,size=3)+theme_bw()+ theme(legend.position="top")
diff_plot<-diff_plot+labs(title="Precision-Recall AUC")
diff_plot<-diff_plot+scale_colour_manual(values = cbp1) +theme(plot.title = element_text( size=12))


diff_boxplot<-ggplot()+stat_boxplot(aes(x='10 CV folds',y=bart_logreg_diff,middle=mean(bart_logreg_diff)),width=0.3,fill='grey')
diff_boxplot<-diff_boxplot+geom_jitter(aes(x='10 CV folds',y=bart_logreg_diff),width=0.15,size=3,col="#0072B2")+theme_bw()
diff_boxplot<-diff_boxplot+ylab('BART PR-AUC - LogReg PR-AUC')+xlab(NULL)
diff_boxplot<-diff_boxplot+geom_hline(aes(yintercept=0),col='red',linetype=2,size=1)
diff_boxplot<-diff_boxplot+labs(title='Individual CV Fold Difference')+theme(plot.title = element_text( size=12))

bart_logreg_plot<-grid.arrange(diff_plot, diff_boxplot,ncol=2, top=textGrob("Classification Performance of BART and Logistic Regression of 10-fold CV", gp=gpar(fontsize=15,font=8)))

ggsave('incMDD/plots/bart_logreg_plot.png', plot = bart_logreg_plot, device = 'png', dpi = 450)



iter<-9
test<-which(folds==(iter+10))
train<-which(folds!=(iter+10))
x<-x_cov[,varbvs_selection[[iter]] ]

glm.fit <- glm(y ~.,  data = cbind(y,as.data.frame(x)), family = binomial, subset = train)
glm.probs <- predict(glm.fit, newdata = as.data.frame(x[test,]), type = "response")
	
x<-cbind(annotation[,covariates],x)
model <- lbart(x[train,], as.numeric(as.character(y[train])), ndpost=2000, power=hyperpar[best_hyperpar,2], ntree=hyperpar[best_hyperpar,1], k=hyperpar[best_hyperpar,3], nskip=2000,printevery=1000L)

pred<-predict(model,x[test,])$prob.test.mean
roc<-roc.curve(scores.class0 = pred[y[test] == 1], scores.class1=pred[y[test]==0], curve=TRUE)
pr<-pr.curve(scores.class0 = pred[y[test] == 1], scores.class1=pred[y[test]==0], curve=TRUE)
roc_logreg<-roc.curve(scores.class0=glm.probs[y[test]==1], scores.class1=glm.probs[y[test]==0], curve=TRUE)
pr_logreg<-pr.curve(scores.class0=glm.probs[y[test]==1], scores.class1=glm.probs[y[test]==0], curve=TRUE)

roc_df<-as.data.frame(rbind(roc$curve,roc_logreg$curve))
roc_df[,3]<-'LogRegression'
roc_df[1:nrow(roc$curve),3]<-'BART'
names(roc_df)<-c('FPR','Recall','Model')
cols<-c("#E69F00","#D55E00")

pr_df<-as.data.frame(rbind(pr$curve,pr_logreg$curve))
pr_df[,3]<-'LogRegression'
pr_df[1:nrow(pr$curve),3]<-'BART'
names(pr_df)<-c('Recall','Precision','Model')

roc_plot<-ggplot(roc_df)+geom_line(aes(x=FPR,y=Recall,col=Model),size=1.5)+ coord_fixed()+theme_linedraw()+theme(legend.position="top")
roc_plot<-roc_plot+scale_colour_manual(values = cols)+labs(title='ROC-AUC curve')
roc_plot

pr_plot<-ggplot(pr_df)+geom_line(aes(x=Recall,y=Precision,col=Model),size=1.5)+ coord_fixed()+theme_linedraw()+theme(legend.position="top")
pr_plot<-pr_plot+scale_colour_manual(values = cols)+labs(title='PR-AUC curve')
pr_plot


auc_both_plot<-grid.arrange(roc_plot, pr_plot,ncol=2)
ggsave('incMDD/plots/auc_both_plot.png', plot = auc_both_plot, device = 'png', dpi = 450)


################

################	VARIABBLE IMPORTANCE

################


pr_change<-readRDS('incMDD/train_valid_test/results_varimp_pr.rds')
pr_change<-apply(pr_change,1:2,mean)
selected<-which((apply(pr_change,1,mean)<0.005 )  &   (rowSums(pr_change==0)<6)    )
selected_mat<-matrix(0,10,length(selected))

for(i in 1:length(selected)){
	selected_mat[,i]<-pr_change[selected[i],]
	selected_mat[,i][order(selected_mat[,i],decreasing=T)[1]]<-mean(selected_mat[,i])
}

selected_mat[selected_mat==0]<-NA
colnames(selected_mat)<-names(selected)
significance<-apply(selected_mat,2,function(x) quantile(x,0.70,na.rm=T))
significant<-names(which(significance< 0))

selected_mat<-selected_mat[,significant]
selected_mat<-selected_mat[,order(colMeans(selected_mat,na.rm=T))]
selected_df<-melt(selected_mat)
names(selected_df)<-c('iter','Var','PR_AUC_change')
selected_df[selected_df==0]<-NA

selected_df$Var_type[which(substr(selected_df[,2], 1, 2)=='cg')]<-'CpG'
selected_df$Var_type[which(substr(selected_df[,2], 1, 2)!='cg')]<-'Covariate'

s <- ggplot(selected_df, aes(x=Var, y=PR_AUC_change,col=Var_type))
s<- s + geom_boxplot(coef = 4) +theme_bw()+geom_jitter()+geom_hline(yintercept=0)
s<- s + theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  		labs(subtitle="Box plot - Decrease of Precision-Recall AUC", 
       		title="Variable Importance Plot",
       		x="Variable",
       		y="Change Precision-Recall AUC")
s


ggsave('incMDD/plots/var_imp_pr.png', plot = s, device = 'png', dpi = 450)


s_df<-list(s1_df,s2_df,s3_df)

s_all<-Reduce(union,list(s1,s2,s3))
selection_rep<-data.frame(Var=s_all,rep1=NA,rep2=NA,rep3=NA)

for(i in 1:3){
	df<-s_df[[i]]
	agg<-aggregate(df$PR_AUC_change, by=list(df$Var), FUN=mean, na.rm=TRUE)

	selection_rep[,i+1]<-agg$x[match(selection_rep$Var,agg$Group.1)]
}

sign_high<-which(rowSums(!is.na(selection_rep))>2)
rowMeans(selection_rep[sign_high,])
sign_high_df<-data.frame(Var=selection_rep[sign_high,1],PR_AUC_change=rowMeans(selection_rep[sign_high,2:4],na.rm=T) )

sign_high_df<-sign_high_df[order(sign_high_df[,2]),1]

sign_high_df$Var <- factor(sign_high_df$Var, levels = sign_high_df[order(sign_high_df[,2]),1])

sign_high_df$Var_type[which(substr(sign_high_df[,1], 1, 2)=='cg')]<-'CpG'
sign_high_df$Var_type[which(substr(sign_high_df[,1], 1, 2)!='cg')]<-'Covariate'


p<-ggplot(data=sign_high_df, aes(x=Var, y=abs(PR_AUC_change),fill=Var_type)) +
  geom_bar(stat="identity")+
  theme_minimal()+ theme(axis.text.x = element_text(angle=65, vjust=0.6))+
  ylab('negative PR-AUC change')+xlab('Variable')
p



#########################  HEATMAP CORRELATION PLOT


u_cpg<-as.character(sign_high_df$Var[sign_high_df$Var_type=='CpG'])
variables<-c('age','sex','smoke1','smoke2','smoke3','bmi','NK')
data_covmat<-cbind(annotation[,variables],readRDS('incMDD/meth_high_sign.rds'))


cormat <- round(cor(data_covmat),2)
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
corr_plot<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+coord_fixed()
corr_plot<-corr_plot + geom_text(aes(Var2, Var1, label = value), color = "black", size = 3.4) +
theme(axis.title.x = element_blank(),  axis.title.y = element_blank(),  panel.grid.major = element_blank(),
  panel.border = element_blank(),  panel.background = element_blank(),  axis.ticks = element_blank(),
  legend.justification = c(1, 0),   legend.position = c(0.6, 0.7),  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
corr_plot

ggsave('incMDD/plots/corr_plot.png', plot = corr_plot, device = 'png', dpi = 450)





setwd('../')
probes<-readRDS('Daniel/EPIC_AnnotationObject_df.rds')
setwd('Imrich/')

u_cpg<-sign_high_df$Var[sign_high_df$Var_type=='CpG']
cpg_selected<-data.frame(name=u_cpg)

for(i in 1:length(u_cpg)){
	cpg_selected$chr[i]<-probes$chr[probes$Name==cpg_selected$name[i]]
	cpg_selected$Regulatory_Feature_Group[i]<-probes$Regulatory_Feature_Group[probes$Name==cpg_selected$name[i]]
	cpg_selected$Gene_Name[i]<-probes$UCSC_RefGene_Name[probes$Name==cpg_selected$name[i]]
	cpg_selected$OpenChromatin[i]<-probes$OpenChromatin_Evidence_Count[probes$Name==cpg_selected$name[i]]
	cpg_selected$DNase_hyper[i]<-probes$DNase_Hypersensitivity_Evidence_Count[probes$Name==cpg_selected$name[i]]
}





meth<-readRDS('meth.rds')
find_correlated<-function(cpg,cpg_list,meth){
	corm<-cor(meth[,c(cpg,cpg_list)])
	return(corm[1,])
}

pg_list<-c('cg03021120' ,'cg16435316', 'cg16743005', 'cg16834187')
cpg<-'cg06576748' 

tmp<-find_correlated(cpg,cpg_list,meth)






