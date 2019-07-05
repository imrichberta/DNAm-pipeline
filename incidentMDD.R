
library(gridExtra)
library(ggplot2)
library(pROC)
library(biglasso)
library(BART)
source('functions.R')
source('classifiers.R')



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

idx<-which(!is.na(annotation$MDDinc) &  annotation$unrelated==1)

meth<-meth[idx,]
meth_std<-apply(meth,2,standardise)
meth<-NA
#saveRDS(meth_std,'incMDD/meth_std.rds')
meth_std<-as.data.frame(meth_std)
meth_std_bm<-as.big.matrix(meth_std)

annotation<-annotation[idx,]

y<-annotation$MDDinc

y<-as.factor(y)



##### TAKE UNION OF SELVARS on DIFFERENT ALPHA LEVELS

### The bigLASSO package offers the user a cross-validation hyperparameter autotunning function. However, the metric choice is only limited to misclassification error (MSE and MAPE for regression). We prefer to evaluate the goodness of fit with area under reciever operator curve and F1 score - harmonic mean of precision and recall. Therefore, we construct the cross-validation step manually, measuring the AUC and F1 score.


### For given level of alpha, too low lambda values results in unstable convergence and costly computations. The bigLASSO package can precompute a feasible reagion of lambdas at given alpha level when it is given number of lambda (nlambda). At one hand, this is very useful for us as users because we don't have to calibrate the min/max lambda, which can be very sensitive with high dimensional data, but when we want to constructy the k-fold crossvalidation on our own the package would use non-consistent lambda values at each k. Therefore we use the inner package function to find the feasible range of lambdas at different alpha levels

lambda_search<-list()
for(i in 1:8){
	res<-elastic_net(x=meth_std_bm,y=y, train_valid=(1:length(y)),nlambda=500
	,alpha=(i/10+0.1),ncores=15)
	lambda_search[[i]]<-res
	saveRDS(lambda_search,'incMDD/lambda_search.rds')
}



lambda_range<-lapply(lambda_search,function(x) x$lambda)


model_enet_CV_bootstrap_nostd<-list()
for( i in 1:9){
	print(' ')
	print(paste0('ALPHA=' ,i/10))
	print(' ')
	res<-elastic_net_CV_bootstrap(x=meth_bm, y=y, kfolds=10, train_valid=(1:length(y)), nlambda=100, alpha=(i/10), ncores=15)
	model_enet_CV_bootstrap_nostd[[i]]<-res
	saveRDS(model_enet_CV_bootstrap_nostd,'incMDD/model_enet_CV_bootstrap_nostd.rds')
}


model_enet_CV_bootstrap<-list()
for( i in 1:9){
	print(' ')
	print(paste0('ALPHA=' ,i/10))
	print(' ')
	res<-elastic_net_CV_bootstrap(x=meth_std_bm, y=y, kfolds=10, train_valid=(1:length(y)), nlambda=100, alpha=(i/10), ncores=15)
	model_enet_CV_bootstrap[[i]]<-res
	saveRDS(model_enet_CV_bootstrap,'incMDD/model_enet_CV_bootstrap.rds')
}

model_enet_CV_nobootstrap<-list()
for( i in 1:9){
	print(' ')
	print(paste0('ALPHA=' ,i/10))
	print(' ')
	res<-elastic_net_CV_nobootstrap(x=meth_std_bm, y=y, kfolds=10, train_valid=(1:length(y)), nlambda=100, alpha=(i/10), ncores=15)
	model_enet_CV_nobootstrap[[i]]<-res
	saveRDS(model_enet_CV_nobootstrap,'incMDD/model_enet_CV_nobootstrap.rds')
}

sel_1<-var_selection_best(res,method='f1_auc',kfolds=10)
sel_2<-var_selection_best(res,method='auc_sensitivity',kfolds=10) 


sel_3<-var_selection_union(res,method='f1_auc',kfolds=10)
sel_4<-var_selection_union(res,method='auc_sensitivity',kfolds=10) 







results_bart<-array(0,dim=c(4,5,10))
for(titer in 1:10){
	test<-sample(1:length(y),round(length(y)/10))
	train_valid<-(1:length(y))[-test]
	model_enet_CV_bootstrap<-list()
	for( i in 1:9){
		print(' ')
		print(paste0('ALPHA=' ,i/10))
		print(' ')
		res<-elastic_net_CV_bootstrap(x=meth_std_bm, y=y, kfolds=10, train_valid=train_valid, nlambda=100, alpha=(i/10), ncores=15)
		model_enet_CV_bootstrap[[i]]<-res
		saveRDS(model_enet_CV_bootstrap,'incMDD/titer/model_enet_CV_bootstrap.rds')
	}
	sel_1<-var_selection_best(model_enet_CV_bootstrap,method='f1_auc',kfolds=10,CVthreshold=1)
	sel_2<-var_selection_best(model_enet_CV_bootstrap,method='auc_sensitivity',kfolds=10,CVthreshold=1) 
	sel_3<-var_selection_union(model_enet_CV_bootstrap,method='f1_auc',kfolds=10,CVthreshold=5)
	sel_4<-var_selection_union(model_enet_CV_bootstrap,method='auc_sensitivity',kfolds=10,CVthreshold=1) 
	sel<-list(sel_1,sel_2,sel_3,sel_4)

	for(l in 1:length(sel) ){
		x<-meth_std[,sel[[l]]]
		bart_hyperpar_sel<-bart_CV(x,y,train_valid)
		bart_fit<-bart_test(x,y,valid_train,test,hyperpar=bart_hyperpar_sel$hyperpar_best)
			
		results_bart[l,,titer]<-bart_fit
		saveRDS(results_bart,'incMDD/results_bart.rds')
	}
}

	
