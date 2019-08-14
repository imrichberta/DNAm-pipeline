
library(gridExtra)
library(ggplot2)
library(pROC)
library(biglasso)
library(BART)



meth<-readRDS('meth.rds')


standardise<-function(data){
	return((data-mean(data))/sd(data))
}


annotation<-readRDS('incMDD/annotationF.rds')


#annotation$smoke[is.na(annotation$smoke)]<-4
annotation$smoke1<-0
annotation$smoke2<-0
annotation$smoke3<-0

annotation$smoke[is.na(annotation$smoke)]<-4

annotation$smoke1[(annotation$smoke)==1]<-1
annotation$smoke2[(annotation$smoke)==2]<-1
annotation$smoke3[(annotation$smoke)==3]<-1


idx<-which(!is.na(annotation$MDD) &  annotation$unrelated==1)

meth<-meth[idx,]
meth_std<-apply(meth,2,standardise)
meth<-NA



annotation<-annotation[idx,]

y<-annotation$MDD

y<-as.factor(y)

source('functions.R')
source('classifiers.R')

covariates<-c('sex','age','Bcell','CD4T','CD8T','Gran','NK','Mono','bmi','smoke1','smoke2','smoke3')

meth_std_cov<-as.data.frame(cbind(annotation[,covariates],meth_std))

meth_std_cov_bm<-as.big.matrix(meth_std_cov)


##### TAKE UNION OF SELVARS on DIFFERENT ALPHA LEVELS

### The bigLASSO package offers the user a cross-validation hyperparameter autotunning function. However, the metric choice is only limited to misclassification error (MSE and MAPE for regression). We prefer to evaluate the goodness of fit with area under reciever operator curve and F1 score - harmonic mean of precision and recall. Therefore, we construct the cross-validation step manually, measuring the AUC and F1 score.


### For given level of alpha, too low lambda values results in unstable convergence and costly computations. The bigLASSO package can precompute a feasible reagion of lambdas at given alpha level when it is given number of lambda (nlambda). At one hand, this is very useful for us as users because we don't have to calibrate the min/max lambda, which can be very sensitive with high dimensional data, but when we want to constructy the k-fold crossvalidation on our own the package would use non-consistent lambda values at each k. Therefore we use the inner package function to find the feasible range of lambdas at different alpha levels


penalty.factor<-numeric(dim(meth_std_cov)[2])
penalty.factor<-penalty.factor+1
penalty.factor[c(1,2,7,10,11,12)]<-0.01




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
		saveRDS(model_enet_CV_bootstrap,'prevMDD/titer/model_enet_CV_bootstrap.rds')
	}
	sel_1<-var_selection_best(model_enet_CV_bootstrap,method='f1',titer)
	sel_2<-var_selection_best(model_enet_CV_bootstrap,method='pr_auc',titer)
	sel_3<-var_selection_best(model_enet_CV_bootstrap,method='pr_auc_f1',titer) 
	sel_4<-var_selection_union(model_enet_CV_bootstrap,method='f1')
	sel_5<-var_selection_union(model_enet_CV_bootstrap,method='pr_auc')
	sel_6<-var_selection_union(model_enet_CV_bootstrap,method='pr_auc_f1') 
	sel_all<-Reduce(union,list(sel_1,sel_2,sel_3,sel_4,sel_5,sel_6))


	x_varbvs<-meth_std_cov[,c(covariates[-(1:2)],sel_all)]
	x_varbvs<-as.matrix(x_varbvs)
	fit.varbvs<-varbvs(X=x_varbvs[train_valid,],y=as.numeric(as.character(y[train_valid])),Z=annotation[train_valid,c('sex','age')],family='binomial')
	sel_varbvs<-names(fit.varbvs$pip[fit.varbvs$pip>0.95])
	sel_varbvs<-c('sex','age',sel_varbvs)
	sel<-list(sel_1,sel_2,sel_3,sel_4,sel_5,sel_6,sel_varbvs)
	
	saveRDS(sel,paste0('prevMDD/titer/sel_titer_',titer,'.rds') )
	
	results_logreg_cov_test[[titer]]<-logreg_cov(annotation,train_valid,test)
	for(l in 7:length(sel) ){
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



		saveRDS(results_bart_test,'prevMDD/test_results/results_bart_test.rds')
		saveRDS(results_rf_test,'prevMDD/test_results/results_rf_test.rds')
		saveRDS(results_knn_test,'prevMDD/test_results/results_knn_test.rds')
		saveRDS(results_svm_test,'prevMDD/test_results/results_svm_test.rds')
		saveRDS(results_logreg_cov_test,'prevMDD/test_results/results_logreg_cov_test.rds')


	}

}

	

#################################

################################# VARIABLE IMPORTANCE

#################################





set.seed(2020)


x_cand<-readRDS('prevMDD/sel_varbvs.rds')

annotation<-readRDS('incMDD/annotationF.rds')
idx<-which(annotation$unrelated==1)


meth<-readRDS('meth.rds')
standardise<-function(data){
	return((data-mean(data))/sd(data))
}
meth<-meth[idx,]
meth_std<-apply(meth,2,standardise)
meth<-NA


x_cov<-meth_std[,colnames(meth_std) %in% x_cand]


set.seed(2022)

annotation<-readRDS('incMDD/annotationF.rds')
idx<-which(annotation$unrelated==1)



x_cov<-readRDS('prevMDD/x_varbvs.rds')



annotation<-annotation[idx,]
annotation$smoke1<-0
annotation$smoke2<-0
annotation$smoke3<-0

annotation$smoke[is.na(annotation$smoke)]<-4
annotation$smoke1[(annotation$smoke)==1]<-1
annotation$smoke2[(annotation$smoke)==2]<-1
annotation$smoke3[(annotation$smoke)==3]<-1

y<-annotation$MDD
y<-as.factor(y)
covariates<-c('sex','age','Bcell','CD4T','CD8T','Gran','NK','Mono','bmi','smoke1','smoke2','smoke3')

x_cov<-cbind(as.matrix(annotation[,covariates]),x_cov)
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

	saveRDS(varbvs_selection,'prevMD/varbvs_selection.rds')

	for(iter in 1:10){
		valid<-which(folds==iter_vec[iter+1])
		x<-x_cov[,varbvs_selection[[iter]] ]
		train<-which((folds!=iter)&(folds!=(iter_vec[iter+1])))
		x<-cbind(annotation[,covariates],x)	
		res_valid<-foreach(h=1:36) %dopar% {

			print('----------------------')
			print(paste0('VALIDATION:',iter))
			print('----------------------')

			model<-lbart(x[train,], as.numeric(as.character(y[train])), ndpost=2000, power=hyperpar[h,2], ntree=hyperpar[h,1],
				k=hyperpar[h,3], nskip=1000,printevery=1000L)
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

		saveRDS(pr_change_array,'prevMDD/train_valid_test/results_22CV_varimp_pr.rds')
		saveRDS(roc_change_array,'prevMDD/train_valid_test/results_22CV_varimp_roc.rds')
		saveRDS(pr_auc_test,'prevMDD/train_valid_test/results_22CV_pr_auc_test.rds')
		saveRDS(roc_auc_test,'prevMDD/train_valid_test/results_22CV_roc_auc_test.rds')
		saveRDS(pr_auc_logreg,'prevMDD/train_valid_test/results_22CV_pr_auc_logreg.rds')
		saveRDS(roc_auc_logreg,'prevMDD/train_valid_test/results_22CV_roc_auc_logreg.rds')

	}






