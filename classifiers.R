require(class)
require(pROC)
require(BART)
library(doParallel)
registerDoParallel(cores=12)

knn_CV<-function(x,y,train_valid){
	auc_matrix<-matrix(0,50,10)	
	cases<-which(y[train_valid]==1)
	controls<-which(!(y[train_valid]==1))
	folds<-numeric(length(train_valid))
	names(folds)<-train_valid
	folds_cases<-cut(sample(1:length(cases),length(cases)),breaks=10,labels=FALSE)
	folds_controls<-cut(sample(1:length(controls), length(controls)), breaks=10, labels=FALSE)
	folds[cases]<-folds_cases
	folds[controls]<-folds_controls

	for(k in 1:10){
		valid<-train_valid[which(folds==k)]
    	train<-train_valid[which(folds!=k)]
		for(kn in 1:50){
			res <- knn(train=x[train,], test=x[valid,], cl=y[train], k=kn,prob=TRUE)
			
			pred<-as.numeric((attributes(res)$prob))                                     
			auc_matrix[kn,k]<-as.numeric(auc(y[valid],pred))
		}
	}
	best<-which.max(rowMeans(auc_matrix))
	
	return(list(auc=auc_matrix,best=best))
}

knn_test<-function(x,y,valid_train,test,kn){
	res <- knn(train=x[train_valid,], test=x[test,], cl=y[train_valid], k=kn)
	png(paste0('incMDD/plots/AUC_knn',as.character(Sys.Date()),'.png')) 
	plot.roc(y[test],as.numeric( knn),
		     main = "Confidence intervals", 
		     percent=TRUE,
		     ci = TRUE,                  # compute AUC (of AUC by default)
		     print.auc = TRUE)
	dev.off()



	return(as.numeric(attributes(res)$prob))
}



bart_hyperpar<-function(h,x_train,y_train,x_valid,y_valid,hyperpar){

	model <- lbart(x_train, y_train, ndpost=1000, power=hyperpar[h,2], ntree=hyperpar[h,1], k=hyperpar[h,3], nskip=1000)
	pred<-predict(model,x_valid)
	auc<-as.numeric(auc(as.numeric(as.character(y_valid)),pred$prob.test.mean))
	return(auc)
}



bart_CV<-function(x,y,train_valid){
	ntrees<-c(20,50,100,200)
	powers<-c(0.5,2,4)
	ks<-c(2,4,6)
	hyperpar<-matrix(0,36,3)
	colnames(hyperpar)<-c('ntree','power','k')
	hyperpar[,1]<-rep(ntrees,each=9)
	hyperpar[,2]<-rep(powers,each=3)
	hyperpar[,3]<-rep(ks,12)
	BART_logit_fit<-matrix(0,36,10)
	folds<-stratified_crossvalidation(y,train_valid,10)
	for(k in 1:10){
 		valid<-train_valid[which(folds==k)]
 		train<-train_valid[which(folds!=k)]

		auc_list<-foreach(i=1:36) %dopar% bart_hyperpar(i,x[train,], as.numeric(as.character(y[train])), x[valid,], as.numeric(as.character(y[valid])), hyperpar)

		BART_logit_fit[,k]<-unlist(auc_list)
		saveRDS(BART_logit_fit,'incMDD/titer/BART_logit_fit.rds')
		
	}
	best<-which.max(rowMeans(BART_logit_fit))
	return(list(auc=BART_logit_fit,hyperpar_best=hyperpar[best,]))
}


bart_test<-function(x,y,valid_train,test,hyperpar){
	model <-lbart(x[train_valid,], y[train_valid], ndpost=2000, power=hyperpar[2], ntree=hyperpar[1], k=hyperpar[3], nskip=2000)

	pred<-(predict(model,x[test,]))$prob.test.mean

	res<-measure_fit(pred,y_test=y[test])
	png(paste0('incMDD/plots/AUC_bart',as.character(Sys.Date()),'.png')) 
	plot.roc(as.numeric(as.character(y[test])),as.numeric( pred),
		     main = "Confidence intervals", 
		     percent=TRUE,
		     ci = TRUE,                
		     print.auc = TRUE)
	dev.off()
	return(res)
}







measure_fit<-function(pred,y_test){
	y_test<-as.numeric(as.character(y_test))
	auc<-as.numeric(auc(y_test,pred))
	precision<-sum(as.numeric((y_test==1) & (pred>0.5)))/sum(as.numeric(pred>0.5))

	precision[is.nan(precision)]<-0
    sensitivity<-sum(as.numeric((y_test==1) & (pred>0.5)))/sum(y_test)
    specificity<-sum(as.numeric((y_test==0) & (pred<0.5)))/sum(!y_test)
    f1<-2*(precision*sensitivity/(precision+sensitivity))
	f1[is.nan(f1)]<-0
	return(c(precision,sensitivity,specificity,f1,auc))
}



