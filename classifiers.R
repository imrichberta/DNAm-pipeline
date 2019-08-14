require(class)
require(PRROC)
require(BART)
library(doParallel)
library(randomForest)
library(e1071)
registerDoParallel(cores=20)
library(varbvs)

hyperpar_list<-function(){
	hyperpar_knn<-1:50

	ntrees<-c(20,50,100,200)
	powers<-c(0.5,2,4)
	ks<-c(4,6,8)
	hyperpar<-matrix(0,36,3)
	colnames(hyperpar)<-c('ntree','power','k')
	hyperpar[,1]<-rep(ntrees,each=9)
	hyperpar[,2]<-rep(powers,each=3)
	hyperpar[,3]<-rep(ks,12)
	hyperpar_bart<-hyperpar

	hyperpar<-matrix(0,180,3)
	ntree<-c(200,500,750,1000)
	mtry<-(seq((1/10),(1/2),length.out=5))
	nodesize<-c(1:9)
	colnames(hyperpar)<-c('ntree','mtry','nodesize')
	hyperpar[,1]<-rep(ntree,each=45)
	hyperpar[,2]<-rep(mtry,each=9)
	hyperpar[,3]<-rep(nodesize,20)
	hyperpar_rf<-hyperpar
	
	kernel<-c('linear','radial','sigmoid','polynomial')
	degree<-3:5
	gamma<-seq(1/(3),3,length.out=5)
	cost<-10^(-2:2)
	hyperpar<-matrix(NaN,length(cost)+2*length(cost)*length(gamma)+length(cost)*length(gamma)*length(degree),4)
	hyperpar<-as.data.frame(hyperpar)
	colnames(hyperpar)<-c('kernel','cost','gamma','degree')
	hyperpar[1:length(cost),1]<-'linear'
	hyperpar[(length(cost)+1):(length(cost)+length(cost)*length(gamma)),1]<-'radial'
	hyperpar[length(cost)*length(gamma)+((length(cost)+1):(length(cost)+length(cost)*length(gamma))),1]<-'sigmoid'
	hyperpar[56:130,1]<-'polynomial'
	hyperpar[,2]<-cost
	hyperpar[hyperpar$kernel!='linear',3]<-rep(gamma,each=length(cost))
	hyperpar[hyperpar$kernel=='polynomial',4]<-rep(degree,each=(length(cost)*length(gamma)))
	hyperpar_svm<-hyperpar

	return(list(hyperpar_bart=hyperpar_bart, hyperpar_rf=hyperpar_rf, hyperpar_knn=hyperpar_knn, hyperpar_svm=hyperpar_svm))
}

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
			auc_matrix[kn,k]<-as.numeric(auc(y[valid],pred,quiet=TRUE))
		}
	}
	best<-which.max(rowMeans(auc_matrix))
	
	return(list(auc=auc_matrix,best=best))
}

knn_test<-function(x,y,train_valid,test,kn){
	res <- knn(train=x[train_valid,], test=x[test,], cl=y[train_valid], k=kn,prob=TRUE)
	prob<-as.numeric(attributes(res)$prob)
	png(paste0('prevMDD/plots/AUC_knn',as.character(Sys.Date()),'.png')) 
	plot.roc(y[test],prob,
		     main = "Confidence intervals", percent=TRUE,ci = TRUE,  print.auc = TRUE)
	dev.off()

	return(list(res=measure_fit(prob,y[test]),pred=prob))
}


knn_hyperpar_test<-function(x,y,train_valid,test){
	roc_aucs<-numeric(50)
	pr_aucs<-numeric(50)
	for(kn in 1:50){
		res <- knn(train=x[train_valid,], test=x[test,], cl=y[train_valid], k=kn,prob=TRUE)
		pred<-as.numeric((attributes(res)$prob))                                     
		roc_aucs[kn]<-(measure_fit(pred,y[test]))$roc_auc
		pr_aucs[kn]<-(measure_fit(pred,y[test]))$pr_auc		
	}
	best<-which.max(roc_aucs)
	return(list(roc_aucs=roc_aucs,pr_aucs=pr_aucs,hyperpar=best))
}

##############

##############	BAYESIAN ADDITIVE REGRESSION TREES

##############






bart_hyperpar<-function(h,x_train,y_train,x_valid,y_valid,hyperpar){
	model <- lbart(x_train, y_train, ndpost=1000, power=hyperpar[h,2], ntree=hyperpar[h,1], k=hyperpar[h,3], nskip=1000)
	pred<-predict(model,x_valid)$prob.test.mean
	roc_auc<-(measure_fit(pred,y_valid))$roc_auc
	pr_auc<-(measure_fit(pred,y_valid))$pr_auc
	return(c(roc_auc,pr_auc))
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
		saveRDS(BART_logit_fit,'prevMDD/titer/BART_logit_fit.rds')
	}
	best<-which.max(rowMeans(BART_logit_fit))
	return(list(auc=BART_logit_fit,hyperpar_best=hyperpar[best,]))
}


bart_test<-function(x,y,train_valid,test,hyperpar){
	model <-lbart(x[train_valid,], as.numeric(as.character(y[train_valid])), ndpost=2000, power=hyperpar[2], ntree=hyperpar[1], k=hyperpar[3], nskip=2000)
	pred<-(predict(model,x[test,]))$prob.test.mean
	res<-measure_fit(pred,y_test=y[test])
	png(paste0('prevMDD/plots/AUC_bart',as.character(Sys.Date()),'.png')) 
	plot.roc(as.numeric(as.character(y[test])),as.numeric( pred),
		     main = "Confidence intervals", percent=TRUE,ci = TRUE,print.auc = TRUE)
	dev.off()
	return(list(res=res,pred=pred))
}


bart_hyperpar_test<-function(x,y,train_valid,test,l){
	ntrees<-c(20,50,100,200)
	powers<-c(0.5,2,4)
	ks<-c(4,6,8)
	hyperpar<-matrix(0,36,3)
	colnames(hyperpar)<-c('ntree','power','k')
	hyperpar[,1]<-rep(ntrees,each=9)
	hyperpar[,2]<-rep(powers,each=3)
	hyperpar[,3]<-rep(ks,12)
	BART_logit_fit<-numeric(36)
	folds<-stratified_crossvalidation(y,train_valid,10)
	
	auc_list<-foreach(i=1:36) %dopar% bart_hyperpar(i,x[train_valid,], as.numeric(as.character(y[train_valid])), x[test,], as.numeric(as.character(y[test])), hyperpar)
	roc_auc<-matrix(unlist(auc_list),nrow=2)[1,]
	pr_auc<-matrix(unlist(auc_list),nrow=2)[2,]

	#saveRDS(pr_auc,paste0('prevMDD/test_results/BART_logit_test_sel_',l,'.rds'))
	
	best<-which.max((pr_auc))
	return(list(roc_auc=roc_auc,pr_auc=pr_auc,hyperpar_best=hyperpar[best,]))
}



##############

##############	RANDOM FOREST

##############




rf_hyperpar<-function(h,x_train,y_train,x_valid,y_valid,hyperpar){
	print(paste0('Training Random Forest, ',c('ntree','mtry','nodesize'),':',hyperpar[h,])) 
	model <- randomForest(x=x_train, y=y_train,xtest=x_valid,ntree=hyperpar[h,1],mtry=hyperpar[h,2],nodesize=hyperpar[h,3])
	pred<-model$test$votes[,2]

	roc_auc<-(measure_fit(pred,y_valid))$roc_auc
	pr_auc<-(measure_fit(pred,y_valid))$pr_auc
	return(c(roc_auc,pr_auc))
}

rf_CV<-function(x,y,train_valid){
	ntree<-c(200,500,750,1000)
	mtry<-round(seq((dim(x)[2]/10),(dim(x)[2]/2),length.out=5))
	nodesize<-c(1:9)
	hyperpar<-matrix(0,180,3)
	colnames(hyperpar)<-c('ntree','mtry','nodesize')
	hyperpar[,1]<-rep(ntree,each=45)
	hyperpar[,2]<-rep(mtry,each=9)
	hyperpar[,3]<-rep(nodesize,20)
	rf_fit<-matrix(0,180,10)
	folds<-stratified_crossvalidation(y,train_valid,10)
	for(k in 1:10){
 		valid<-train_valid[which(folds==k)]
 		train<-train_valid[which(folds!=k)]
		auc_list<-foreach(i=1:180) %dopar% rf_hyperpar(i,x[train,], y[train], x[valid,], as.numeric(as.character(y[valid])), hyperpar)

		rf_fit[,k]<-unlist(auc_list)
		saveRDS(rf_fit,'prevMDD/titer/rf_fit.rds')
	}
	best<-which.max(rowMeans(rf_fit))
	return(list(auc=rf_fit,hyperpar_best=hyperpar[best,]))
}



rf_test<-function(x,y,train_valid,test,hyperpar){
	model <- randomForest(x=x[train_valid,], y=y[train_valid], xtest=x[test,], ntree=hyperpar[1],mtry=hyperpar[2],nodesize=hyperpar[3])
	pred<-model$test$votes[,2]


	res<-measure_fit(pred,y_test=y[test])
	png(paste0('prevMDD/plots/AUC_rf',as.character(Sys.Date()),'.png')) 
	plot.roc(as.numeric(as.character(y[test])),as.numeric( pred),
		     main = "Confidence intervals", 
		     percent=TRUE,
		     ci = TRUE,                
		     print.auc = TRUE)
	dev.off()
	return(list(res=res,pred=pred))
}


rf_hyperpar_test<-function(x,y,train_valid,test,l){
	ntree<-c(200,500,750,1000)
	mtry<-round(seq((dim(x)[2]/10),(dim(x)[2]/2),length.out=5))
	nodesize<-c(1:9)
	hyperpar<-matrix(0,180,3)
	colnames(hyperpar)<-c('ntree','mtry','nodesize')
	hyperpar[,1]<-rep(ntree,each=45)
	hyperpar[,2]<-rep(mtry,each=9)
	hyperpar[,3]<-rep(nodesize,20)
	rf_fit<-numeric(180)

	auc_list<-foreach(i=1:180) %dopar% rf_hyperpar(i,x[train_valid,], y[train_valid], x[test,], as.numeric(as.character(y[test])), hyperpar)

	roc_auc<-matrix(unlist(auc_list),nrow=2)[1,]
	pr_auc<-matrix(unlist(auc_list),nrow=2)[2,]


	#saveRDS(pr_auc,paste0('prevMDD/test_results/rf_fit_test_sel_',l,'.rds'))
	
	best<-which.max((pr_auc))
	return(list(roc_auc=roc_auc,pr_auc=pr_auc,hyperpar_best=hyperpar[best,]))
}

##############

##############	SUPPORT VECTOR MACHINE

##############

svm_hyperpar<-function(h,x_train,y_train,x_valid,y_valid,hyperpar){
		print(paste0('Training Support Vector Machine, ',c('kernel','cost','gamma','degree'),':',hyperpar[h,])) 
	y_valid<-as.numeric(as.character(y_valid))
	model<-svm(x=x_train, y=y_train, type = 'nu-classification', kernel=hyperpar[h,1], degree=hyperpar[h,4], gamma=hyperpar[h,3], cost=hyperpar[h,2],probability=TRUE,nu=0.2)
	pred<-as.numeric(attributes(predict(model,x_valid,probability=TRUE))$probabilities[,2])
	roc_auc<-(measure_fit(pred,y_valid))$roc_auc
	pr_auc<-(measure_fit(pred,y_valid))$pr_auc
	return(c(roc_auc,pr_auc))

}

svm_CV<-function(x,y,train_valid){
	p<-dim(x)[2]
	kernel<-c('linear','radial','sigmoid','polynomial')
	degree<-3:5
	gamma<-seq(1/(3*p),3/p,length.out=5)
	cost<-10^(-2:2)
	hyperpar<-matrix(0,length(cost)+2*length(cost)*length(gamma)+length(cost)*length(gamma)*length(degree),4)
	hyperpar<-as.data.frame(hyperpar)
	colnames(hyperpar)<-c('kernel','cost','gamma','degree')
	hyperpar[1:length(cost),1]<-'linear'
	hyperpar[(length(cost)+1):(length(cost)+length(cost)*length(gamma)),1]<-'radial'
	hyperpar[length(cost)*length(gamma)+((length(cost)+1):(length(cost)+length(cost)*length(gamma))),1]<-'sigmoid'
	hyperpar[hyperpar$kernel==0,1]<-'polynomial'
	hyperpar[,2]<-cost
	hyperpar[hyperpar$kernel!='linear',3]<-rep(gamma,each=length(cost))
	hyperpar[hyperpar$kernel=='polynomial',4]<-rep(degree,each=(length(cost)*length(gamma)))

	svm_fit<-matrix(0,130,10)
	folds<-stratified_crossvalidation(y,train_valid,10)
	for(k in 1:10){
 		valid<-train_valid[which(folds==k)]
 		train<-train_valid[which(folds!=k)]
		auc_list<-foreach(i=1:130) %dopar% svm_hyperpar(h,x[train,], y[train], x[valid,], y[valid], hyperpar)

		svm_fit[,k]<-unlist(auc_list)
		saveRDS(svm_fit,'prevMDD/titer/svm_fit.rds')
	}
	best<-which.max(rowMeans(svm_fit))
	return(list(auc=svm_fit,hyperpar_best=hyperpar[best,]))
}



svm_test<-function(x,y,train_valid,test,hyperpar){
	model<-svm(x=x[train_valid,], y=y[train_valid], type = 'nu-classification', kernel=hyperpar[1], degree=hyperpar[4], gamma=hyperpar[3], cost=hyperpar[2],nu = 0.2,probability=TRUE)
	pred<-as.numeric(attributes(predict(model,x[test,],probability=TRUE))$probabilities[,2])

	res<-measure_fit(pred,y_test=y[test])
	png(paste0('prevMDD/plots/AUC_svm',as.character(Sys.Date()),'.png')) 
	plot.roc(as.numeric(as.character(y[test])),as.numeric( pred),
		     main = "Confidence intervals", 
		     percent=TRUE,
		     ci = TRUE,                
		     print.auc = TRUE)
	dev.off()
	return(list(res=res,pred=pred))
}

svm_hyperpar_test<-function(x,y,train_valid,test,l){
	p<-dim(x)[2]
	kernel<-c('linear','radial','sigmoid','polynomial')
	degree<-3:5
	gamma<-seq(1/(3*p),3/p,length.out=5)
	cost<-10^(-2:2)
	hyperpar<-matrix(0,length(cost)+2*length(cost)*length(gamma)+length(cost)*length(gamma)*length(degree),4)
	hyperpar<-as.data.frame(hyperpar)
	colnames(hyperpar)<-c('kernel','cost','gamma','degree')
	hyperpar[1:length(cost),1]<-'linear'
	hyperpar[(length(cost)+1):(length(cost)+length(cost)*length(gamma)),1]<-'radial'
	hyperpar[length(cost)*length(gamma)+((length(cost)+1):(length(cost)+length(cost)*length(gamma))),1]<-'sigmoid'
	hyperpar[hyperpar$kernel==0,1]<-'polynomial'
	hyperpar[,2]<-cost
	hyperpar[hyperpar$kernel!='linear',3]<-rep(gamma,each=length(cost))
	hyperpar[hyperpar$kernel=='polynomial',4]<-rep(degree,each=(length(cost)*length(gamma)))
	svm_fit<-numeric(130)

	auc_list<-foreach(i=1:130) %dopar% svm_hyperpar(i,x[train_valid,], y[train_valid], x[test,], y[test], hyperpar)

	roc_auc<-matrix(unlist(auc_list),nrow=2)[1,]
	pr_auc<-matrix(unlist(auc_list),nrow=2)[2,]


	#saveRDS(pr_auc,paste0('prevMDD/test_results/vsm_fit_test_sel_',l,'.rds'))
	best<-which.max((pr_auc))
	return(list(pr_auc=pr_auc,roc_auc=roc_auc,hyperpar_best=hyperpar[best,]))
}




##############

##############	COVARIATE LOGISTIC REGRESSION

##############


logreg_cov<-function(annotation,train_valid,test){
	glm.fit <- glm(MDD ~ age + sex + smoke1 + smoke2 + smoke3,  data = annotation, family = binomial, subset = train_valid)
	glm.probs <- predict(glm.fit, newdata = annotation[test,c('age','sex','smoke1','smoke2','smoke3')], type = "response")
	return(measure_fit(glm.probs,y[test]))
}



##############

##############	VARIABLE IMPORTANCE BART

##############



bart_permute<-function(x_train,y_train,x_test,y_test,hyperpar){
	y_train<-as.numeric(as.character(y_train))
	model <- lbart(x_train, y_train, ndpost=2000, power=hyperpar[2], ntree=hyperpar[1], k=hyperpar[3], nskip=1000,printevery=1000L)

	y_test<-as.numeric(as.character(y_test))


	pred<-predict(model,x_test)$prob.test.mean
	roc_auc_fit<-(measure_fit(pred,y_test))$roc_auc
	pr_auc_fit<-(measure_fit(pred,y_test))$pr_auc

	permute_roc_auc<-numeric(ncol(x_train))
	permute_pr_auc<-numeric(ncol(x_train))
	for(var in 1:ncol(x_train)){
		x<-x_test
		x[,var]<-sample(x_test[,var],replace=FALSE)
		pred_permute<-predict(model,x)$prob.test.mean
		permute_roc_auc[var]<-(measure_fit(pred_permute,y_test))$roc_auc
		permute_pr_auc[var]<-(measure_fit(pred_permute,y_test))$pr_auc

	}
	return(list(roc_auc_fit=roc_auc_fit, pr_auc_fit=pr_auc_fit, permute_roc_auc=(permute_roc_auc-roc_auc_fit), permute_pr_auc=(permute_pr_auc-pr_auc_fit)))
}





