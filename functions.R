require(PRROC)

stratified_crossvalidation<-function(y,train_valid,breaks){
	cases<-which(y[train_valid]==1)
	controls<-which(!(y[train_valid]==1))
	folds<-numeric(length(train_valid))
	names(folds)<-train_valid
	folds_cases<-cut(sample(1:length(cases),length(cases)),breaks=breaks,labels=FALSE)
	folds_controls<-cut(sample(1:length(controls), length(controls)), breaks=breaks, labels=FALSE)
	folds[cases]<-folds_cases
	folds[controls]<-folds_controls
	return(folds)
}


train_bootstrap<-function(train,y){
	train_id_cases<-train[(y[train]==1)]
    add<-(length(train)-2*length(train_id_cases))
    train_add<-sample(train_id_cases,add,replace = TRUE)
    train_bootstrapped<-c(train,train_add)
	return(train_bootstrapped)
}

elastic_net_CV_bootstrap<-function(x,y,kfolds=10,nlambda,alpha, train_valid,ncores){
  folds <- stratified_crossvalidation(y,train_valid,kfolds)
  sensitivity<-matrix(0,nlambda,kfolds)
  specificity<-matrix(0,nlambda,kfolds)
  precision<-matrix(0,nlambda,kfolds)
  f1<-matrix(0,nlambda,kfolds)
  roc_auc<-matrix(0,nlambda,kfolds)
  pr_auc<-matrix(0,nlambda,kfolds)
  #nfeatures<-matrix(0,nlambda,kfolds)  
  #features<-list()
  #allfeatures<-list()
  train_valid_bootstrapped<-train_bootstrap(train_valid,y)
  get_lambdas<- biglasso(x, y, row.idx=train_valid_bootstrapped, penalty='enet',penalty.factor=penalty.factor, nlambda=nlambda,alpha=alpha, family='binomial',ncores = ncores)
  lambda_max<-get_lambdas$lambda[1]*1.1
  lambda<-c(seq(1.5*lambda_max,lambda_max,length.out=11)[-11], seq(lambda_max, get_lambdas$lambda[70], length.out=90))
  main_model<-biglasso(x, y, row.idx=train_valid_bootstrapped, penalty.factor=penalty.factor,penalty='enet', lambda=lambda,alpha=alpha, family='binomial',ncores = ncores)
  for(k in 1:kfolds){
    print(paste0(k,'th fold'))
    valid<-train_valid[which(folds==k)]
    train<-train_valid[which(folds!=k)]
    train_bootstrapped<-train_bootstrap(train,y)

    print(paste0('Proportion of cases in the train set is ', sum(y[train_bootstrapped]==1)/length(train_bootstrapped)))
    model <- biglasso(x, y,row.idx=train_bootstrapped,penalty.factor=penalty.factor,penalty='enet',lambda=lambda,alpha=alpha, family='binomial',ncores = ncores)
    print('TRAINED')
    pred<-as.matrix(predict(model,x,row.idx=valid,type='response'))
    roc_auc[,k]<-as.numeric(lapply(1:nlambda,function(z) measure_fit(pred[,z],y[valid])$roc_auc ))
    pr_auc[,k]<-as.numeric(lapply(1:nlambda,function(z) measure_fit(pred[,z],y[valid])$pr_auc ))
    precision[,k]<-as.numeric(lapply(1:nlambda,function(z) measure_fit(pred[,z],y[valid])$precision ))
	precision[is.nan(precision[,k]),k]<-0
    sensitivity[,k]<-as.numeric(lapply(1:nlambda,function(z) measure_fit(pred[,z],y[valid])$sensitivity ))
    specificity[,k]<-as.numeric(lapply(1:nlambda,function(z) measure_fit(pred[,z],y[valid])$specificity ))
    f1[,k]<-2*(precision[,k]*sensitivity[,k]/(precision[,k]+sensitivity[,k]))
	f1[is.nan(f1[,k]),k]<-0
	}
  nfeatures<-dim(x)[2]-main_model$rejections
  features<-rowSums(as.matrix(abs(main_model$beta)>1e-7))
  allfeatures<-(as.matrix(abs(main_model$beta)>1e-7))
  lambdas<-lambda
  return(list(lambda=lambda,nfeatures=nfeatures
              ,features=features,sensitivity=sensitivity,specificity=specificity,pr_auc=pr_auc
				,roc_auc=roc_auc,precision=precision,f1=f1,allfeatures=allfeatures))
}



elastic_net_CV_nobootstrap<-function(x,y,kfolds=10,nlambda,alpha, train_valid,ncores){
  folds <- stratified_crossvalidation(y,train_valid,kfolds)
  accuracy<-matrix(0,nlambda,kfolds)
  sensitivity<-matrix(0,nlambda,kfolds)
  specificity<-matrix(0,nlambda,kfolds)
  precision<-matrix(0,nlambda,kfolds)
  f1<-matrix(0,nlambda,kfolds)
  auc<-matrix(0,nlambda,kfolds)
  nfeatures<-matrix(0,nlambda,kfolds)  
  lambdas<-matrix(0,nlambda,kfolds)
  features<-list()
  allfeatures<-list()
  get_lambdas<- biglasso(x, y, row.idx=train_valid, penalty='enet', nlambda=nlambda,alpha=alpha, family='binomial',ncores = ncores)
  lambda_max<-get_lambdas$lambda[1]*1.1
  lambda<-c(seq(1.5*lambda_max,lambda_max,length.out=11)[-11], seq(lambda_max, get_lambdas$lambda[50], length.out=90))
  for(k in 1:kfolds){
    print(paste0(k,'th fold'))
    valid<-train_valid[which(folds==k)]
    train<-train_valid[which(folds!=k)]


    print(paste0('Bootstraped increased proportion of positive subjects in the train set is ', sum(y[train_bootstrapped]==1)/length(train_bootstrapped)))
    model <- biglasso(x, y,row.idx=train,penalty='enet',lambda=lambda,alpha=alpha, family='binomial',ncores = ncores)
    print('TRAINED')
    pred<-as.matrix(predict(model,x,row.idx=valid,type='response'))
    valid_matrix<-matrix(rep(as.numeric(as.character(y[valid])),nlambda),ncol=nlambda)
    auc[,k]<-as.numeric(lapply(1:nlambda,function(z) auc(valid_matrix[,z],pred[,z],quiet=TRUE) ))
    accuracy[,k]<-colMeans(as.matrix(valid_matrix==(pred>0.5)))
    lambdas[,k]<-model$lambda
    nfeatures[,k]<-dim(x)[2]-model$rejections
    features[[k]]<-rowSums(as.matrix(abs(model$beta)>1e-6))
    allfeatures[[k]]<-(as.matrix(abs(model$beta)>1e-6))
    precision[,k]<-colSums(as.matrix((valid_matrix==1) & (pred>0.5)))/colSums(as.matrix(pred>0.5))

	precision[is.nan(precision[,k]),k]<-0
    sensitivity[,k]<-colSums(as.matrix((valid_matrix==1) & (pred>0.5)))/colSums(valid_matrix)
    specificity[,k]<-colSums(as.matrix((valid_matrix==0) & (pred<0.5)))/colSums(!valid_matrix)
    f1[,k]<-2*(precision[,k]*sensitivity[,k]/(precision[,k]+sensitivity[,k]))
	f1[is.nan(f1[,k]),k]<-0

	}
  return(list(accuracy=accuracy,lambda=lambda,nfeatures=nfeatures
              ,features=features,sensitivity=sensitivity,specificity=specificity
				,auc=auc,precision=precision,f1=f1,allfeatures=allfeatures))
}




elastic_net<-function(x,y,lambda,alpha,train_valid,test,ncores){
  model <- biglasso(x, y,row.idx=train_valid, penalty='enet', lambda=lambda,alpha=alpha, family='binomial', ncores = ncores)
  pred<-as.matrix(predict(model,x,row.idx=test,type='response'))
  m<-measure_fit(pred,y[test])
  return(m)
}


#red<-list()
#for(i in 2:8){
#	red[[i]]<-Reduce(union,readRDS(paste0('incMDD/titer/sel_titer_',i,'.rds')))
#}
#big_union<-Reduce(union,red)
#big_union<-big_union[substr(big_union, 1, 2)=='cg']
#saveRDS(big_union,'incMDD/CpG_candidates.rds')
#x_varbvs<-meth[annotation$Sample_Sentrix_ID,big_union]
#saveRDS(x_varbvs,'incMDD/x_varbvs.rds')

measure_fit<-function(pred,y_test){
	y_test<-as.numeric(as.character(y_test))
	fg <- pred[y_test == 1]
	bg <- pred[y_test == 0]
	if(all(unique(fg) %in% unique(bg))){
		return(list(precision=0,sensitivity=0,specificity=0,f1=0,roc_auc=0,pr_auc=0))
	}
	roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

	pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
	m<-roc$curve[which.min(roc$curve[,1]^2+(1-roc$curve[,2])^2),]
	thr<-m[3]
	precision<-sum(as.numeric((y_test==1) & (pred>thr)))/sum(as.numeric(pred>thr))
	precision[is.nan(precision)]<-0
    sensitivity<-m[2]
    specificity<-1-m[1]
    f1<-2*(precision*sensitivity/(precision+sensitivity))
	f1[is.nan(f1)]<-0
	measures<-list(precision=precision,sensitivity=sensitivity,specificity=specificity,f1=f1,roc_auc=roc$auc,pr_auc=pr$auc.integral)
	return(measures)
}







var_selection_best<-function(res,method,titer){
	if(method=='roc_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$roc_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$roc_auc)))
	}
	if(method=='roc_auc_sensitivity'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$roc_auc*x$sensitivity)))
		max<-sapply(res,function(x) max(rowMeans(x$roc_auc*x$sensitivity)))
	}
	if(method=='f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$f1)))
	}
	if(method=='pr_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc)))
		saveRDS(list(best_lambda=best_lambda,max=max),paste0('prevMDD/titer/elnet_best_titer',titer,'.rds'))
	}
	if(method=='pr_auc_roc_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc*x$roc_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc*x$roc_auc)))
	}
	if(method=='pr_auc_f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc*x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc*x$f1)))
	}
	print(titer)
	best_alpha<-which.max(max)
	dim1<-dim(res[[best_alpha]]$allfeatures)[1]
	dim2<-dim(res[[best_alpha]]$allfeatures)[2]
	
	cpg<-((res[[best_alpha]])$allfeatures)[,best_lambda[best_alpha]]

	cpg<-cpg[-1]
	selected<-names(cpg)[which(cpg>(0))]
	return(selected)

}






var_selection_union<-function(res,method){
	if(method=='roc_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$roc_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$roc_auc)))
	}
	if(method=='roc_auc_sensitivity'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$roc_auc*x$sensitivity)))
		max<-sapply(res,function(x) max(rowMeans(x$roc_auc*x$sensitivity)))
	}
	if(method=='f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$f1)))
	}
	if(method=='pr_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc)))
	}
	if(method=='specificity'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$specificity)))
		max<-sapply(res,function(x) max(rowMeans(x$specificity)))
	}
	if(method=='pr_auc_roc_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc*x$roc_auc)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc*x$roc_auc)))
	}
	if(method=='pr_auc_f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$pr_auc*x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$pr_auc*x$f1)))
	}
	dim1<-dim(res[[1]]$allfeatures)[1]
	dim2<-dim(res[[1]]$allfeatures)[2]

	
	selected_union<-numeric(0)
	for(k in 1:length(res)){
		cpg<-numeric(length(dim1))	
		
		cpg<-((res[[k]])$allfeatures)[,best_lambda[k]]
		

		cpg<-cpg[-1]
		selected<-names(cpg)[which(cpg>(0))]
		selected_union<-union(selected_union,selected)	
	}
	return(selected_union)
}

#library(varbvs)
#fit.varbvs<-varbvs(X=x_varbvs,y=y,Z=annotation[,c('sex','age','smoke12','bmi','NK')],family='binomial')




find_correlated<-function(cpg,cpg_list,meth){
	corm<-cor(cbind(meth[,cpg],meth[,cpg_list]))
	return(corm[1,])
}



