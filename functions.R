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
  train_valid_bootstrapped<-train_bootstrap(train_valid,y)
  get_lambdas<- biglasso(x, y, row.idx=train_valid_bootstrapped, penalty='enet', nlambda=nlambda,alpha=alpha, family='binomial',ncores = ncores)
  lambda_max<-get_lambdas$lambda[1]*1.1
  lambda<-c(seq(1.5*lambda_max,lambda_max,length.out=11)[-11], seq(lambda_max, get_lambdas$lambda[50], length.out=90))
  for(k in 1:kfolds){
    print(paste0(k,'th fold'))
    valid<-train_valid[which(folds==k)]
    train<-train_valid[which(folds!=k)]
    train_bootstrapped<-train_bootstrap(train,y)

    print(paste0('Bootstraped increased proportion of positive subjects in the train set is ', sum(y[train_bootstrapped]==1)/length(train_bootstrapped)))
    model <- biglasso(x, y,row.idx=train_bootstrapped,penalty='enet',lambda=lambda,alpha=alpha, family='binomial',ncores = ncores)
    print('TRAINED')
    pred<-as.matrix(predict(model,x,row.idx=valid,type='response'))
    valid_matrix<-matrix(rep(as.numeric(as.character(y[valid])),nlambda),ncol=nlambda)
    auc[,k]<-as.numeric(lapply(1:nlambda,function(z) auc(valid_matrix[,z],pred[,z]) ))
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
    auc[,k]<-as.numeric(lapply(1:nlambda,function(z) auc(valid_matrix[,z],pred[,z]) ))
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




elastic_net<-function(x,y,nlambda,alpha, train_valid,ncores){
  model <- biglasso(x, y,row.idx=train_valid, penalty='enet', nlambda=nlambda,alpha=alpha, family='binomial', ncores = ncores)
  print('TRAINED')
  lambda<-model$lambda
  nfeatures<-dim(x)[2]-model$rejections
  features<-rowSums(abs(model$beta)>1e-6)
  return(list(lambda=lambda,nfeatures=nfeatures
              ,features=features))
}







var_selection_best<-function(res,method,kfolds,CVthreshold){
	if(method=='auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$auc)))
		max<-sapply(res,function(x) max(rowMeans(x$auc)))
	}
	if(method=='auc_sensitivity'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$auc*x$sensitivity)))
		max<-sapply(res,function(x) max(rowMeans(x$auc*x$sensitivity)))
	}
	if(method=='f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$f1)))
	}
	if(method=='f1_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1*x$auc)))
		max<-sapply(res,function(x) max(rowMeans(x$f1*x$auc)))
	}

	best_alpha<-which.max(max)
	dim1<-dim(res[[best_alpha]]$allfeatures[[1]])[1]
	dim2<-dim(res[[best_alpha]]$allfeatures[[1]])[2]

	
	cpg<-numeric(length(dim1))	

	for(i in 1:kfolds){
		cpg<-cpg+((res[[best_alpha]])$allfeatures[[i]])[,best_lambda[best_alpha]]
	}

	cpg<-cpg[-1]
	selected<-names(cpg)[which(cpg>(CVthreshold))]
	return(selected)

}






var_selection_union<-function(res,method,kfolds,CVthreshold){
	if(method=='auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$auc)))
		max<-sapply(res,function(x) max(rowMeans(x$auc)))
	}
	if(method=='auc_sensitivity'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$auc*x$sensitivity)))
		max<-sapply(res,function(x) max(rowMeans(x$auc*x$sensitivity)))
	}
	if(method=='f1'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1)))
		max<-sapply(res,function(x) max(rowMeans(x$f1)))
	}
	if(method=='f1_auc'){
		best_lambda<-sapply(res,function(x) which.max(rowMeans(x$f1*x$auc)))
		max<-sapply(res,function(x) max(rowMeans(x$f1*x$auc)))
	}


	dim1<-dim(res[[1]]$allfeatures[[1]])[1]
	dim2<-dim(res[[1]]$allfeatures[[1]])[2]

	
	selected_union<-numeric(0)
	for(k in 1:9){
		cpg<-numeric(length(dim1))	
		for(i in 1:kfolds){
			cpg<-cpg+((res[[k]])$allfeatures[[i]])[,best_lambda[k]]
		}

		cpg<-cpg[-1]
		selected<-names(cpg)[which(cpg>(CVthreshold))]
		selected_union<-union(selected_union,selected)	
	}
	return(selected_union)
}



