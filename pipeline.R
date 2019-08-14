meth<-readRDS('meth.rds')
annotation<-readRDS('annotation.rds')

library(pROC)
library(biglasso)
library(BART)

meth<-as.data.frame(meth[which(annotation$unrelated==1),])

y<-annotation$MDD[which(annotation$unrelated==1)]
meth <- as.big.matrix(meth)

elastic_net<-function(x,y,kfolds=10,nlambda,alpha, train_valid,prob=0.5,shuffle=TRUE,ncores){
  folds <- cut(sample(1:length(train_valid),length(train_valid)),breaks=kfolds,labels=FALSE)
  accuracy<-matrix(0,nlambda,kfolds)
  sensitivity<-matrix(0,nlambda,kfolds)
  specificity<-matrix(0,nlambda,kfolds)
  precision<-matrix(0,nlambda,kfolds)
  f1<-matrix(0,nlambda,kfolds)
  auc<-matrix(0,nlambda,kfolds)
  nfeatures<-matrix(0,nlambda,kfolds)  
  lambda<-matrix(0,nlambda,kfolds)
  features<-list()
  for(k in unique(folds)){
    print(paste0(k,'th fold'))
    valid<-train_valid[which(folds==k)]
    train_id<-train_valid[which(folds!=k)]
    train_id_MDD<-train_id[(y[train_id]==1)]
    add<-(length(train_id)-2*length(train_id_MDD))
    train_id_add<-sample(train_id_MDD,add,replace = TRUE)
    train<-c(train_id,train_id_add)
    print(paste0('Bootstraped increased proportion of positive subjects in the train set is ', length(c(train_id_MDD,train_id_add))/length(train)))
    model <- biglasso(x, y,row.idx=train,penalty='enet',nlambda=nlambda,alpha=alpha, family='binomial',ncores = ncores)
    print('TRAINED')
    pred<-predict(model,x,row.idx=valid,type='response')
    valid_matrix<-matrix(rep(y[valid],nlambda),ncol=nlambda)
    auc[,k]<-as.numeric(lapply(1:nlambda,function(x) auc(valid_matrix[,x],pred[,x]) ))
    accuracy[,k]<-colMeans(valid_matrix==(pred>prob))
    lambda[,k]<-model$lambda
    nfeatures[,k]<-dim(x)[2]-model$rejections
    features[[k]]<-rowSums(abs(model$beta)>1e-6)
    precision[,k]<-colSums((valid_matrix==1) & (pred>prob))/colSums(pred>prob)
    sensitivity[,k]<-colSums((valid_matrix==1) & (pred>prob))/colSums(valid_matrix)
    specificity[,k]<-colSums((valid_matrix==0) & (pred<prob))/colSums(!valid_matrix)
    f1[,k]<-2*(precision[,k]*sensitivity[,k]/(precision[,k]+sensitivity[,k]))
	}
  return(list(accuracy=accuracy,lambda=lambda,nfeatures=nfeatures
              ,features=features,sensitivity=sensitivity,specificity=specificity
				,auc=auc,precision=precision,f1=f1))
}


annotation<-readRDS('annotationF.rds')
annotation$smoke[is.na(annotation$smoke)]<-4
annotation$smoke12<-1*(annotation$smoke < 3)
annotation<-annotation[annotation$unrelated==1,]
annotation[which(is.na(annotation),arr.ind=T)]<-mean(annotation[,14],na.rm=T)

res_test<-elastic_net(x=meth.bm,y=y,kfolds=10,nlambda=300, length.out=100),nlambda=100, alpha=0.4,shuffle=TRUE,ncores=15)


BART_logit_auc1<-matrix(0,20,10)
BART_logit_auc2<-matrix(0,20,10)
for(titer in 1:10){
  test<-sample(1:length(y),round(length(y)/10))
	train_valid<-(1:length(y))[-test]
	results_enet<-list()
	for( i in 1:4){
		print(' ')
		print(paste0('ALPHA=' ,i/10))
		print(' ')
		res<-elastic_net(x=meth,y=y,kfolds=10, train_valid=train_valid,nlambda=300
		,alpha=i/10,shuffle=TRUE,ncores=15)
		results_enet[[i]]<-res
	}
	idcpg<-colnames(meth)
	selected<-data.frame(matrix(0,length(idcpg),5))
	names(selected)[1]<-c('name')
	selected$name<-idcpg
	for(d in 1:4){
		for(i in 1:10){
			selected[,d+1]<-selected[,d+1]+(results_enet[[d]])$features[[i]][-1]
		}
	}
	cpg_selected<-selected[which(rowSums(selected[,2:5]>1000)>1),]
	names(cpg_selected)[1]<-'name'


  X2<-cbind(annotation[,c('age','sex','bmi','smoke12')],meth[,cpg_selected$name])
	X1<-meth[,cpg_selected$name]
	X2<-as.matrix(X2)
	X1<-as.matrix(X1)

	for(tree in 1:20){

		model1 <- mc.lbart(X1[train_valid,], y[train_valid],ndpost=1000,base=0.95, ntree=20*tree,nskip=2000,mc.cores=10)
		model2 <- mc.lbart(X2[train_valid,], y[train_valid],ndpost=1000,base=0.95, ntree=20*tree,nskip=2000,mc.cores=10)
		
		pred1<-predict(model1,X1[test,])
		auc1<-as.numeric(auc(y[test],pred1$prob.test.mean))

		pred2<-predict(model2,X2[test,])
		auc2<-as.numeric(auc(y[test],pred2$prob.test.mean))


		BART_logit_auc1[tree,titer]<-auc1
		BART_logit_auc2[tree,titer]<-auc2
		saveRDS(BART_logit_auc1,'results/bart_logit_reduced1.rds')
		saveRDS(BART_logit_auc2,'results/bart_logit_reduced2.rds')
	}
}


annotation$smoke<-smoke$ever_smoke[match(annotation$Sample_Name,smoke$Sample_Name)]





idcpg<-colnames(meth)
selected<-data.frame(matrix(0,length(idcpg),10))
names(selected)<-c('name','alpha010','alpha020','alpha030','alpha040','alpha050','alpha060','alpha070','alpha080','alpha090')
selected$name<-idcpg
for(d in 1:9){
	for(i in 1:10){
		selected[,d+1]<-selected[,d+1]+(res[[d]])$features[[i]][-1]
	}
}

selected <- selected[order(selected$alpha040,decreasing=T),] 

cpg_1500_6<-selected[which(rowSums(selected[,2:10]>1200)>6),]
cpg_1200_4<-selected[which(rowSums(selected[,2:10]>1200)>4),]
cpg_1000_2<-selected[which(rowSums(selected[,2:10]>1000)>2),]
cpg_500_1<-selected[which(rowSums(selected[,2:10]>500)>1),]
names(cpg_1500_6)[1]<-'name'
names(cpg_1200_4)[1]<-'name'
names(cpg_1000_2)[1]<-'name'
names(cpg_500_1)[1]<-'name'

saveRDS(meth[,cpg_500_1$name],'reduced/x_GS_cpg_500_1.rds')
saveRDS(meth[,cpg_1000_2$name],'reduced/x_GS_cpg_1000_2.rds')
saveRDS(meth[,cpg_1200_4$name],'reduced/x_GS_cpg_1200_4.rds')
saveRDS(meth[,cpg_1500_6$name],'reduced/x_GS_cpg_1500_6.rds')
saveRDS(annotation[which(annotation$unrelated==1),],'reduced/annotation_GS.rds')


####################################### Incident DEPRESSION ####################################


meth<-readRDS('meth.rds')
annotation<-readRDS('annotation.rds')

library(pROC)
library(biglasso)

meth<-as.data.frame(meth[which(annotation$unrelated==1),])

y<-annotation$MDD[which(annotation$unrelated==1)]
meth <- as.big.matrix(meth)

tmp<-readRDS('DNAm_incidence_cases_controls.rds')                                 
annotation$MDD_w3inc<-tmp$dep_status[(match(annotation$Sample_Sentrix_ID,tmp$ID))]
saveRDS(annotation,'annotation_MDD_w3i.rds')







