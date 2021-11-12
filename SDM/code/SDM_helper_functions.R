# helper functions for SDM course 1-5 March 2021

# Select 07
# Select weakly correlated variables based on univariate importance based on Dormann et al. (2013). Univariate variable importance is based on AIC. Variable importance can also be pre-defined by hand.


select07 <- function(X, y, family="binomial",univar="glm2", threshold=0.7, method="spearman",sequence=NULL)
{
  require(mgcv)
  # selects variables based on removing correlations > 0.7, retaining those
  # variables more important with respect to y
  # Order of importance can be provided by the character vector 'sequence'
  
  # 1. step: cor-matrix
  # 2. step: importance vector
  # 3. step: identify correlated pairs
  # 4. step: in order of importance: remove collinear less important variable,
  #           recalculate correlation matrix a.s.f.
  
  var.imp <- function (variable, response, univar=univar,family="gaussian")
  {
    # calculates the univariate (=marginal) importance of a variable for a response
    if(univar=="glm1")
    {
      fm.glm <- glm(response ~ variable, family=family)
      summary(fm.glm)$aic
    } else
      if(univar=="glm2")
      {
        fm.glm <- glm(response ~ poly(variable,2), family=family)
        summary(fm.glm)$aic
      } else  
        if(univar=="gam")
        {
          fm.gam <- mgcv::gam(response ~ s(variable,k=4), family=family)
          AIC(fm.gam)
        } else return(F)
  }
  
  cm <- cor(X, method=method)
  
  if (is.null(sequence)) {
    a<-try(var.imp(X[,1],y,univar=univar))
    if (is.numeric(a)!=1) {stop("invalid univar method")}
    
    imp <- apply(X, 2, var.imp, response=y, family=family,univar=univar) #importance as AIC: the lower the better!
    sort.imp<-names(sort(imp)) 
  } else
  { sort.imp <- sequence }
  
  pairs <- which(abs(cm)>= threshold, arr.ind=T) # identifies correlated variable pairs
  index <- which(pairs[,1]==pairs[,2])           # removes entry on diagonal
  pairs <- pairs[-index,]                        # -"-
  
  exclude <- NULL
  for (i in 1:length(sort.imp))
  {
    if ((sort.imp[i] %in% row.names(pairs))&
        ((sort.imp[i] %in% exclude)==F)) {
      cv<-cm[setdiff(row.names(cm),exclude),sort.imp[i]]
      cv<-cv[setdiff(names(cv),sort.imp[1:i])]
      exclude<-c(exclude,names(which((abs(cv)>=threshold)))) }
  }
  
  pred_sel <- sort.imp[!(sort.imp %in% unique(exclude)),drop=F]
  return(list(AIC=sort(imp), cor_mat=cm, pred_sel=pred_sel))
}


#---------------------------

# Explained Deviance

expl_deviance <- function(obs, pred, family='binomial'){
  require(dismo)
  pred <- ifelse(pred<.00001,.00001,ifelse(pred>.9999,.9999,pred))
  
  null_pred <- rep(mean(obs), length(obs))
  
  1 - (dismo::calc.deviance(obs, pred, family=family) / 
         dismo::calc.deviance(obs, null_pred, family=family))
}


#---------------------

# True skill statistic

TSS = function(cmx){
  require(PresenceAbsence)
  PresenceAbsence::sensitivity(cmx, st.dev=F) + 
    PresenceAbsence::specificity(cmx, st.dev=F) - 1
}


#---------------------

# evaluation statistics
evalSDM <- function(observation, predictions, thresh.method='MaxSens+Spec'){
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  thresh <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1)
  cmx.opt <- PresenceAbsence::cmx(DATA= thresh.dat, threshold=thresh[thresh$Method==thresh.method,2])
  
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
             TSS = TSS(cmx.opt), 
             Kappa = PresenceAbsence::Kappa(cmx.opt, st.dev=F),
             Sens = PresenceAbsence::sensitivity(cmx.opt, st.dev=F),
             Spec = PresenceAbsence::specificity(cmx.opt, st.dev=F),
             PCC = PresenceAbsence::pcc(cmx.opt, st.dev=F),
             D2 = expl_deviance(observation, predictions),
             thresh = thresh[thresh$Method==thresh.method,2])
}

#-------------------------

# Make predictions

predictSDM <- function(model, newdata) {
  require(dismo)
  require(rpart)
  require(randomForest)
  require(maxnet)
  require(gam)
  
  switch(class(model)[1],
         Bioclim = predict(model, newdata),
         Domain = predict(model, newdata),
         glm = predict(model, newdata, type='response'),
         Gam = predict(model, newdata, type='response'),
         rpart = predict(model, newdata),
         randomForest = switch(model$type,
                               regression = predict(model, newdata, type='response'),
                               classification = predict(model, newdata, type='prob')[,2]),
         gbm = predict.gbm(model, newdata, 
                           n.trees=model$gbm.call$best.trees, type="response"),
         maxnet = predict(model, newdata, type="logistic"))
}

#------------------------

# cross validation

crossvalSDM <- function(model, kfold=5, traindat, colname_species, colname_pred,
                        env_r=NULL, colname_coord=NULL) {
  
  require(dismo)
  require(rpart)
  require(randomForest)
  require(maxnet)
  require(gam)
  
  if (length(kfold)==1) {
    # Make k-fold data partitions
    ks <- dismo::kfold(traindat, k = kfold, by = traindat[,colname_species])
  } else {
    ks <- kfold
    kfold <- length(unique(kfold))
  }
  
  cross_val_preds = numeric(length = nrow(traindat))
  
  for(i in seq_len(kfold)){
    cv_train <- traindat[ks!=i,]
    cv_test <- traindat[ks==i,]
    
    # Because we used the gbm.step() for BRTs, we need a small work-around:
    if (class(model)[1]=='gbm') {
      cv_train_gbm <- cv_train;
      names(cv_train_gbm)[names(cv_train_gbm)==colname_species] <- 
        model$response.name
    }
    
    # We update the model for the new training data
    modtmp <- switch(class(model)[1],
                     Bioclim = dismo::bioclim(env_r[[colname_pred]], cv_train[cv_train[, colname_species]==1, colname_coord]),
                     Domain = dismo::domain(env_r[[colname_pred]], cv_train[cv_train[, colname_species]==1, colname_coord]),
                     glm = update(model, data=cv_train),
                     Gam = update(model, data=cv_train),
                     rpart = rpart::update(model, data=cv_train),
                     randomForest = update(model, data=cv_train),                
                     gbm = gbm::gbm(model$call, 'bernoulli', data=cv_train_gbm[,c(colname_pred,model$response.name)], n.trees=model$gbm.call$best.trees, shrinkage=model$gbm.call$learning.rate, bag.fraction=model$gbm.call$bag.fraction, interaction.depth=model$gbm.call$tree.complexity),
                     maxnet = maxnet::maxnet(p= cv_train[,colname_species], data= cv_train[,colname_pred, drop=F]))
    
    # We make predictions for k-fold:
    if (class(model)[1]=='gbm') {
      cross_val_preds[ks==i] <- 
        dismo::predict.gbm(modtmp, cv_test[, colname_pred, drop=F], n.trees=model$gbm.call$best.trees, type="response")
    } else {
      cross_val_preds[ks==i] <- predictSDM(modtmp, cv_test[, colname_pred, drop=F])
    }
  }
  cross_val_preds
}

#--------------------------

# Inflated response curves

inflated_response=function(object,predictors,select.columns=NULL,label=NULL,len=50,lhsample=100,lwd=1,
                           ylab="Occurrence probabilities",method="stat3",disp="all",overlay.mean=T,
                           col.curves='grey',col.novel='grey',col.mean='black',lwd.known=2,lwd.mean=2,...){
  
  require(lhs)
  
  if (is.null(select.columns)) select.columns=seq_len(ncol(predictors))
  
  for (i in select.columns)
  {
    summaries=data.frame(matrix(0,6,ncol(predictors)))
    for (iz in 1:ncol(predictors)) {
      summaries[,iz]=summary(predictors[,iz])
    }
    if (method=="stat3") {
      summaries.j=as.matrix(summaries[c(1,4,6),-i],ncol=(ncol(predictors)-1));comb=min(lhsample,3^(ncol(predictors)-1));nc=3
    } else
      if (method=="stat6") {
        summaries.j=as.matrix(summaries[,-i],ncol=(ncol(predictors)-1));comb=min(lhsample,6^(ncol(predictors)-1));nc=6
      } else
        if (method=="mean") {
          summaries.j=as.matrix(summaries[4,-i],ncol=(ncol(predictors)-1));comb=1;nc=1;overlay.mean=F
        }
    
    dummy.j=as.matrix(predictors[1:len,-i],ncol=(ncol(predictors)-1))
    
    if (comb<lhsample) {
      mat=vector("list",ncol(dummy.j))
      for (m in 1:ncol(dummy.j)) mat[[m]]=1:nc
      mat=expand.grid(mat)
    } else {
      mat=round(qunif(lhs::randomLHS(lhsample,ncol(dummy.j)),1,nrow(summaries.j)),0)
    }
    
    if (is.null(label)) {
      label=names(predictors)
    }
    
    for (r in 1:nrow(mat))
    {
      for (j in 1:ncol(dummy.j))
      {
        dummy.j[,j]=as.vector(rep(summaries.j[mat[r,j],j],len))
      }
      
      dummy=data.frame(seq(min(predictors[,i]),max(predictors[,i]),length=len),dummy.j)
      names(dummy)[-1]=names(predictors)[-i]
      names(dummy)[1]=names(predictors)[i]
      
      curves <- predictSDM(object, dummy)
      
      # display all lines in same type
      if (disp=='all')
      {
        if (r==1)
        {
          if (i==1) plot(dummy[,names(predictors)[i]],
                         curves,type="l",ylim=c(0,1),xlab=label[i],ylab=ylab,
                         lwd=lwd,col=col.curves,...)
          else plot(dummy[,names(predictors)[i]],
                    curves,type="l",ylim=c(0,1),xlab=label[i],ylab="",lwd=lwd,col=col.curves,...)
        }
        else lines(dummy[,names(predictors)[i]],
                   curves,lwd=lwd,col=col.curves,...)
      }
      
      # highlight extrapolation to novel environmental conditions
      if (disp=='eo.mask')
      {
        novel=eo.mask(predictors,dummy)
        curves.known=curves
        curves.known[novel==1]=NA
        curves.novel=curves
        curves.novel[novel==0]=NA
        
        if (r==1)
        {
          if (i==1) {plot(dummy[,names(predictors)[i]],
                          curves.known,type="l",ylim=c(0,1),xlab=label[i],ylab=ylab,
                          lwd=lwd.known,col=col.curves,...)
            lines(dummy[,names(predictors)[i]],
                  curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
          else {plot(dummy[,names(predictors)[i]],
                     curves.known,type="l",ylim=c(0,1),xlab=label[i],ylab="",lwd=lwd.known,col=col.curves,...)
            lines(dummy[,names(predictors)[i]],
                  curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
        }
        else {lines(dummy[,names(predictors)[i]],
                    curves.known,lwd=lwd.known,col=col.curves,...)
          lines(dummy[,names(predictors)[i]],
                curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
      }
    }
    
    #-------------------------------------------------
    # now, this is for overlaying mean response curve
    if (overlay.mean==T)
    {
      dummy=predictors[1:len,]
      dummy[,i]=seq(min(predictors[,i]),max(predictors[,i]),length=len)
      for (j in 1:ncol(predictors))
      {
        if (j!=i) 
        {
          dummy[,j]=rep(mean(predictors[,j]),len)
        }
      }
      
      curves <- predictSDM(object, dummy)
      
      lines(dummy[,names(predictors)[i]],
            curves,lwd=lwd.mean,col=col.mean,...)
    }    
  }}


#----------------------

# partial response plots

partial_response=function(object,predictors,select.columns=NULL, label=NULL, len=50,
                          ylab="Occurrence probabilities", col='black',...){
  
  inflated_response(object,predictors,select.columns,label,len,method='mean',col.curves=col, ...)
}


#----------------------

# environmental overlap masks

eo_mask=function(traindata,newdata,nbin=5,type="EO")
{
  train.minima=apply(traindata,2,min)
  train.maxima=apply(traindata,2,max)
  
  train.ids=apply(apply(ceiling(apply(round(
    sweep(sweep(traindata, 2, train.minima, "-"), 2, train.maxima - train.minima, "/")*nbin,4),
    c(1,2),FUN=function(x){if(x==0)x=1 else x=x})),
    c(1,2),FUN=function(x){if(x<1)x=0 else if(x>nbin)x=nbin+1 else x=x}),1,paste,collapse=".")
  
  new.ids=apply(apply(ceiling(apply(round(
    sweep(sweep(newdata[,names(train.minima)], 2, train.minima, "-"), 2, train.maxima - train.minima, "/")*nbin,4),
    c(1,2),FUN=function(x){if(x==0)x=1 else x=x})),
    c(1,2),FUN=function(x){if(x<1)x=0 else if(x>nbin)x=nbin+1 else x=x}),1,paste,collapse=".")
  
  if (type=="ID") return(new.ids)
  else if (type=="EO") return(sapply(new.ids%in%train.ids,FUN=function(x){if(x==T) x=0 else if(x==F)x=1}))    
}  


