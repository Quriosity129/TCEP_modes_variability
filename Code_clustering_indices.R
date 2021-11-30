
require(glmnet)
require(lubridate)

P__disp <- function(f){
  return(pchisq(sum(residuals(f,type="pearson")**2),df=f$df.residual,lower.tail=FALSE))
}

myDecluster <- function(v,x,r){
  
  # v: data series
  # x: binary time series of event occurrences (raw)
  # r: run length for declustering
  y <- x
  w <- rle(x)
  # Set runs of zeros with length<r to 1
  final_pos <- cumsum(w$lengths)
  for (uu in 1:(length(w$lengths)-1)){
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  if (x[length(x)]==1){
    uu = w$lengths[length(w$lengths)]
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  w <- rle(y)
  # maximum element of each run of '1's
  final_pos <- cumsum(w$lengths)
  out_max_idx = out_max = out_clust_len = out_clust_evts = out_clust_tot = 0*1:sum(w$values==1)
  index <- 1
  for (i in 1:length(final_pos)){
    if (w$values[i]==1){
      vec = v[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      vec2 = x[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      out_max[index] = max(vec)
      out_max_idx[index] = which.max(vec)+final_pos[i]-w$lengths[i]
      out_clust_len[index] = length(vec)
      out_clust_evts[index] = sum(vec2)
      out_clust_tot[index] = sum(vec)
      index = index+1
    }
  }
  return(list(out_max_idx,out_max,out_clust_len,out_clust_tot))
}

################################################################################
## STEPWISE ##
###############################################################################

fitGLM <- function(pr_data,d,season,mat_covs,L=21,r){

  ## INPUTS ##
  
  # pr_data: daily precipitation data (typically mm/day)
  # d: date arraycorresponding to the precipitation data
  # season: selected season to model extreme precipitation counts (1: DJF; 2: MAM; 3: JJA and 4: SON)
  # mat_covs = matrix of covariates (daily, same time coverage as pr_data)
  # L: averaging time window (days)
  # r: declustering window (days)

  mth_sel <- list(c(12,1,2),3:5,6:8,9:11)[[season]]
  NN <- floor(length(d)/L)
  stepVec <- c(rep(1:NN,each=L),rep(NN+1,length(d)%%L))
  monthVec <- aggregate(month~step,data.frame(step=stepVec,month=month(d)),median)$month

  # Determine extreme precipitation events
  ww <- 0*pr_data
  for (mm in 1:12){
    ww[month(d)==mm] <- 1*(pr_data>=quantile(pr_data[month(d)==mm],0.99,na.rm=T))[month(d)==mm]
  }
  res <- myDecluster(pr_data,ww,r=r)
  X <- 0*pr_data
  X[res[[1]]] <- 1

  # Aggregate data over selected time window
  data <- as.data.frame(cbind(stepVec,X,mat_covs))
  colnames(data) <- c('step','cnt',colnames(mat_covs))
  data_final <- aggregate(.~step,data,mean)
  data_final$cnt <- aggregate(.~step,data,sum)$cnt

  # Fit GLM
  data_train <- data[monthVec%in%mth_sel,]
  listAIC <- vector()
  modelVars <- vector()
  remainingVars <- 2:ncol(data_train)
  currentAIC <- AIC(glm(data_train$cnt~1,family=poisson("log")))
  listAIC[1] <- AIC(glm(data_train$cnt~1,family=poisson("log")))
  foundBetterModel = TRUE

  while (foundBetterModel){

    # Loop on all remaining variables
    foundBetterModel <- F
    minAIC <- tail(listAIC,1) # current minimum AIC value

    for (var in remainingVars){

      for (varBis in modelVars){
        fit <- glm(formula=cnt~.,data=data_train[,c(1,modelVars[!modelVars==varBis],var)],family=poisson("log"))
        if (AIC(fit)<minAIC){
          foundBetterModel <- T
          whichModel <- c(modelVars[!modelVars==varBis],var)
          minAIC <- AIC(fit)
        }
      }
      fit <- glm(formula=cnt~.,data=data_train[,c(1,modelVars,var)],family=poisson("log"))
      if (is.infinite(AIC(fit))){
        break
      }
      # Select variable with best improvement
      if (AIC(fit)<minAIC){
        foundBetterModel <- T
        whichModel <- c(modelVars,var)
        selVarAIC <- var
        minAIC <- AIC(fit)
      }
    }

    # Update
    if (foundBetterModel){
      listAIC <- c(listAIC,minAIC)
      modelVars <- whichModel
      remainingVars <- remainingVars[!remainingVars%in%whichModel]
    } else{
      break
    }
  }
  preds = colnames(data_final)[modelVars]
  if (length(preds)==0){
    model.formula <- as.formula("cnt~1")
    fit.train <- glm(formula=model.formula,data=data_train,family=poisson("log"))
  }else{
    model.formula <- as.formula(paste("cnt~",paste(preds,collapse="+"),sep=""))
    fit.train <- glm(formula=model.formula,data=data_train,family=poisson("log"))
  }
  x <- fit.train$coefficients[-1] # remove intercept
  ii <- 1
  coeffs_result <- 0*1:(ncol(data_train)-1)
  for (pp in preds){
    coeffs_result[which(colnames(data_train)==pp)-1] <- x[ii]
    ii <- ii+1
  }
  colnames(coeffs_result) <- colnames(data_train)
  dev_ratio <- 1-fit.train$deviance/fit.train$null.deviance
  p_value <- P__disp(fit.train)

  return(list(coeffs_result,dev_ratio,p_value))
}