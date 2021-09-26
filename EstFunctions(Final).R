RUN <- function(filepath, this.date, mobility, pre_policy, relativeGrowth,
                DC, NY, screen, SL.library, A, Y){
  
  DT <- get.data(filepath = filepath, this.date=this.date, mobility=mobility, pre_policy=pre_policy,DC=DC, NY=NY)
  dt <- DT$dt
  
  dt[dt$everImplement==0,A] <- 1
  set0 <- DT$set0
  
  if(!DC){
    dt <- dt[dt$State_abb!='DC',]
  }
  
  # if excluding NY
  if(!NY){
    dt <- dt[dt$State_abb!='NY',]
  }
  
  if(screen){
    # exclude NYC from screening process
    temp <- dt[dt$State!='New York City',]
    set1 <- set0[screen.corP(Y=temp[,Y], X=temp[,set0] ) ]
    
  }
  
  pscore_unadj <- fit.g(dt=dt, W=set0, A=A, SL.library=SL.library)
  pscore_TMLE <- fit.g(dt=dt, W=set1, A=A, SL.library=SL.library)
  pscore_all<- as.data.frame(cbind(dt %>% select(State),dt %>% select(A), pscore_unadj,pscore_TMLE))
  
  table_unadj <- data.frame(cbind(
    all = get.table1(W=dt[,set0]),
    txt = get.table1(W=dt[dt[,A]==1,set0]),
    con = get.table1(W=dt[dt[,A]==0,set0])
  ))
  colnames(table_unadj) <- c('All', 'Early', 'Delayed')
  
  table_TMLE <- data.frame(cbind(
    all = get.table1(W=dt[,set1]),
    txt = get.table1(W=dt[dt[,A]==1,set1]),
    con = get.table1(W=dt[dt[,A]==0,set1])
  ))
  colnames(table_TMLE) <- c('All', 'Early', 'Delayed')
  
  result_unadj <- run.ltmle(dt=dt,  W=NULL, A=A, Y=Y,pscore=NULL, SL.library=NULL)
  result_TMLE <- run.ltmle(dt=dt, W=set1, A=A, Y=Y, pscore=pscore_TMLE, SL.library=SL.library)
  result_all <- as.data.frame(rbind(result_unadj, result_TMLE),row.names = c("Unadjusted","TMLE"))
  
  outcomeDays_suffix = sub(".*_", "", Y)
  
  output <- list(table_unadj, table_TMLE, pscore_all, result_unadj, result_TMLE)
  names(output) <- c("Covariate_Summary_Unadjusted",
                     "Covariate_Summary_TMLE",
                     "Propensity_Score",
                     paste0("Result_Unadj","_",outcomeDays_suffix),
                     paste0("Result_TMLE","_",outcomeDays_suffix))
  
  output
}

# ONLY ESTIMATE THE PROPENSITY SCORE ONCE

fit.g <- function(dt, W, A, SL.library){
  set.seed(1)
  fit.pscore <- SuperLearner(Y=dt[,A], X=dt[,W], family=binomial(), SL.library=SL.library,
               cvControl=list(V=10))
  print('SL for pscore')
  print(fit.pscore)
  pscore <- fit.pscore$SL.predict
  print('Summary pscore dist')
  print(table(pscore) )
  pscore
}

do.analyses <- function(dt, W, A, Y, pscore, SL.library, NY.sens=F, NAME, TYPE){
  #primary adjusted
  primary <- run.ltmle(dt=dt, W=W, A=A, Y=Y, pscore=pscore, SL.library=SL.library,
                       do.print=(NAME=='All' &TYPE=='Rate'))
  primary <- c(primary, who=NAME, type=TYPE, est='adj')
  # unadjusted 
  unadj <- run.ltmle(dt=dt,  W=NULL, A=A, Y=Y,pscore=NULL, SL.library=NULL)
  unadj <- c(unadj, who=NAME,  type=TYPE, est='unadj')
  
  # sensitivity - drop NY
  if(NY.sens){
    dt.noNY <- dt[dt$State_abb!='NY',]
    set.seed(1)
    sens <- run.ltmle(dt=dt.noNY, W=W, A=A, Y=Y, SL.library=SL.library)
    OUT <- data.frame(rbind(primary=primary,  unadjusted=unadj, sens.no.NY=sens))
    
  }else {
    OUT <- data.frame(rbind(adj=primary,  unadj=unadj)) 
    
  }
  rownames(OUT) <- paste( NAME, rownames(OUT), sep='.')
  OUT
}


run.ltmle <- function(dt, W=NULL, A, Y, pscore, SL.library=NULL,  do.print=F){
  
  if(is.null(W)){
    dt.est <- dt[,c(A,Y)]
  } else{
    dt.est <- dt[,c(W,A,Y)]
  }
  set.seed(1)
  txt <- ltmle(data=dt.est, Anodes=A, Ynodes=Y, abar=1, 
               gform=pscore,
               SL.library=SL.library ,  SL.cvControl = list(V=10), estimate.time=F)
  set.seed(1)
  con <- ltmle(data=dt.est, Anodes=A, Ynodes=Y, abar=0, 
               gform=pscore,
               SL.library=SL.library ,  SL.cvControl = list(V=10), estimate.time=F)
  
  # est <-  ltmle(data=dt.est, Anodes=A, Ynodes=Y, abar=list(1,0), 
  #               gform=as.array( array(pscore, dim=c(nrow(dt.est),2))) , #variance.method='ic',
  #               SL.library=SL.library ,  SL.cvControl = list(V=10), estimate.time=F)
  if(do.print){
    print('Primary approach: estimate of expected outcome under txt')
    print(txt$fit$Q[[1]])
    print('Primary approach: estimate of expected outcome under control')
    print(con$fit$Q[[1]])
    
    # print('RESULTS')
  }
  out <- get.inference(txt=txt, con=con)
  
  out
}

get.inference <- function(txt, con){
  
  psi.1 <- txt$estimates['tmle']
  psi.0 <- con$estimates['tmle']
  
  IC.1 <- txt$IC$tmle
  IC.0 <- con$IC$tmle
  IC.diff <- IC.1 - IC.0
  # going after aRR, then get IC estimate on log scale
  #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
  IC.ratio <- 1/psi.1*IC.1 - 1/psi.0*IC.0
  
  txt <- get.CI(psi.hat=psi.1, IC=IC.1,  start='txt')
  con <- get.CI(psi.hat=psi.0, IC=IC.0,  start='con')
  RD <- get.CI(psi.hat=(psi.1-psi.0), IC=IC.diff,  start='RD')
  RR <- get.CI(psi.hat=log(psi.1/psi.0), IC=IC.ratio, do.relative=T, start='RR')
  out <- data.frame(cbind(txt, con, RD, RR))
  out
}
get.CI <- function(psi.hat, IC, do.relative=F, start){
  
  J <- length(IC)
  var.IC <- var(IC)/J
  
  # cutoff based on t-dist for testing and CI	
  cutoff <- qt(0.05/2, df=(J-2), lower.tail=F)
  
  # standard error (square root of the variance)
  se<- sqrt(var.IC)
  
  # test statistic (if goal=aRR then on the transformed scale)
  tstat <- psi.hat/se
  pval<- 2*pt(abs(tstat), df=(J-2), lower.tail=F) 
  
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  if(do.relative){
    psi.hat<- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  out<- data.frame(pt=psi.hat,  CI.lo, CI.hi ,se, pval)
  colnames(out) <- paste(start, colnames(out), sep='.')
  out
}







do.race.calcs <- function(dt.temp, GROUP, set0, A, SL.library){
  
  
  if(GROUP=='Black'){
    Y <- dt.temp[,'Deaths.Race.Known.Black.']
    perc <- dt.temp[,'Black_pct']
  }else if(GROUP=='White'){
    Y <- dt.temp[,'Deaths.Race.Known.White.']
    perc <- dt.temp[,'White_pct']
  } else if(GROUP=='LatinX'){
    Y <- dt.temp[,'Deaths.Ethnicity.Known.Hispanic']
    perc <- dt.temp[,'Hispanic_pct']
  }
  dt.temp <- cbind(dt.temp, Y, perc)
  
  keep <- !is.na(Y)
  dt.temp <- dt.temp[keep, ]
  
  # total reported divided by proportion of the population 
  Y.rate <-  dt.temp$Y/(dt.temp$perc/100*dt.temp$Total.Population)*100000
  #note cannot do residual analyses, bc most states did not report race-specific death early on

  # if PT, then sensitivitiy analysis accounting for days since first COVID outcome
  if(!is.na(dt.temp[1, 'PT'] )){
    Y.rate <- Y.rate/dt.temp$PT*100
  }
  
  
  dt.temp <- cbind(dt.temp, Y.rate)
  rate <- do.analyses(dt=dt.temp, W=set0, A=A, Y='Y.rate',pscore=as.matrix(dt.temp$pscore), 
                      SL.library=SL.library, NAME=GROUP, TYPE='Rate')

  list(rate=rate, states= as.character(dt.temp$State_abb), Nstates=length(dt.temp$State_abb))
}





#SL.earth functions

.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

SL.earth <- function(Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, nk = max(21, 2*ncol(X) + 1), pmethod = "backward", nfold = 0, ncross = 1, minspan = 0, endspan = 0,...) {
  .SL.require('earth')
  if(family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan)
  }
  if(family$family == "binomial") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan, glm = list(family = binomial))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}

# 
predict.SL.earth <- function(object, newdata,...) {
  .SL.require('earth')
  pred <- predict(object$object, newdata = newdata, type = "response")
  return(pred)
}
