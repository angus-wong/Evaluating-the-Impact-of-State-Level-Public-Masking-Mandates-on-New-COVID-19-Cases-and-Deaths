get.data <- function(filepath,this.date, mobility=F,pre_policy=F,relativeGrowth=F, DC=F, NY=T){
  
  dt <- read.csv(filepath)  
  dt <- dt[dt$date==this.date,]
  
  # if excluding DC
  if(!DC){
    dt <- dt[dt$State_abb!='DC',]
  }
  
  # if excluding NY
  if(!NY){
    dt <- dt[dt$State_abb!='NY',]
  }
  
  dt <- cbind(dt, Percent.of.People.Below.Poverty=
                dt$Percentage.of.people.whose.income.in.the.past.12.months.is.below.the.poverty.level)

    # get covariates
  set0 <- get.W(mobility=mobility, pre_policy=pre_policy, relativeGrowth=relativeGrowth)
  
  dt<-cbind(dt, trump=as.numeric(dt$poli_party=='republican'))
  list(dt=dt, set0=set0)
}


get.W <- function(mobility=F, pre_policy=F, relativeGrowth=F){
  # List of covariates before screening
  set0 <- c('Over65_pct', 
            'Black_pct', 
            'Hispanic_pct',
            'Asian_pct', 
            'Mixed_pct', 
            'White_pct', 
            "MedianAge_pct",
            'trump',
            'Population.Density..people.per.square.kilometer.',
            'Urban_pct_2010',
            'Average.Household.Size',
            'MedianIncome',
            'Households..Income.Below.Poverty.Level....',
            'Percent.of.People.Below.Poverty',
            'PercentSmoke',
            'YesDiabete',
            'Drove.alone._pct',
            'Worked.at.home._pct',
            "Public.transportation._pct",
            "Taxi..motorcycle..or.other._pct",
            'Bicycle._pct',
            'Walked._pct',
            "Total.Population",
            "positive_normalized.30day.before",
            "positive_normalized.14day.before", 
            "positive_normalized.7day.before",
            "death_normalized.30day.before",
            "death_normalized.14day.before",
            "death_normalized.7day.before",                                       
            "totaltest_normalized.30day.before",
            "totaltest_normalized.14day.before",
            "totaltest_normalized.7day.before"
  )
  
  # List of Prior Policy covariates
  pp <- c("ever_SAHpolicy",
                      "ever_GRpolicy",  
                      "ever_RRpolicy",                                                                    
                      "ever_NEBCpolicy",                                                                   
                      "ever_OBCpolicy",                                                                    
                      "ever_BMpolicy",                                                                    
                      "ever_SMpolicy" )
  
  # List of Mobility covariates
  mob <- c("residential_percent_change_from_baseline.14day.before",
           "residential_percent_change_from_baseline.7day.before")
  
  rg <- c("relativeChange_7days.before",
          "relativeChange_14days.before",
          "relativeChange_30days.before"
          )
  
  if(pre_policy){
    set0 <- c(set0, pp)  
  }
  if(mobility){
    set0 <- c(set0, mob)  
  }
  if(relativeGrowth){
    set0 <- c(set0, rg)  
  }
  set0
}

get.table1 <- function(W){
  summarize.me <- function(X){
    q2 <- median(X)
    q1 <- quantile(X, .25)
    q3 <- quantile(X, .75)
    paste( round(q2,1), ' (', round(q1,1), ', ' , round(q3,1), ')', sep='')
  }
  out <- NULL
  for(w in 1:ncol(W)){
    out <- cbind(out, summarize.me(W[,w]))
  }
  colnames(out) <- names(W)
  if("trump" %in% colnames(out))
    out[,'trump'] <- round(mean(W$trump)*100,1)
  if("ever_SAHpolicy" %in% colnames(out))
    out[,'ever_SAHpolicy'] <- round(mean(W$ever_SAHpolicy)*100,1)
  if("ever_GRpolicy" %in% colnames(out))
    out[,'ever_GRpolicy'] <- round(mean(W$ever_GRpolicy)*100,1)
  if("ever_RRpolicy" %in% colnames(out))
    out[,'ever_RRpolicy'] <- round(mean(W$ever_RRpolicy)*100,1)
  if("ever_NEBCpolicy" %in% colnames(out))
    out[,'ever_NEBCpolicy'] <- round(mean(W$ever_NEBCpolicy)*100,1)
  if("ever_OBCpolicy" %in% colnames(out))
    out[,'ever_OBCpolicy'] <- round(mean(W$ever_OBCpolicy)*100,1)
  if("ever_BMpolicy" %in% colnames(out))
    out[,'ever_BMpolicy'] <- round(mean(W$ever_BMpolicy)*100,1)
  if("ever_SMpolicy" %in% colnames(out))
    out[,'ever_SMpolicy'] <- round(mean(W$ever_SMpolicy)*100,1)
  out <- cbind( N=nrow(W), out )
  t(out)
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

  if(do.print){
    print('Primary approach: estimate of expected outcome under txt')
    print(txt$fit$Q[[1]])
    print('Primary approach: estimate of expected outcome under control')
    print(con$fit$Q[[1]])

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
