library(MSGARCH)
library(e1071)
library(quantmod)
library(GAS)
library(TTR)
library(PerformanceAnalytics)
library(tictoc)
library(Rcpp)
library(purrr)
library(foreach)
library(doParallel)
library(furrr)
library(tidyr)
library(dplyr)
library(writexl)


parallel::detectCores()
my.cluster = parallel::makeCluster(6, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

##Get Data and Clean and De-mean Data
getSymbols('^GSPC', src = 'yahoo', from = '2002-12-26', to = ' 2016-11-18')
close_sp = GSPC$GSPC.Adjusted
sum(is.na(close_sp))
rsp = diff(log(close_sp), lag = 1)*100
AR_filter = arima(rsp, order = c(1,0,0))
mean = AR_filter$coef[2]
de_rsp = rsp-as.numeric(mean)
rsp = as.numeric(rsp)
rsp[is.na(rsp)]=0


##Creating Mesh
zl = -100
zu = 100
M = 1000
x <- seq(zl,zu,(zu-zl)/M)
dz = 1000/(zu-zl)



##Setup Different Models
combinations = crossing(models = c("sGARCH","gjrGARCH"),
                        distributions = c("norm","std","sstd"))
plan(multisession, workers = 6)
len_m = length(combinations[[1]])



##Iterate through Data to get VaR and ES

windowsize = 1500
periodlength = 500
k.update = 10
alpha = c(0.01, 0.05)

VaR.ML <-array(data=NA,dim=c(periodlength,2*len_m,2))
VaR.Bayes <-array(data=NA,dim=c(periodlength,2*len_m,2))
ES.ML <-array(data=NA,dim=c(periodlength,2*len_m,2))
ES.Bayes <-array(data=NA,dim=c(periodlength,2*len_m,2))
modelML2.fit <- vector(mode = "list", length = len_m)
modelML1.fit <- vector(mode = "list", length = len_m)
modelBayes2.fit <- vector(mode = "list", length = len_m)
modelBayes1.fit <- vector(mode = "list", length = len_m)
crps.ML<-array(data=NA,dim=c(periodlength,(2*len_m)))
crps.MCMC<-array(data=NA,dim=c(periodlength,(2*len_m)))
y.ots <- matrix(NA, nrow = periodlength, ncol = 1)

tic()
for(i in 1:periodlength){
  tic()
  cat("Backtest - Iteration: ", i, "\n")
  tryCatch({
    s.ML=0
    s.MCMC=0
    y.its<-rsp[i:(windowsize+i-1)]
    y.ots[i]<-rsp[windowsize+i]
    if(i==1||i%%k.update==1){
      ## K=2 MS model
      modelML2.fit <- future_map2(combinations$models,combinations$distributions, function(x,y){
        MSgarch_settings <- CreateSpec(variance.spec=list(model = x),
                                       distribution.spec=list(distribution = y),
                                       switch.spec = list(K = 2))
        try(FitML(spec = MSgarch_settings,data = y.its,ctr = list(do.se = FALSE)))
      },.options = furrr_options(seed = T)
      )
      #respec ML
      for (j in 1:len_m){
        modelML2.fit[[j]][["spec"]] = CreateSpec(variance.spec=list(model = combinations[[1]][[j]]),
                                                distribution.spec=list(distribution = combinations[[2]][[j]]),
                                                switch.spec = list( K = 2))
      }
      #fit MC
      modelBayes2.fit <- future_map2(combinations$models,combinations$distributions, function(x,y){
        MSgarch_settings <- CreateSpec(variance.spec=list(model = x),
                                       distribution.spec=list(distribution = y),
                                       switch.spec = list( K = 2))
        try(FitMCMC(spec = MSgarch_settings,data = y.its,ctr = list(do.se = FALSE)))
      },.options = furrr_options(seed = T)
      )
      #respec MC
      for (j in 1:len_m){
        modelBayes2.fit[[j]][["spec"]] = CreateSpec(variance.spec=list(model = combinations[[1]][[j]]),
                                                   distribution.spec=list(distribution = combinations[[2]][[j]]),
                                                   switch.spec = list( K = 2))
      }
      
      
      ## K = 1 Single regime 
      modelML1.fit <- future_map2(combinations$models,combinations$distributions, function(x,y){
        MSgarch_settings <- CreateSpec(variance.spec=list(model = x),
                                       distribution.spec=list(distribution = y),
                                       switch.spec = list(K = 1))
        try(FitML(spec = MSgarch_settings,data = y.its,ctr = list(do.se = FALSE)))
      },.options = furrr_options(seed = T)
      )
      #respec ML
      for (j in 1:len_m){
        modelML1.fit[[j]][["spec"]] = CreateSpec(variance.spec=list(model = combinations[[1]][[j]]),
                                                distribution.spec=list(distribution = combinations[[2]][[j]]),
                                                switch.spec = list( K = 1))
      }
      #fit MC
      modelBayes1.fit <- future_map2(combinations$models,combinations$distributions, function(x,y){
        MSgarch_settings <- CreateSpec(variance.spec=list(model = x),
                                       distribution.spec=list(distribution = y),
                                       switch.spec = list( K = 1))
        try(FitMCMC(spec = MSgarch_settings,data = y.its,ctr = list(do.se = FALSE)))
      },.options = furrr_options(seed = T)
      )
      #respec MC
      for (j in 1:len_m){
        modelBayes1.fit[[j]][["spec"]] = CreateSpec(variance.spec=list(model = combinations[[1]][[j]]),
                                                   distribution.spec=list(distribution = combinations[[2]][[j]]),
                                                   switch.spec = list( K = 1))
      }
      modelBayes.fit <-c(modelBayes1.fit, modelBayes2.fit)
      modelML.fit <-c(modelML1.fit, modelML2.fit)
    }
    
    ## estimate predict pdf for CRSP and estimate Risk function
    for(j in 1:(2*len_m)){
      predpdfML <- PredPdf(object = modelML.fit[[j]], x = x, log = FALSE, do.its = FALSE, nahead = 1, ctr = list(nsim = 500L))
      predpdfMCMC <- PredPdf(object = modelBayes.fit[[j]], x = x, log = FALSE, do.its = FALSE, nahead = 1, ctr = list(nsim = 500L))
      Risk1<-Risk(modelML.fit[[j]]$spec, par=modelML.fit[[j]]$par,data=y.its,n.ahead=1,alpha=alpha, do.its=FALSE)
      Risk2<-Risk(modelBayes.fit[[j]]$spec, par=modelBayes.fit[[j]]$par,data=y.its,n.ahead=1,alpha=alpha, do.its=FALSE)
      
      VaR.ML[i,j,]<-Risk1$VaR
      VaR.Bayes[i,j,]<-Risk2$VaR
      ES.ML[i,j,]<-Risk1$ES
      ES.Bayes[i,j,]<-Risk2$ES
      for (k in 1:1000){
        zm = zl + k * (zu - zl) / M
        F.ML = sum(predpdfML[1:round((zm-zl)*dz)])/dz
        F.MCMC = sum(predpdfMCMC[1:round((zm-zl)*dz)])/dz
        s.ML = s.ML +  (1-pnorm(zm)) * (F.ML - 1*(y.ots[i] < zm))*2
        s.MCMC = s.MCMC +  (1-pnorm(zm)) * (F.MCMC - 1*(y.ots[i] < zm))*2
      }
      crps.ML[i,j] = (zu-zl)/(M-1) * s.ML
      crps.MCMC[i,j] = (zu-zl)/(M-1) * s.MCMC

    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  toc()
}
toc()


## NA check and solve
crps1.ML = crps.ML
crps1.MCMC = crps.MCMC
for(j in 1:(2*len_m)){
  a = mean(crps1.ML[,j][!is.na(crps1.ML[,j])])
  b = mean(crps1.MCMC[,j][!is.na(crps1.MCMC[,j])])
  crps1.ML[,j] = crps1.ML[,j] %>% replace(is.na(crps1.ML[,j]), a)
  crps1.MCMC[,j] = crps1.MCMC[,j] %>% replace(is.na(crps1.MCMC[,j]), b)
}

VaR1.ML = VaR.ML
VaR1.Bayes =VaR.Bayes
ES1.ML = ES.ML
ES1.Bayes = ES.Bayes

for(j in 1:(2*len_m)){
  for (k in 1:2){
    a = mean(VaR1.ML[, j, k][!is.na(VaR1.ML[, j, k])])
    b = mean(VaR1.Bayes[, j, k][!is.na(VaR1.Bayes[, j, k])])
    c = mean(ES1.ML[, j, k][!is.na(ES1.ML[, j, k])])
    d = mean(ES1.Bayes[, j, k][!is.na(ES1.Bayes[, j, k])])
    VaR1.ML[, j, k] = VaR1.ML[, j, k] %>% replace(is.na(VaR1.ML[, j, k]), a)
    VaR1.Bayes[, j, k] = VaR1.Bayes[, j, k] %>% replace(is.na(VaR1.Bayes[, j, k]), b)
    ES1.ML[, j, k] = ES1.ML[, j, k] %>% replace(is.na(ES1.ML[, j, k]), a)
    ES1.Bayes[, j, k] = ES1.Bayes[, j, k] %>% replace(is.na(ES1.Bayes[, j, k]), b)
    
  }
}


#wCRPS test
sumCRPS.ML<-vector(length=length(2*len_m))
sumCRPS.MCMC<-vector(length=length(2*len_m))
for(j in 1:(2*len_m)){
  sumCRPS.ML[j]=sum(crps1.ML[,j])
  sumCRPS.MCMC[j]=sum(crps1.MCMC[,j])
}



##test
UC.pval.ML = array(data = NA, dim = c((2*len_m), 2))
UC.pval.MC = array(data = NA, dim = c((2*len_m), 2))
DQ.pval.ML = array(data = NA, dim = c((2*len_m), 2))
DQ.pval.MC = array(data = NA, dim = c((2*len_m), 2))
QL.ML = array(data = NA, dim = c((2*len_m), 2))
QL.MC = array(data = NA, dim = c((2*len_m), 2))

FZL.ML = array(data = NA, dim = c((2*len_m), 2))
FZL.MCMC = array(data = NA, dim = c((2*len_m), 2))

##for alpha = 0.05




test_ML_5 <- foreach(j=1:(2*len_m)) %dopar% GAS::BacktestVaR(data = y.ots, VaR = VaR1.ML[, j, 2],alpha = 0.05)
test_MC_5 <- foreach(j=1:(2*len_m)) %dopar%  GAS::BacktestVaR(data = y.ots, VaR = VaR1.Bayes[, j, 2],alpha = 0.05)
test_ML_1 <- foreach(j=1:(2*len_m)) %dopar% GAS::BacktestVaR(data = y.ots, VaR = VaR1.ML[, j, 1],alpha = 0.01)
test_MC_1 <- foreach(j=1:(2*len_m)) %dopar%  GAS::BacktestVaR(data = y.ots, VaR = VaR1.Bayes[, j, 1],alpha = 0.01)

for (j in 1:length(modelML.fit)){

  UC.pval.ML[j,2] = as.numeric(test_ML_5[[j]]$LRuc[2])
  DQ.pval.ML[j,2] = as.numeric(test_ML_5[[j]]$DQ[2])
  QL.ML[j,2] = as.numeric(test_ML_5[[j]]$Loss[1])
  
  UC.pval.MC[j,2] = as.numeric(test_MC_5[[j]]$LRuc[2])
  DQ.pval.MC[j,2] = as.numeric(test_MC_5[[j]]$DQ[2])
  QL.MC[j,2] = as.numeric(test_MC_5[[j]]$Loss[1])
  
  UC.pval.ML[j,1] = as.numeric(test_ML_1[[j]]$LRuc[2])
  DQ.pval.ML[j,1] = as.numeric(test_ML_1[[j]]$DQ[2])
  QL.ML[j,1] = as.numeric(test_ML_1[[j]]$Loss[1])
  
  UC.pval.MC[j,1] = as.numeric(test_MC_1[[j]]$LRuc[2])
  DQ.pval.MC[j,1] = as.numeric(test_MC_1[[j]]$DQ[2])
  QL.MC[j,1] = as.numeric(test_MC_1[[j]]$Loss[1])
  
}

 
for(i in 1:length(modelML.fit)){
  FZL.ML[,1] = mean(unlist(FZLoss(data = y.ots, VaR = VaR1.ML[,j,1], ES = ES1.ML[, j,1],
                           alpha = .01)))
  FZL.MCMC[,1] = mean(unlist(FZLoss(data = y.ots, VaR = VaR1.Bayes[, j,1], ES = ES1.Bayes[, j,1],
                             alpha = .01)))
  FZL.ML[,2] = mean(unlist(FZLoss(data = y.ots, VaR = VaR1.ML[,j,2], ES = ES1.ML[, j,2],
                           alpha = .05)))
  FZL.MCMC[,2] = mean(unlist(FZLoss(data = y.ots, VaR = VaR1.Bayes[, j,2], ES = ES1.Bayes[, j,2],
                             alpha = .05)))
}


dfcrps = dplyr::bind_cols(as.data.frame(sumCRPS.ML), as.data.frame(sumCRPS.MCMC))

dfFZL = dplyr::bind_cols(as.data.frame(FZL.ML),as.data.frame(FZL.MCMC))

dftest = dplyr::bind_cols(as.data.frame(UC.pval.ML),as.data.frame(DQ.pval.ML),as.data.frame(QL.ML),

                 as.data.frame(UC.pval.MC),as.data.frame(DQ.pval.MC),as.data.frame(QL.MC))

write_xlsx(dfcrps, "dfcrps.xlsx")
write_xlsx(dfFZL, "dfFZL.xlsx")
write_xlsx(dftest, "dftest.xlsx")
















