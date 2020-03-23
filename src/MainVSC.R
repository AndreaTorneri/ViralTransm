args <- commandArgs(trailingOnly = TRUE)

wd = args[1]
cores = args[2]
scen = args[3]
rho = args[4]
as.prop = args[5]
eta = args[6]
rq = args[7]
tfailtest = args[8]

setwd(wd)

#load packages
library(foreach)
library(doParallel)
# load functions
source("SimulationFunctions.R")
# Parameters
n<-500 #Population size
lambda<-12 #daily rate of social contacts
#rho<- 0.6 # proportion of symptomatic individual
#as.prop<-0.5
sigma<- 0.16 # proportion of severe cases (computed the 15.2.2020 from: https://docs.google.com/spreadsheets/d/1Z7VQ5xlf3BaTx_LBBblsW4hLoGYWnZyog3jqsS9Dbgc/edit#gid=957283529)
avg.inc<-5.2 # NEJM Li
mu.IP<-18 #28.6 #2.3 NEJM Li (SerialInt-Inc) -Anne Cori - SARS
#R.s<-2.5 # reproduction number symptomatic (mild)
#R.a<-2.5 #reproduction number asymptomatic
avg.dt<-1.6
#eta<-0.5
#tfailtest<-2
tdelaytest<-tfailtest
lambda.dec<-rq*lambda


cl<-makeCluster(as.numeric(cores))
registerDoParallel(cl)

nSim<-5000
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)

if (scen=="IAS"){
  VShed<-scen
  epi.outbreak<-foreach(i = 1:nSim) %dopar%{
    epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.aftersympt (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, as.prop = as.prop)
  }
}
if (scen=="IBS"){
  VShed<-scen
  epi.outbreak<-foreach(i = 1:nSim) %dopar%{
     epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.beforesympt.testing (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, as.prop = as.prop, tfailtest = tfailtest, tdelaytest = tdelaytest)  
  }
}
if (scen=="ITRB"){
  VShed<-scen
  epi.outbreak<-foreach(i = 1:nSim) %dopar%{
  epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.testing (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, as.prop = as.prop, tfailtest = tfailtest, tdelaytest = tdelaytest)  
  }
}

stopCluster(cl)


# Compute the final size and peak incidence for each simulation in which  >10% of individuals are ultimately infected
finalSize<-NULL
not.extinct<-NULL
for (i in 1:nSim){
  finalSize[i]<-length(c(which(epi.outbreak[[i]]$time.events[,2]==1.0),which(epi.outbreak[[i]]$time.events[,2]==1.1),which(epi.outbreak[[i]]$time.events[,2]==1.2)))
  if (finalSize[i]>round(n*0.1)){not.extinct<-c(not.extinct,i)}
}
casesPeak<-NULL
for (j in not.extinct) {
  time.events<-epi.outbreak[[j]]$time.events  
  epi.curve<-0
  for (i in 1:length(time.events[,1]) ){
    epi.curve[i]<-length(c(which(time.events[1:i,2]==1.1),which(time.events[1:i,2]==1.2),which(time.events[1:i,2]==1.0)))-length(which(time.events[1:i,2]==0.1))
  }
  casesPeak<-c(casesPeak, max(epi.curve))
}

FinSize<-finalSize[not.extinct]
PeakInc<-casesPeak
Ext.prob<-(1-length(not.extinct)/nSim)



name<-paste("Asymptomatic",n,"_propUnd",rho,"_Rs",R.s,"_Ra",R.a,"_lambda",lambda,"_Scenario",VShed,"_eta",eta,"_lambda_dec",lambda.dec,"_tfail",tfailtest,"_.RData", sep = "")

save(epi.outbreak, FinSize, PeakInc, Ext.prob, file = name)

