args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
scen = args[3] #scenario: IAS, IBS  or IBTBS
cat(",scen=",scen)
rho = as.numeric(args[4]) #proportion of diagnosed
cat(",rho=",rho)
as.prop = as.numeric(args[5]) #proportion of asymptomatic (among the undiagnosed)
cat(",as.prop=",as.prop)
eta = as.numeric(args[6]) # probability of successfuly trace-back an individual
cat(",eta=",eta)
rq = as.numeric(args[7]) # decrease in contact rate due to quarantine
cat(",rq=",rq)
tfailtest = as.numeric(args[8]) # time since infection at which the test results positive
cat(",tfailtest=",tfailtest)
R.a = as.numeric(args[9]) # reproduction number asymptomatic
cat(",R.a=",R.a)
R.s = as.numeric(args[10]) # reproduction number symptomatic
cat(",R.s=",R.s)

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
sigma<- 0.16 # proportion of severe cases 
avg.inc<-5.2 # NEJM Li
mu.IP<-18 #mean length symptomatic period
#R.s<-2.5 # reproduction number symptomatic (mild)
#R.a<-2.5 #reproduction number asymptomatic
avg.dt<-1.6 # mean length time to diagnosed (after symptoms onset)
#eta<-0.5
#tfailtest<-2
tdelaytest<-tfailtest # time at which the second test takes place, if the first one fails
lambda.dec<-rq*lambda # contact rate in quarantine


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
if (scen=="IBTBS"){
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
  finalSize[i]<-length(c(which(epi.outbreak[[i]]$time.events[,2]==1.01),which(epi.outbreak[[i]]$time.events[,2]==1.02),which(epi.outbreak[[i]]$time.events[,2]==1.1),which(epi.outbreak[[i]]$time.events[,2]==1.2)))
  if (finalSize[i]>round(n*0.1)){not.extinct<-c(not.extinct,i)}
}
casesPeak<-NULL
for (j in not.extinct){
  time.events<-epi.outbreak[[j]]$time.events  
  epi.curve<-0
  for (i in 1:length(time.events[,1]) ){
    epi.curve[i]<-length(c(which(time.events[1:i,2]==1.1),which(time.events[1:i,2]==1.2),which(time.events[1:i,2]==1.01),which(time.events[1:i,2]==1.02)))-length(which(time.events[1:i,2]==0.1))
  }
  casesPeak<-c(casesPeak, max(epi.curve))
}

FinSize<-finalSize[not.extinct]
PeakInc<-casesPeak
Ext.prob<-(1-length(not.extinct)/nSim)



name<-paste("Asymptomatic",n,"_propUnd",rho,"_Rs",R.s,"_Ra",R.a,"_lambda",lambda,"_Scenario",VShed,"_eta",eta,"_lambda_dec",lambda.dec,"_tfail",tfailtest,"_.RData", sep = "")

setwd(out)
save(epi.outbreak, FinSize, PeakInc, Ext.prob, file = name)

