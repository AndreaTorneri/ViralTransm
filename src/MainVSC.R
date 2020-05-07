args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
scen = args[3] #scenario: IAS, IBS  or IBTBS
cat(",scen=",scen)
rho = as.numeric(args[4]) #proportion of symptomatic
cat(",rho=",rho)
s.undiagn = as.numeric(args[5]) #proportion of mild-symptomatic that are undiagnosed
cat(",s.undiagn=",s.undiagn)
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
n = as.numeric(args[11]) # population size
cat(",n=",n)
seeds = as.numeric(args[12]) # initials infectives
cat(",seeds=",seeds)
infm.t = as.numeric(args[13]) # infectiousness decrease type due to antiviral injection
cat(",infm.t=",infm.t)



#print(paste0("scen",scen,"_rho",rho,"_sund",s.undiagn, "_eta",eta,"_rq",rq, "_tfail", tfailtest, "_R.a",R.a,"_R.s",R.s,"_n",n,"_inSeed",seeds))

#load packages
library(foreach)
library(doParallel)
# load functions
source("SimulationFunctions.R")
# Fixed Parameters
lambda<-12 #daily rate of social contacts
sigma<- 0.16 # proportion of severe cases 
avg.inc<-5.2 # NEJM Li
mu.IP<-18 #mean length symptomatic period
avg.dt<-3.8 # mean length time to diagnosed (after symptoms onset)
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
    epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.aftersympt (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, s.undiagn = s.undiagn, seeds=seeds)
  }
}
if (scen=="IBS"){
  VShed<-scen
  epi.outbreak<-foreach(i = 1:nSim) %dopar%{
     epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.beforesympt.testing (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, s.undiagn = s.undiagn, tfailtest = tfailtest, tdelaytest = tdelaytest, seeds=seeds)  
  }
}
if (scen=="IBTBS"){
  if (infm.t==1){
    VShed<-paste0(scen,"_INFMT_Exp")
    epi.outbreak<-foreach(i = 1:nSim) %dopar%{
      epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.testing (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, s.undiagn = s.undiagn, tfailtest = tfailtest, tdelaytest = tdelaytest, seeds=seeds)  
    }
  }
  if (infm.t==2){
    VShed<-paste0(scen,"_INFMT_Log")
    epi.outbreak<-foreach(i = 1:nSim) %dopar%{
      epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.testing.log (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, s.undiagn = s.undiagn, tfailtest = tfailtest, tdelaytest = tdelaytest, seeds=seeds)  
    }
  }
  if (infm.t==3){
    VShed<-paste0(scen,"_INFMT_K")
    epi.outbreak<-foreach(i = 1:nSim) %dopar%{
      epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.testing.k (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta, s.undiagn = s.undiagn, tfailtest = tfailtest, tdelaytest = tdelaytest, seeds=seeds)  
    }
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
for (j in not.extinct){
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



name<-paste(VShed,n,"_propAs",rho,"_SU",s.undiagn,"_Rs",R.s,"_Ra",R.a,"_lambda",lambda,"_eta",eta,"_lambda_dec",lambda.dec,"_tfail",tfailtest,"InSeeds",seeds,"_.RData", sep = "")

setwd(out)
save(epi.outbreak, FinSize, PeakInc, Ext.prob,not.extinct, file = name)

