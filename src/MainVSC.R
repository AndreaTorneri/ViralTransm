args <- commandArgs(trailingOnly = TRUE)

wd = args[1]
cores = args[2]

setwd(wd)

#load packages
library(foreach)
library(doParallel)
# load functions
source("SimulationFunctions.R")
# Parameters
n<-500 #Population size
lambda<-12 #daily rate of social contacts
pbq.s<-1 #probability that a (mild) symptomatic individual is quarantined
pbq.ss<-1 #probability that a severe symptomatic individual is quarantined
rho<- 0.6 # proportion of symptomatic individual
sigma<- 0.16 # proportion of severe cases (computed the 15.2.2020 from: https://docs.google.com/spreadsheets/d/1Z7VQ5xlf3BaTx_LBBblsW4hLoGYWnZyog3jqsS9Dbgc/edit#gid=957283529)
avg.inc<-5.2 # NEJM Li
mu.IP<-18 #28.6 #2.3 NEJM Li (SerialInt-Inc) -Anne Cori - SARS
R.s<-2.5 # reproduction number symptomatic (mild)
R.ss<-2.5 # reproduction number severe symptomatic
R.a<-R.s #reproduction number asymptomatic
avg.dt<-1.6
t.zero<-6
eta<-0.5
tfailtest<-2
tdelaytest<-2
lambda.dec<-0.25*lambda
VShed<-"After"

cl<-makeCluster(as.numeric(cores))
registerDoParallel(cl)

nSim<-4
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)

epi.outbreak<-foreach(i = 1:nSim) %dopar%{
  #epi.outbreak[[i]]<-nCov.simulator.ExpQuar (n=n, lambda = lambda, pbq.s = pbq.s, pbq.ss = pbq.ss, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s, R.ss=R.ss,tqs=tqs,tqss = tqss, rq=rq)
  #epi.outbreak[[i]]<-nCov.simulator.GamIso(n=n, lambda = lambda, pbq.s = pbq.s, pbq.ss = pbq.ss, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s, R.ss=R.ss,tqs=tqs,tqss = tqss)
  #epi.outbreak[[i]]<-nCov.simulator.ConstQuar(n=n, lambda = lambda, pbq.s = pbq.s, pbq.ss = pbq.ss, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s, R.ss=R.ss,tqs=tqs,tqss = tqss)
  epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.cttracing.aftersympt.asymptm (n=n, lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s,avg.dt=avg.dt, rate.quar = lambda.dec, eta=eta)  
  #epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.cttracing.beforesympt.delaytesting.asymptm(n=n, lambda = lambda, pbq.s = pbq.s, pbq.ss = pbq.ss, rho = rho, sigma = sigma, avg.inc = avg.inc, mu.IP = mu.IP, R.a = R.a, R.s=R.s, R.ss=R.ss,tqs=tqs,tqss = tqss, rate.quar = lambda.dec, eta=eta, t.zero = t.zero, tfailtest = tfailtest, tdelaytest = tdelaytest)  
  
  #print(i)
}
stopCluster(cl)


name<-paste("Asymptomatic",n,"_propUnd",rho,"_Rs",R.s,"_Ra",R.a,"_lambda",lambda,"_Scenario",VShed,"_eta",eta,"_lambda_dec",lambda.dec,"_tfail",tfailtest,"_.RData", sep = "")

save(epi.outbreak, file = name)

