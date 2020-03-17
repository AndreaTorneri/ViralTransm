

#load packages
# load functions
source("SimulationFunctions.R")
# Parameters
n<-500 #Population size
lambda<-12 #daily rate of  contacts
rho<- 1 # proportion of symptomatic individual
sigma<- 0.16 #proportion severe cases
mu.IP<-18  #mean length incubation period 
R<-2.5 #Reproduction number
avg.dt<-1.6 #Mean time to diagnosis
avg.inc<-5.2
eta<-0.75 #probability of tracing correctly an individual
tfailtest<-2 # time since infection at which the test detects positive an infectious individual
tdelaytest<-2 #time of re-test
rate.quar<-0.1*lambda #contact rate in quarantine
scenario<-"IAS"


nSim<-5000
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)


for (i in 1:nSim){
  epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.aftersympt(n=n, lambda=lambda, sigma=sigma,avg.inc=avg.inc, mu.IP=mu.IP,R,avg.dt=avg.dt, eta=eta, rate.quar=rate.quar)
# epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIso.beforesympt.delaytesting(n=n, lambda=lambda, sigma=sigma,avg.inc=avg.inc, mu.IP=mu.IP,R=R,avg.dt=avg.dt, eta=eta, rate.quar=rate.quar, tfailtest = tfailtest, tdelaytest = tdelaytest)
 # epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoAfter.Treatment.beforesympt.delaytesting(n=n, lambda=lambda, sigma=sigma,avg.inc=avg.inc, mu.IP=mu.IP,R=R,avg.dt=avg.dt, eta=eta, rate.quar=rate.quar, tfailtest = tfailtest, tdelaytest = tdelaytest)
  #epi.outbreak[[i]]<-nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.delaytesting(n=n, lambda=lambda, sigma=sigma,avg.inc=avg.inc, mu.IP=mu.IP,R=R,avg.dt=avg.dt, eta=eta, rate.quar=rate.quar, tfailtest = tfailtest, tdelaytest = tdelaytest)
}


name<-paste("SimNCovid",n,"_R",R.s,"_lambda",lambda,"_Scenario",scenario,"_eta",eta,"_lambda_dec",lambda.dec,"_tfail",tfailtest,"_.RData", sep = "")

save(epi.outbreak, file = name)

