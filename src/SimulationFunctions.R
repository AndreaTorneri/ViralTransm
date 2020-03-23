#####################################################################
#
# simFunctions.R  
#
#
# This script contains a collection of functions used to simulate an epidemic dynamic based on the notion 
# of effective contact process. The Infection dynamic is described by contact rates accepted/rejected according to an infectiousness measure
# based on the viral load. Control strategies that decrease infectiousness or contact rate are included in the script.
# Precisely, we simulate:
# Isolation - infected individuals in isolation are assumed to do not make contacts anymore, i.e. contact rate is set to zero
# Quarantine - infected individuals in isolation are assumed to decrease their contact rate
# Antiviral drug - infected individuals injected with an antiviral decrease their infectiousness measure
#
# In the first part there are functions used in the simulators, as the simulation of the incubation period, infectiousness measure
# and symptomatic period. In a second part we report the simulators used to simulate the epidemic outbreak 
#
# For the function reported hereunder is specified the input and the output parameters.
# 
#
# Author: Andrea Torneri
# Last version: 16/03/2020
#############################################################################################

#############################################################################################
# nCov.diagnosis: Function to simulate the time to diagnosis (Exponentially distributed)
#
# The time to diagnosis is intended as the time from symptoms onset to the diagnosis.
#
#INPUT:
# avg.dt - average time to diagnosis (in day)
#
#OUTPUT:
# t.d - time to diagnosis
#


nCov.diagnosis<-function(avg.dt){
  t.d<-rexp(1,rate = 1/avg.dt)
  return(t.d)
}


#############################################################################################
# nCov.IncubationPeriod: Function to simulate the Incubation Period
#
# The time to incubation is the time to symptom onset since infection
# This is here simulated via a Weibull distribution with shape 2. This because, together with 
# the mean value of 5.2, this represents the incubation period for COVID-19 (Zhang et al. 2020 MedRxiv).
# The function can be easily modified for other distributions
#
#
#INPUT:
# avg.inc - average incubation time (in day)
#
#OUTPUT:
#  - incubation period
#

nCov.IncubationPeriod<-function(avg.inc){
  shape<-2
  scale<-avg.inc*shape/gamma(1/shape)
  return(rweibull(1,shape = shape, scale = scale))
}


######################################################################################################
# nCov.InfMeasure: Function evaluate the infectiousness measure at a precise time point since symptoms onset
#
# The infectiousness measure is set to have a similar shape to the viral load of the specific pathogen.
# Here we want to represent the viral load for COVID-19. Therefore we set the curve to represent the few
# data available to date. (Pan et al. 2020, Lancet; Zhou et al. 2020 Nature; Zou et al. 2020 NEJM)
#
# We choose to use a gamma distribution. When more information would be available, the script can be easily adapted
# The infectiousness measure is scaled to have the same shape for different length of the infectious period.
# The curve is defined over the incubation and symptomatic period, or alternatively, over the exposed and infectious period.
#
#
# INPUT:
# t - time since infection
# lengthI - length of the symptomatic plus incubation period 
#
# OUTPUT:
# value of the infectiousness measure at the specific time since infection
# 

nCov.InfMeasure<-function(t, lengthI){
  vload.comp<-(25/lengthI)*dgamma(25/lengthI*t, shape =7 , scale =1.43 )
  return(vload.comp)
  
}

######################################################################################################
# nCov.InfMeasure.treatment: Function evaluate the infectiousness measure at a precise time point since symptoms onset (after injection)
#
# This function is used to model the decrease in the infectiousness measure due to the administration of an antiviral
# The decrease is model via a Malthus model.
# The rate of decrease is set to represent the decrease in viral load due to the administration of Remdesivir as reported in Sheahan et al. 2020 for MERS 
#
# INPUT:
# t - time since infection
# lengthI - length of the symptomatic plus incubation period 
#
# OUTPUT:
# value of the infectiousness measure at the specific time since infection (after treatment)
# 

nCov.InfMeasure.treatment<-function(t, vload.tt, t.zero){
  y<-seq(0,30,0.01)
  maxVl<-max(dgamma(y,shape = 7, scale = 1.43))
  VlAt<-maxVl*log(10^2)/log(10^7)
  r<--log(VlAt/maxVl)/4
  vload.comp<-vload.tt*exp(-r*t)
  return(vload.comp)
}

######################################################################################################
# nCov.SymptmPeriod: Function to simulate the infectious period length
#
# This is modeled to follow an exponential distribution.
#
# INPUT:
# mu.IP: mean length of the symptomatic period
#
# OUTPUT:
# continous value for the length of the symptomatic period
# 

#Function to simulate Infectious period. Assumed to be exponentially distributed
nCov.SymptmPeriod<-function(mu.IP){
  return(rexp(1,1/mu.IP))
  
  #return(rgamma(1,shape = 37.60435, scale = 0.2473118))
  
}



####################################################################################
####################################################################################
# Simulators
# Here we report the scripts used to simulate the epidemic outbreak. Scripts are different based on the control strategy 
# they represent. More precisely:
#
# - nCov.simulator.ExpHosp.QuarIso.cctracing.aftersympt : individual traced are controlled. When they show symptoms are isolated/quarantined 
# - nCov.simulator.ExpHosp.QuarIso.cctracing.beforesympt : individual traced immediately isolated/quarantined 
# - nCov.simulator.ExpHosp.QuarIsoAfter.Treatmentbefore.cctracing : individual traced are immediately injected with an antiviral. They are put in isolation/quarantine when symptoms onset
# - nCov.simulator.ExpHosp.QuarIsoTreat.cttracing.beforesympt: individual traced are immediately injected with the antiviral and put in quarantine/isolation
#
# We fully comment the first function while for the others only the parts that differs are commented
#
# We report here the list of input parameters needed in the simulator
# n : population size
# lambda: contact rate
# sigma: probability of severe symptoms
# avg.inc: mean incubation period 
# avg.dt: mean time to diagnosis
# eta: probability of succesful tracing 
# rate.quar: contact rate in quarantine
# tfailtest: days since infection at which the test gives positive result when the tested individual has been previously infected
 
# OUTPUT
# Time events: matrix with three columns indicating the time at which an infection/recover event happens (first column) that time at which this event happens (second column) and the individual id (third column)
# Status.matrix: matrix that keep track of different features of the individaul during the epidemic, e.g. status (infected/recovered/susceptible), time of infection, time of symptom onsets, time of diagnosis, time of treatment.



nCov.simulator.ExpHosp.QuarIso.cttracing.aftersympt.asymptm<-function(n, lambda, rho, sigma,avg.inc, mu.IP,R.a,R.s,avg.dt,eta, rate.quar){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 9) #matrix containing information about the state of the individuals
  col.name<-c("infected","time.of.infection","infector", "severity", "TimeSymptomOnset","TimeQuarantine", "Treatment", "TimeTreatment","Traced")
  # infected can be: 0 (susceptible), 1 (mild-simptomatic), 2 (severe-simptomatic) / severity can be 1.0 (asymptomatic), 1.1 (mild symptoms), 1.2 (severe symptoms)
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  status.matrix[,7]<-0 #no one treated
  status.matrix[,8]<-0 
  status.matrix[,9]<-0 #no one traced
  
  recovery.vector<-rep(NA,n) #vector of recovery times
  quarantine.day<-rep(Inf,n) #day at which individuals will be quarantined
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  quarantine<-rep(0,n) #vector that says who is at in quarantine at the current time. 0 means no quarantine, 1 quarantine
  current.time<-0
  
  #transmission parameter dataframe: each line is an individual, the first colum is the transmsission coeffficient and the third the length of IP (needed to re-scale Viral load curve)
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=lambda)   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first<-sample(which(status.matrix[,1]==0), 1) #initial case
  status.matrix[first,1] <- 1 
  status.matrix[first,2] <- 0
  status.matrix[first,5]<-current.time+nCov.IncubationPeriod(avg.inc =avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
  recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.IP) # the total length since infection (Exposed+IP) 
  transmission.parameters$total_infectionPeriod[first]<-recovery.vector[first]-current.time
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  
  if (runif(1)<rho){
    if (runif(1)<sigma){ #check whether are severe
      status.matrix[first,4]<-2
      transmission.parameters$q[first]<-R.s/lambda
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
      time.events[1,]<-c(current.time,1.2,first)
      status.matrix[first,6]<-quarantine.day[first]
    }else{
      status.matrix[first,4]<-1
      short<-FALSE
      transmission.parameters$q[first]<-R.s/lambda
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
      status.matrix[first,6]<-quarantine.day[first]
      time.events[1,]<-c(current.time,1.1,first)
    }
  }else{
    status.matrix[first,4]<-0
    short<-FALSE
    transmission.parameters$q[first]<-R.a/lambda
    status.matrix[first,5]<-Inf # no symptoms
    quarantine.day[first]<-Inf # no diagnosed
    status.matrix[first,6]<-quarantine.day[first]
    time.events[1,]<-c(current.time,1.0,first)
  }
  
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  proposed.individual<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  index.contact[first]<-1 #when 1 individual proposes a contact 
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  recovered<-0
  contact.list<-list()
  for (i in 1:n){
    contact.list[[i]]<-0
  }
  
  
  err<-0
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals who has to, propose a new social contact
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      #contact.time$pr.infectee[i] <-sample(setdiff(1:n,c(i,which(quarantine!=0))),1) # individual can make contact only with people that are not quarantined
      contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
      #if (temp.contact.time>quarantine.day[i]){contact.time$pr.ctc[i]<-NA} # he does not make contact contact till is quarantine
    }
    
    #Phase 2: identify the next event: possible infection, someone is quarantined or a recovery
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf) # among all the proposed social contact between houeholds we select the minimum
    ifelse(length(which(!is.infinite(quarantine.day)))>0,Q<-min(quarantine.day),Q<-Inf ) #minimum quarantine pathogen 1
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    
    
    if (Q<=min(T_NextCtc,R_a)){
      current.time<-Q
      person.quarantined<-which(Q==quarantine.day)
      quarantine[person.quarantined]<-1
      quarantine.day[person.quarantined]<-Inf
      if (status.matrix[person.quarantined,4]==2){
        contact.time$pr.ctc[person.quarantined]<-NA #if severe isolated
      }else{
        transmission.parameters$contact_rate[person.quarantined]<-rate.quar
      }
      #Contact-tracing: with a certain probability, select a previously contacted individual and treat him/her 
      if (status.matrix[person.quarantined,9]==0){
        contacted.individuals<-NULL
        ifelse(length(contact.list[[person.quarantined]])>0, contacted.individuals<-contact.list[[person.quarantined]][-1], contacted.individuals<-NULL)
        contacted.individuals<-setdiff(contacted.individuals,which(quarantine==1))
        if (length(contacted.individuals)>0){
          for (t in unique(contacted.individuals)) {
            #print(t)
            if (rbinom(1,1,eta)==1){ #succesful tracing 
              if (status.matrix[t,1]==1){
                quarantine.day[t]<-status.matrix[t,5]
                status.matrix[t,6]<-quarantine.day[t]
              }else{
                if (status.matrix[t,9]==0 | (status.matrix[t,9]==1 & (status.matrix[t,8] < current.time ))){
                  status.matrix[t,8]<-current.time+14
                }
              }
              status.matrix[t,9]<-1
            }
          }
        }
      }
    }else{
      if (T_NextCtc<R_a){
        current.time<-T_NextCtc
        infector<-which(contact.time$pr.ctc ==T_NextCtc) 
        infectee<-sample(setdiff(1:n,infector),1)
        acc.rate<-nCov.InfMeasure(t=(current.time-status.matrix[infector,2]), lengthI = transmission.parameters$total_infectionPeriod[infector])*transmission.parameters$q[infector]
        if (acc.rate>1){err<-err+1}
        contact.list[[infector]]<-c(contact.list[[infector]],infectee)
        if (status.matrix[infectee,1]==0 & runif(1)<acc.rate){
          status.matrix[infectee,1]<-1
          status.matrix[infectee,2]<-current.time
          status.matrix[infectee,3]<-infector
          status.matrix[infectee,5]<-current.time+nCov.IncubationPeriod(avg.inc=avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
          recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.IP) # the total length since infection (Exposed+IP) 
          transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
          infectives[infectee]<-1
          if (status.matrix[infectee,9]==1){
            if ((status.matrix[infectee,5])<status.matrix[infectee,8]){
              quarantine.day[infectee]<-status.matrix[infectee,5]
            }else{
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt) 
              status.matrix[infectee,9]<-0
            }
          }else{
            quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt) 
          }
          
          if (runif(1)<rho){
            if (runif(1)<sigma){ # whether are severe
              status.matrix[infectee,4]<-2
              transmission.parameters$q[infectee]<-R.s/lambda
              status.matrix[infectee,6]<-quarantine.day[infectee]
              time.events<-rbind(time.events,c(current.time,1.2,infectee))
            }else{
              status.matrix[infectee,4]<-1
              status.matrix[infectee,6]<-quarantine.day[infectee]
              transmission.parameters$q[infectee]<-R.s/lambda
              time.events<-rbind(time.events,c(current.time,1.1,infectee))
            }
          }else{
            status.matrix[infectee,4]<-1
            transmission.parameters$q[infectee]<-R.a/lambda
            quarantine.day[infectee]<-Inf
            status.matrix[infectee,6]<-quarantine.day[infectee]
            status.matrix[infectee,5]<-Inf
            time.events<-rbind(time.events,c(current.time,1.0,infectee))
          }
          index.contact[infectee]<-1
        }
        index.contact[infector]<-1
        contact.time$pr.ctc[infector]<-NA
        #Phase 2.3 a recovery occurs
      }else{
        current.time<-R_a
        recovered<-which(recovery.vector==R_a)
        recovery.vector[recovered]<-NA
        status.matrix[recovered,1]<--1
        time.events<-rbind(time.events,c(current.time,0.1,recovered))
        contact.time[recovered,2:3]<-rep(NA,4)
        infectives[recovered]<-NA
        quarantine.day[recovered]<-Inf
        quarantine[recovered]<-0
      }
    }
  }
  #When also the other pathogen is present.
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix,err=err))
}


nCov.simulator.ExpHosp.QuarIso.aftersympt<-function(n, lambda, sigma,avg.inc, mu.IP,R,avg.dt, eta, rate.quar){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 9) #matrix containing information about the state and other features of individuals
  col.name<-c("infected","time.of.infection","infector", "severity", "TimeSymptomOnset","TimeQuarantine", "Treatment", "TimeTreatment","Traced")
  # infected can be: 0 (susceptible), 1 (mild-simptomatic), 2 (severe-simptomatic), -1 (removed) / severity can be 1.1 (mild symptoms), 1.2 (severe symptoms)
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  status.matrix[,7]<-0 #no one treated
  status.matrix[,9]<-0 #no one traced
  status.matrix[,8]<-0 

  
  recovery.vector<-rep(NA,n) #vector of recovery times
  quarantine.day<-rep(Inf,n) #day at which individuals will be diagnosed
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  quarantine<-rep(0,n) #vector that says who is at in quarantine at the current time. 0 means no quarantine, 1 quarantine
  current.time<-0
  
  #transmission parameter dataframe: each line is an individual, the first colum is the transmsission coeffficient q (see "TheoreticalDescription.pdf") and the third the length of incubation+symptomatic period (needed to re-scale Viral load curve)
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=lambda)  # the last column is the contact rate of the specific individual. When quarantined or isolated this changes.
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first<-sample(which(status.matrix[,1]==0), 1) #initial case
  status.matrix[first,1] <- 1 
  status.matrix[first,2] <- 0
  status.matrix[first,5]<-current.time+nCov.IncubationPeriod(avg.inc =avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
  recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.IP) # the total length since infection (Exposed+IP) 
  transmission.parameters$total_infectionPeriod[first]<-recovery.vector[first]-current.time
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
    if (runif(1)<sigma){ #check whether is severe
      status.matrix[first,4]<-2
      transmission.parameters$q[first]<-R/lambda
      short<-TRUE
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
      time.events[1,]<-c(current.time,1.2,first)
      status.matrix[first,6]<-quarantine.day[first]
    }else{
      status.matrix[first,4]<-1
      short<-FALSE
      transmission.parameters$q[first]<-R/lambda
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
      status.matrix[first,6]<-quarantine.day[first]
      time.events[1,]<-c(current.time,1.1,first)
    }
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) 
  #initializing variables
  proposed.individual<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new social contact - 1 yes 0 no
  index.contact[first]<-1 # the index cases starts with making contacts 
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  recovered<-0
  contact.list<-list() # list that keeps track of all the contacts an infected makes. This is used in the contact tracing
  # when an individual is diagnosed we check all the contacts he/she had
  for (i in 1:n){
    contact.list[[i]]<-0
  }
  
  err<-0
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals who has to, propose a new contact
    for (i in which(index.contact==1) ){ 
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next contact for individual i
      index.contact[i]<-0 # until individuals i resolves this contact we do not simulate other contact time for he/she
      contact.time$pr.ctc[i]<-temp.contact.time # We save the contact
    }
    
    #Phase 2: identify the next event: possible infection, someone is quarantined or a recovery
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf) # among all the proposed contact we select the minimum
    ifelse(length(which(!is.infinite(quarantine.day)))>0,Q<-min(quarantine.day),Q<-Inf ) # minimum of quarantine/isolation time
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    
    #Phase 3: resolve the next event
    if (Q<=min(T_NextCtc,R_a)){ #the next event is a quarantine event
      current.time<-Q # update the current time
      person.quarantined<-which(Q==quarantine.day) #identify the person that goes in quarantine
      quarantine[person.quarantined]<-1 # the selected person is in quarantine
      quarantine.day[person.quarantined]<-Inf # the selected person is not going to quarantine anymore. We assume stays in quarantine till the end of the infectious period.
      if (status.matrix[person.quarantined,4]==2){
        contact.time$pr.ctc[person.quarantined]<-NA # if severe isolated. To do so, we deleted the next contact. The vector index contact is at value 0, so this individual will not make other contacts
      }else{
        transmission.parameters$contact_rate[person.quarantined]<-rate.quar # we change the contact rate of this individual
      }
      #Contact-tracing: with a certain probability, select a previously contacted individual and treat him/her 
      if (status.matrix[person.quarantined,9]==0){ # if the quarantine individual has not been traced yet. (we assume only 1-level of tracing)
        contacted.individuals<-NULL
        ifelse(length(contact.list[[person.quarantined]])>0, contacted.individuals<-contact.list[[person.quarantined]][-1], contacted.individuals<-NULL)
        contacted.individuals<-setdiff(contacted.individuals,which(quarantine==1)) # We select the individuals that have been in contact with the diagnosed, considering the one that are not already in quarantine  
        if (length(contacted.individuals)>0){
          for (t in unique(contacted.individuals)) {
            #print(t)
            if (rbinom(1,1,eta)==1){ #succesful tracing 
              if (status.matrix[t,1]==1){
                quarantine.day[t]<-status.matrix[t,5]
                status.matrix[t,6]<-quarantine.day[t]
              }else{
                if (status.matrix[t,9]==0 | (status.matrix[t,9]==1 & (status.matrix[t,8] < current.time ))){
                  status.matrix[t,8]<-current.time+14
                }
              }
              status.matrix[t,9]<-1
            }
          }
        }
      }
    }else{
      if (T_NextCtc<R_a){ # the next event is a contact event
        current.time<-T_NextCtc
        infector<-which(contact.time$pr.ctc ==T_NextCtc) 
        infectee<-sample(setdiff(1:n,infector),1)
        acc.rate<-nCov.InfMeasure(t=(current.time-status.matrix[infector,2]), lengthI = transmission.parameters$total_infectionPeriod[infector])*transmission.parameters$q[infector] #acceptance rate
        if (acc.rate>1){err<-err+1}
        contact.list[[infector]]<-c(contact.list[[infector]],infectee)
        if (status.matrix[infectee,1]==0 & runif(1)<acc.rate){ # if the contact is accepted (infectee susceptible + acc. rate) -> infection event
          status.matrix[infectee,1]<-1 #update the state of the infectee
          status.matrix[infectee,2]<-current.time
          status.matrix[infectee,3]<-infector
          status.matrix[infectee,5]<-current.time+nCov.IncubationPeriod(avg.inc=avg.inc) 
          recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.IP)  
          transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
          infectives[infectee]<-1 
           if (status.matrix[infectee,9]==1){ #if the individual that has been infected is under tracing we monitor him/her. If show symptoms while being traced; he/she is
             #going to be diagnosed when symptoms onset
            if ((status.matrix[infectee,5])<status.matrix[infectee,8]){
              quarantine.day[infectee]<-status.matrix[infectee,5]
            }else{
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt)
              status.matrix[infectee,9]<-0
            }
          }else{
            quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt)
          } 
          if (runif(1)<sigma){ #check whether are severe
              status.matrix[infectee,4]<-2
              transmission.parameters$q[infectee]<-R/lambda
              status.matrix[infectee,6]<-quarantine.day[infectee]
              time.events<-rbind(time.events,c(current.time,1.2,infectee))
          }else{
              status.matrix[infectee,4]<-1
              status.matrix[infectee,6]<-quarantine.day[infectee]
              transmission.parameters$q[infectee]<-R/lambda
              time.events<-rbind(time.events,c(current.time,1.1,infectee))
            }
          index.contact[infectee]<-1
        }
        index.contact[infector]<-1 
        contact.time$pr.ctc[infector]<-NA
      }else{ #next event is a recovery
        current.time<-R_a
        recovered<-which(recovery.vector==R_a)
        recovery.vector[recovered]<-NA
        status.matrix[recovered,1]<--1
        time.events<-rbind(time.events,c(current.time,0.1,recovered))
        contact.time[recovered,2:3]<-rep(NA,4)
        infectives[recovered]<-NA
        quarantine.day[recovered]<-Inf
        quarantine[recovered]<-0
      }
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix,err=err))
}

nCov.simulator.ExpHosp.QuarIso.beforesympt.delaytesting<-function(n, lambda, rho, sigma,avg.inc, mu.IP,R,avg.dt, eta, rate.quar, tfailtest, tdelaytest){
  status.matrix <- matrix(NA,nrow = n,ncol = 9) 
  col.name<-c("infected","time.of.infection","infector", "severity", "TimeSymptomOnset","TimeQuarantine", "Treatment", "TimeTreatment","Traced")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  status.matrix[,7]<-0 #no one treated
  status.matrix[,9]<-0 #no one traced
  
  recovery.vector<-rep(NA,n) #vector of recovery times
  quarantine.day<-rep(Inf,n) #day at which individuals will be diagnosed
  delay.testing<-rep(Inf,n)
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  quarantine<-rep(0,n) #vector that says who is at in quarantine at the current time. 0 means no quarantine/isolation, 1 quarantine/isolation
  current.time<-0
  
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=lambda)
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first<-sample(which(status.matrix[,1]==0), 1) #initial case
  status.matrix[first,1] <- 1 
  status.matrix[first,2] <- 0
  status.matrix[first,5]<-current.time+nCov.IncubationPeriod(avg.inc =avg.inc) 
  recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.IP)  
  transmission.parameters$total_infectionPeriod[first]<-recovery.vector[first]-current.time
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  if (runif(1)<sigma){ #check whether are severe
    status.matrix[first,4]<-2
    transmission.parameters$q[first]<-R/lambda
    short<-TRUE
    quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
    time.events[1,]<-c(current.time,1.2,first)
    status.matrix[first,6]<-quarantine.day[first]
  }else{
    status.matrix[first,4]<-1
    short<-FALSE
    transmission.parameters$q[first]<-R/lambda
    quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
    status.matrix[first,6]<-quarantine.day[first]
    time.events[1,]<-c(current.time,1.1,first)
  }
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n))   
  proposed.individual<-0
  index.contact<-rep(0,n) 
  index.contact[first]<-1 
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  recovered<-0
  t.testing<-Inf # we keep track of the time at which individuals are possibly re-tested
  contact.list<-list()
  for (i in 1:n){
    contact.list[[i]]<-0
  }
  err<-0
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    for (i in which(index.contact==1) ){ 
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time
      index.contact[i]<-0
      contact.time$pr.ctc[i]<-temp.contact.time 
    }
    
    #Phase 1: identify the next event: possible infection, someone is quarantined or a recovery
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf) 
    ifelse(length(which(!is.infinite(quarantine.day)))>0,Q<-min(quarantine.day),Q<-Inf ) 
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    ifelse(length(which(!is.infinite(delay.testing)))>0,t.testing<-min(delay.testing),t.testing<-Inf ) 
    
    
    if (t.testing<min(Q,T_NextCtc,R_a)){ #time of re-testing: at this time infectives result positive (for the assumption of the test)
      current.time<-t.testing
      tested.individual<-which(delay.testing==t.testing)
      for (t.temp in tested.individual){
        quarantine[t.temp]<-1 #individuals are put in quarantine/isolation based on the symptoms
        if (status.matrix[t.temp,4]==2){
          contact.time$pr.ctc[t.temp]<-NA 
        }else{
          transmission.parameters$contact_rate[t.temp]<-rate.quar
        }
        quarantine.day[t.temp]<-Inf
        delay.testing[t.temp]<-Inf
      }
    }else{
      if (Q<=min(T_NextCtc,R_a)){
        current.time<-Q
        person.quarantined<-which(Q==quarantine.day)
        quarantine[person.quarantined]<-1
        quarantine.day[person.quarantined]<-Inf
        if (status.matrix[person.quarantined,4]==2){
          contact.time$pr.ctc[person.quarantined]<-NA #if severe isolated
        }else{
          transmission.parameters$contact_rate[person.quarantined]<-rate.quar
        }
        if (status.matrix[person.quarantined,9]==0){ 
          contacted.individuals<-NULL
          ifelse(length(contact.list[[person.quarantined]])>0, contacted.individuals<-contact.list[[person.quarantined]][-1], contacted.individuals<-NULL)
          contacted.individuals<-setdiff(contacted.individuals,which(quarantine==1))
          if (length(contacted.individuals)>0){
            for (t in unique(contacted.individuals)){
              #print(t)
              if (rbinom(1,1,eta)==1){ 
                if (status.matrix[t,1]==1){
                  if ((current.time -status.matrix[t,2])<(tfailtest)){ #infectives that failed the test 
                    delay.testing[t]<-min((current.time+tfailtest),status.matrix[t,6]) # are re-tested, if not showing symptoms before
                  }else{
                    quarantine[t]<-1 #if positive quarantine/isolation
                    quarantine.day[t]<-Inf
                    if (status.matrix[t,4]==2){
                      contact.time$pr.ctc[t]<-NA #if severe isolated
                    }else{
                      transmission.parameters$contact_rate[t]<-rate.quar
                    }
                  }
                  status.matrix[t,9]<-1
                }
              }
            }
          }
        }
      }else{
        if (T_NextCtc<R_a){
          current.time<-T_NextCtc
          infector<-which(contact.time$pr.ctc ==T_NextCtc) 
          infectee<-sample(setdiff(1:n,infector),1)
          acc.rate<-nCov.InfMeasure(t=(current.time-status.matrix[infector,2]), lengthI = transmission.parameters$total_infectionPeriod[infector])*transmission.parameters$q[infector]
          if (acc.rate>1){err<-err+1}
          contact.list[[infector]]<-c(contact.list[[infector]],infectee)
          if (status.matrix[infectee,1]==0 & runif(1)<acc.rate){
            status.matrix[infectee,1]<-1
            status.matrix[infectee,2]<-current.time
            status.matrix[infectee,3]<-infector
            status.matrix[infectee,5]<-current.time+nCov.IncubationPeriod(avg.inc=avg.inc) 
            recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.IP)  
            transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
            infectives[infectee]<-1
            if (runif(1)<sigma){ 
              status.matrix[infectee,4]<-2
              transmission.parameters$q[infectee]<-R/lambda
              short<-TRUE
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt) 
              status.matrix[infectee,6]<-quarantine.day[infectee]
              time.events<-rbind(time.events,c(current.time,1.2,infectee))
            }else{
              status.matrix[infectee,4]<-1
              short<-FALSE
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt) 
              status.matrix[infectee,6]<-quarantine.day[infectee]
              transmission.parameters$q[infectee]<-R/lambda
              time.events<-rbind(time.events,c(current.time,1.1,infectee))
            }
            index.contact[infectee]<-1
          }
          index.contact[infector]<-1
          contact.time$pr.ctc[infector]<-NA
        }else{
          current.time<-R_a
          recovered<-which(recovery.vector==R_a)
          recovery.vector[recovered]<-NA
          status.matrix[recovered,1]<--1
          time.events<-rbind(time.events,c(current.time,0.1,recovered))
          contact.time[recovered,2:3]<-rep(NA,4)
          infectives[recovered]<-NA
          quarantine.day[recovered]<-Inf
          quarantine[recovered]<-0
        }
      }
    }
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix,err=err))
}



nCov.simulator.ExpHosp.QuarIsoTreat.beforesympt.delaytesting<-function(n, lambda, rho, sigma,avg.inc, mu.IP,R,avg.dt, eta, rate.quar, tfailtest, tdelaytest){
  status.matrix <- matrix(NA,nrow = n,ncol = 9) #matrix containing information about the state of the individuals
  col.name<-c("infected","time.of.infection","infector", "severity", "TimeSymptomOnset","TimeQuarantine", "Treatment", "TimeTreatment","Traced")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  status.matrix[,7]<-0 #no one treated
  status.matrix[,9]<-0 #no one traced
  
  recovery.vector<-rep(NA,n) #vector of recovery times
  quarantine.day<-rep(Inf,n) #day at which individuals will be diagnosed
  delay.testing<-rep(Inf,n)
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  quarantine<-rep(0,n) #vector that says who is in quarantine/isolation at the current time. 0 means no quarantine/isolation, 1 quarantine/isolation
  current.time<-0
  
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=lambda)   
  
  # first infected: randomly chosen in the population (among the susceptibles)
  first<-sample(which(status.matrix[,1]==0), 1) #initial case
  status.matrix[first,1] <- 1 
  status.matrix[first,2] <- 0
  status.matrix[first,5]<-current.time+nCov.IncubationPeriod(avg.inc =avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
  recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.IP) # the total length since infection (Exposed+IP) 
  transmission.parameters$total_infectionPeriod[first]<-recovery.vector[first]-current.time
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  if (runif(1)<sigma){ #check whether are severe
    status.matrix[first,4]<-2
    transmission.parameters$q[first]<-R/lambda
    short<-TRUE
    quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
    time.events[1,]<-c(current.time,1.2,first)
    status.matrix[first,6]<-quarantine.day[first]
  }else{
    status.matrix[first,4]<-1
    short<-FALSE
    transmission.parameters$q[first]<-R/lambda
    quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt)
    status.matrix[first,6]<-quarantine.day[first]
    time.events[1,]<-c(current.time,1.1,first)
  }
  
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n))  
  proposed.individual<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new contact - 1 yes 0 no
  index.contact[first]<-1 #when 1 individual proposes a contact 
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  recovered<-0
  t.testing<-Inf
  contact.list<-list()
  for (i in 1:n){
    contact.list[[i]]<-0
  }
  err<-0
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      contact.time$pr.ctc[i]<-temp.contact.time
    }
    
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf)
    ifelse(length(which(!is.infinite(quarantine.day)))>0,Q<-min(quarantine.day),Q<-Inf ) 
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    ifelse(length(which(!is.infinite(delay.testing)))>0,t.testing<-min(delay.testing),t.testing<-Inf ) 
    
    
    if (t.testing<min(Q,T_NextCtc,R_a)){ #Time at re-test. Isolation/quarantine and treatment are both performed
      current.time<-t.testing
      tested.individual<-which(delay.testing==t.testing)
      quarantine[tested.individual]<-1
      quarantine.day[tested.individual]<-Inf
      for(t.temp in tested.individual){
        if (status.matrix[t.temp,4]==2){
          contact.time$pr.ctc[t.temp]<-NA #if severe isolated
        }else{
          transmission.parameters$contact_rate[t.temp]<-rate.quar
          if (status.matrix[t.temp,7]==0){
            status.matrix[t.temp,7]<-1
            status.matrix[t.temp,8]<-current.time
          }
        }
        delay.testing[t.temp]<-Inf
      }
    }else{
      if (Q<=min(T_NextCtc,R_a)){ # at time of diagnosis isolation/quarantine and treatment are performed
        current.time<-Q
        person.quarantined<-which(Q==quarantine.day)
        quarantine[person.quarantined]<-1
        quarantine.day[person.quarantined]<-Inf
        if (status.matrix[person.quarantined,4]==2){
          contact.time$pr.ctc[person.quarantined]<-NA #if severe isolated
        }else{
          transmission.parameters$contact_rate[person.quarantined]<-rate.quar
          if (status.matrix[person.quarantined,7]==0){
            status.matrix[person.quarantined,7]<-1
            status.matrix[person.quarantined,8]<-current.time
          }
        }
        if (status.matrix[person.quarantined,9]==0){
          contacted.individuals<-NULL
          ifelse(length(contact.list[[person.quarantined]])>0, contacted.individuals<-contact.list[[person.quarantined]][-1], contacted.individuals<-NULL)
          contacted.individuals<-setdiff(contacted.individuals,which(quarantine==1))
          if (length(contacted.individuals)>0){
            for (t in unique(contacted.individuals)){
              #print(t)
              if (rbinom(1,1,eta)==1){ #succesful tracing 
                if (status.matrix[t,1]==1){
                  if ((current.time -status.matrix[t,2])<(tfailtest)){
                    delay.testing[t]<-min((current.time+tdelaytest),status.matrix[t,6])
                  }else{ #positevely tested individuals are isolated/quarantined and treated
                    quarantine[t]<-1
                    quarantine.day[t]<-Inf
                    if(status.matrix[t,4]==2){
                      contact.time$pr.ctc[t]<-NA #if severe isolated
                    }else{
                      transmission.parameters$contact_rate[t]<-rate.quar
                      if (status.matrix[t,7]==0){
                        status.matrix[t,7]<-1
                        status.matrix[t,8]<-current.time
                      }
                    }
                  }
                  status.matrix[t,9]<-1
                }
              }
            }
          }
        }
      }else{
        if (T_NextCtc<R_a){
          current.time<-T_NextCtc
          infector<-which(contact.time$pr.ctc ==T_NextCtc) 
          infectee<-sample(setdiff(1:n,infector),1)
          if (status.matrix[infector,7]==1){
            acc.rate<-nCov.InfMeasure.treatment(t=(current.time-status.matrix[infector,8]),vload.tt = nCov.InfMeasure(status.matrix[infector,8],lengthI =transmission.parameters$total_infectionPeriod[infector] ), t.zero = t.zero)*transmission.parameters$q[infector]
          }else{
            acc.rate<-nCov.InfMeasure(t=(current.time-status.matrix[infector,2]), lengthI = transmission.parameters$total_infectionPeriod[infector])*transmission.parameters$q[infector]
          }
          if (acc.rate>1){err<-err+1}
          contact.list[[infector]]<-c(contact.list[[infector]],infectee)
          if (status.matrix[infectee,1]==0 & runif(1)<acc.rate){
            status.matrix[infectee,1]<-1
            status.matrix[infectee,2]<-current.time
            status.matrix[infectee,3]<-infector
            status.matrix[infectee,5]<-current.time+nCov.IncubationPeriod(avg.inc=avg.inc) 
            recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.IP)  
            transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
            infectives[infectee]<-1
            if (runif(1)<sigma){ #check whether are severe
              status.matrix[infectee,4]<-2
              transmission.parameters$q[infectee]<-R/lambda
              short<-TRUE
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt) 
              status.matrix[infectee,6]<-quarantine.day[infectee]
              time.events<-rbind(time.events,c(current.time,1.2,infectee))
            }else{
              status.matrix[infectee,4]<-1
              short<-FALSE
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt)              
              status.matrix[infectee,6]<-quarantine.day[infectee]
              transmission.parameters$q[infectee]<-R/lambda
              time.events<-rbind(time.events,c(current.time,1.1,infectee))
            }
            index.contact[infectee]<-1
          }
          index.contact[infector]<-1
          contact.time$pr.ctc[infector]<-NA
        }else{
          current.time<-R_a
          recovered<-which(recovery.vector==R_a)
          recovery.vector[recovered]<-NA
          status.matrix[recovered,1]<--1
          time.events<-rbind(time.events,c(current.time,0.1,recovered))
          contact.time[recovered,2:3]<-rep(NA,4)
          infectives[recovered]<-NA
          quarantine.day[recovered]<-Inf
          quarantine[recovered]<-0
        }
      }
      
    }
    
    
  }
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix,err=err))
}


