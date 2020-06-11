epi <- read.csv("/Users/gcgibson/tmp.csv",header = T)
# cast variables to correct types 
epi$infected <- ifelse(epi$infected == "true",T,F)
epi$susceptible <- ifelse(epi$susceptible == "true",T,F)
epi$cured <- ifelse(epi$cured == "true",T,F)

# create serial interval distribution based on data 

## get a list of unique agents
unique_agents <- unique(epi$who)
# define empty serial interval samples
serial_interval_samples <- c()
# iterate through agent list to 
# 1) Ask if they are infected
# 2) if they are infected, who infected them
# 3) When was infector infected 
# 4) When was agent infected
# 5) Compute difference in times
for (agent in unique_agents){
  # ask if agent is infected
  is_infected <- any(epi[epi$who == agent,]$infected)
  # if infected 
  if (is_infected){
    # get infector
    infector <- epi[epi$who == agent & epi$infected == T,]$infector[1]
    # remove intital infected people
    if (infector != 0){
      time_infector_was_infected <- get_time_of_infection(epi[epi$who == infector,])
      time_infectee_was_infected <- get_time_of_infection(epi[epi$who == agent,])
    
      serial_interval_samples <- c(serial_interval_samples,time_infectee_was_infected-time_infector_was_infected)
    }
  }
}
hist(serial_interval_samples)


### get epi Curve
infected <- c()

unique_times <- unique(epi$ticks)

for (time in unique_times){
  current_infected <- sum(epi[epi$ticks == time,]$infected)
  infected <- c(infected,current_infected)
}
plot(infected)


## Load deSolve package
library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = .95, gamma = .15)
## Time frame
times      <- seq(0, 70, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)


r_t_hat <-estimate_R(out$I*1000,method=c("parametric_si"),
                     config = make_config(list(mean_si=1/parameters[2],
                                               std_si=sqrt(1/parameters[2]))))
plot(r_t_hat)

## Use EpiEtim

library(EpiEstim)


r_t_hat <-estimate_R(infected,method=c("parametric_si"),
           config = make_config(list(mean_si=mean(serial_interval_samples),
                                     std_si=sqrt(var(serial_interval_samples)),
                                                 t_start=2,t_end=9)))
plot(r_t_hat)


## Constant testing
