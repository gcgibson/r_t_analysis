


## Load deSolve package
library(deSolve)
## Globals
times      <- seq(0, 70, by = 1)
## Create an asymtpomatic SIR function
sir_asymptomatic <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -  S * (beta1*I+ beta2*A )
    dA <- S*(beta1*I+ beta2*A )*percent_asymptomatic - gamma*A
    dI <-   S * (beta1*I+ beta2*A ) *(1-percent_asymptomatic) - gamma * I
    dI_new <-    S * (beta1*I) 
    
    dR <- gamma * I + gamma*A
    
    return(list(c(dS,dA,dI_new, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init_asymptomatic       <- c(S = 1-1e-6,A = 1e-6, I_new =1e-6,I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters_asmyptomatic <- c(beta1 = .95,beta2=4.01, gamma = .75,percent_asymptomatic=.5)
## Time frame

## Solve using ode (General Solver for Ordinary Differential Equations)
out_asymptomatic <- ode(y = init_asymptomatic, times = times, func = sir_asymptomatic, parms = parameters_asmyptomatic)
## change to data frame
out_asymptomatic <- as.data.frame(out_asymptomatic)
## Delete time variable
## multiply precentage to some population size
out_asymptomatic$I_new <- round(c(1e-6,diff(out_asymptomatic$I_new))*600e7)

## Show data
plot(out_asymptomatic$I_new,type='l')


sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -  S * (beta1*I+ beta2*A )
    dA <- S*(beta1*I+ beta2*A )*percent_asymptomatic - gamma*A
    dI <-   S * (beta1*I+ beta2*A ) *(1-percent_asymptomatic) - gamma * I
    dI_new <-    S * (beta1*I+ beta2*A ) 
    
    dR <- gamma * I + gamma*A
    
    return(list(c(dS,dA,dI_new, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6,A = 1e-6, I_new =1e-6,I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta1 = .95,beta2=4.01, gamma = .75,percent_asymptomatic=.5)
## Time frame
times      <- seq(0, 70, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## multiply precentage to some population size
out$I_new <- round(c(1e-6,diff(out$I_new))*600e7)

lines(out$I_new,col='red')



## Estimate under complete testing 
library(EpiEstim)
r_t_hat <-estimate_R(out$I_new,method=c("parametric_si"),
                     config = make_config(list(mean_si=1/parameters[3],
                                               std_si=sqrt(1/parameters[3]))))
plot (r_t_hat$R$`Mean(R)`,type='l')



r_t_hat <-estimate_R(out_asymptomatic$I_new,method=c("parametric_si"),
                     config = make_config(list(mean_si=1/parameters[3],
                                               std_si=sqrt(1/parameters[3]))))
lines(r_t_hat$R$`Mean(R)`,type='l',col='red')

