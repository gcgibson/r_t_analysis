
run_analysis_at_percentage <- function(percentage,beta1,beta2){
      times      <- seq(0, 70, by = 1)
      print (beta1)
      parameters <- c(beta1 = beta1,beta2=beta2, gamma1 = .75,gamma2=.75)
      percent_asymptomatic  =  percentage
      init       <- c(S = 1-1e-3 -  percent_asymptomatic*1e-3,A = percent_asymptomatic*1e-3, I_new =1e-3,I = 1e-3, R = 0.0)
      
      ## Create an asymtpomatic SIR function
      sir_observe_symptoms <- function(time, state, parameters) {
        
        with(as.list(c(state, parameters)), {
          
          dS <- -  S * (beta1*I+ beta2*A )
          dA <- S*beta2*A  - gamma1*A
          dI <-   S * beta1*I - gamma2 * I
          dI_new <-    S * (beta1*I) 
          
          dR <- gamma1 * I + gamma2*A
          
          return(list(c(dS,dA,dI_new, dI, dR)))
        })
      }
      
      ### Set parameters
      ## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
      ## beta: infection parameter; gamma: recovery parameter
      ## Time frame
      
      ## Solve using ode (General Solver for Ordinary Differential Equations)
      out_observe_symptoms <- ode(y = init, times = times, func = sir_observe_symptoms, parms = parameters)
      ## change to data frame
      out_observe_symptoms <- as.data.frame(out_observe_symptoms)
      ## Delete time variable
      ## multiply precentage to some population size
      out_observe_symptoms$I_new <- round(c(1e-6,diff(out_observe_symptoms$I_new))*1e8)
      
      ## Show data
      plot(out_observe_symptoms$I_new,type='l')
      
      
      ################################
      ## Create an asymtpomatic SIR function
      sir_asymptomatic <- function(time, state, parameters) {
        
        with(as.list(c(state, parameters)), {
          
          dS <- -  S * (beta1*I+ beta2*A )
          dA <- S*beta2*A  - gamma1*A
          dI <-   S * beta1*I - gamma2 * I
          dI_new <-    S * (beta1*I + beta2*A) 
          
          dR <- gamma1 * I + gamma2*A
          
          return(list(c(dS,dA,dI_new, dI, dR)))
        })
      }
      
      ### Set parameters
      
      ## Solve using ode (General Solver for Ordinary Differential Equations)
      out_asymptomatic <- ode(y = init, times = times, func = sir_asymptomatic, parms = parameters)
      ## change to data frame
      out_asymptomatic <- as.data.frame(out_asymptomatic)
      ## Delete time variable
      ## multiply precentage to some population size
      out_asymptomatic$I_new <- round(c(1e-6,diff(out_asymptomatic$I_new))*1e8)
      plot(out_asymptomatic$I_new,type='l')
      ## Show data
      
      r_t_hat <-estimate_R(out_asymptomatic$I_new,method=c("parametric_si"),
                           config = make_config(list(mean_si=1/parameters_asmyptomatic[3],
                                                     std_si=sqrt(1/parameters_asmyptomatic[3]))))

      r_t_hat_symptomatic <-estimate_R(out_observe_symptoms$I_new,method=c("parametric_si"),
                           config = make_config(list(mean_si=1/parameters_asmyptomatic[3],
                                                     std_si=sqrt(1/parameters_asmyptomatic[3]))))
      
      data_for_ggplot <- data.frame(t=seq(1,length((r_t_hat$R$`Mean(R)`))),
                                          r_t_symptomatic_only=r_t_hat_symptomatic$R$`Mean(R)`,
                                        r_t_true = r_t_hat$R$`Mean(R)`)
    
      p <- ggplot(data_for_ggplot,aes(x=t,y=r_t_true,col='True')) + geom_line() +
           geom_line(aes(x=t,y=r_t_symptomatic_only,col='Only Observing Symptomatic')) + ylab("R_t") + xlab("t") + theme_bw()
      return (p)
      

}

library(shiny)

ui <- fluidPage(
  
  # App title ----
  titlePanel("Exploration of effect of asymptomatic cases on R_t estimation"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      numericInput(inputId = "pa",
                   label = "Percentage Asymptomatic",
                   value = 0),
      
      numericInput(inputId = "beta1",
                   label = "Symptomatic Transmission",
                   value = 1),
      numericInput(inputId = "beta2",
                   label = "Asymptomatic Transmission",
                   value = 1)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
   
    
    
    p <- run_analysis_at_percentage(input$pa/100,input$beta1,input$beta2)
    print (p)
    
  })
  
}
shinyApp(ui, server)


