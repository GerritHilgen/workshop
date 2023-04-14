#### HH model for workshop AP0610/14
## modified from  Lyoubi-Idrissi, datascience+ 26/5/2021
## parameters from Giant Squid Axon (HH, 1952)

# install.packages("deSolve")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("gridExtra")

require("deSolve")
require("ggplot2")
require("dplyr")
require("gridExtra")

#### Defining model parameters and their values
## Vna, Vk, Vl are ionic reversal potentials for Na+, K+ and leak respectively
Vna =115 
Vk = -12 
Vl = 10.6 

## gna, gk, gl are maximum conductances of sodium, potassium and leak respectively
gna = 120
gk = 36
gl= 0.3

## C is the membrane capacitance
C = 1

## I is the applied current to the membrane
I = 1

##Resting membrane potential
V = -45


################################################################################
################################################################################

## DO NOT CHANGE 
## m,n,h are gating variables (ion channel opening probabilities)
m = 0.052 
h = 0.596 
n = 0.317 

parameters = c(Vna = Vna, Vk = Vk, Vl = Vl, gna = gna, gk = gk, gl= gl, C = C)
yini_2 <- c(V = V, m = m, h = h, n = n)

## Set initial state
sHH <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    ## Ion channels functions/Gating functions (don't change)
    alpha_m <- function(v) 0.1*(25-v)/(exp((25-v)/10)-1)
    beta_m  <- function(v) 4*exp(-v/18)
    alpha_h <- function(v) 0.07*exp(-v/20)
    beta_h  <- function(v) 1/(exp((30-v)/10)+1)
    alpha_n <- function(v) 0.01*(10-v)/(exp((10-v)/10)-1)
    beta_n  <- function(v) 0.125*exp(-v/80)
    # I <- 10*sin(0.5*time)
    # I <- 10*exp( -0.125*(time - 50)^2)
    ## Derivatives
    dV <- ( I - gna*h*(V-Vna)*m^3-gk*(V-Vk)*n^4-gl*(V-Vl))/C
    dm <- alpha_m(V)*(1-m)-beta_m(V)*m
    dh <- alpha_h(V)*(1-h)-beta_h(V)*h
    dn <- alpha_n(V)*(1-n)-beta_n(V)*n 
    list(c(dV, dm, dh, dn))
  })
}
################################################################################
################################################################################

# Model application
##Set integration times
time = seq(from =0, to=50, by = 0.01)

## Running the model
print(system.time(
  out_sHH <- ode(yini_2 , func = sHH , times = time, parms = parameters)))
#Print the summary 
summary(out_sHH)


## plot parameters
s_HH_df <- as.data.frame(out_sHH)
pv <- ggplot(s_HH_df, aes(time,V)) + geom_line(color = "darkblue") 
pn <- ggplot(s_HH_df, aes(time,n)) + geom_line(color = "darkblue") 
pm <- ggplot(s_HH_df, aes(time,m)) + geom_line(color = "darkblue") 
ph <- ggplot(s_HH_df, aes(time,h)) + geom_line(color = "darkblue") 
grid.arrange(pv, pn, pm, ph, ncol=2, nrow =2)

################################################################################

## Question 1: What exactly is plotted on the four panels? 

################################################################################

## Question 2: Which parameter do you need to change to elicit more than one spike? 
# Why is that affecting the spike rate?

################################################################################

## Question 3: What is limiting the maximum spike number?

################################################################################

## Optional task: HH models can be used for studying ion channel diseases 
# (channelopathies). Neuromyotonia is a clinical condition characterized 
# functionally by hyperexcitability  of peripheral  nerves  manifesting as 
# continuous  muscle  fibre  activity. Patch Clamp experiments showed that 
# the potassium current is reduced by 25-80% in neuromyotonia. 
# Use the above HH model to simulate neuromyotonia and explain the outputs. 



