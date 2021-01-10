## IPM for Gizzard Shad (Spring 2021)

library(tidyverse)  #used to make pretty graphs
library(drc) # used for logit model fitting


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Gizzard shad (gs) at LaGrange Ranch(LG) before carp invasion 1983-1999

m.par <- c(
  ## survival probability assumed to be similar to asian carp (EE 2017)
  ## doubled inflection point since gizzard shad are reportedly twice the length
  surv.min  =  0.10, # min survival - 
  surv.max    =  0.70, # max survival - ???
  surv.alpha = 80, # inflection point
  surv.beta = -5, # slope
  ## growth from Catalano & Allen (2010) paper
  grow.max  =   394.30, # maximum length in mm: L_inf = 394.30
  grow.rate =   0.60, # growth rate: 
  grow.sd   =   10,  # growth sd
  ## recruit - From Bodoga (1955) estimated from Figure 14
  recruit.mean = 112,
  recruit.sd = 40,
  egg_viable = 0.002,
  #### Spawning Probability
  prob_spawn = 0.90
#  max_spawn_prob = 0.90,  ## reference?  to match AC?
#  spawn_decay = 0.05 # chosen to keep prob. of reproducing relatively constant
)

#### FOR SIMULATIONS
N <- 50 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty
delta.z <- (U.shad-L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
#####################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Estimate parameters from data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### EGG PRODUCTION
# Using the batch fecundity vs length (Figure 1a) from Jons & Miranda (1997)
egg.size.data <- read.csv('./EggSizeData.csv',header=TRUE,sep=",")
# assume zero eggs below lowest length in data
zero.eggs <- data.frame(x = seq(floor(min(egg.size.data$x))), 
                        EggLengthData=rep(0,floor(min(egg.size.data$x))))
# assume max eggs for length greater than recorded in data
max.eggs <-  data.frame(x = seq(ceiling(max(egg.size.data$x)), floor(max(zmesh))), 
                        EggLengthData=max(egg.size.data$EggLengthData))
egg.size.data.ext <- rbind(zero.eggs, egg.size.data, max.eggs)

### USE LL3. We fit a 3 parameter (max, slope, intercept) logit model to the eggs produced data 
# since biological observations suggest that eggs are not produced by females below size 140mm.
egg.ext.m3 <- drm(formula=EggLengthData~x, data=egg.size.data.ext,
                  fct = LL.3(names=c("slope", "max","intercept")))

eggs_z <- function(z) # Eggs produced (in thousands)
{
  return(egg.ext.m3$coefficients[2]/(1+exp(egg.ext.m3$coefficients[1]*(log(z)-log(egg.ext.m3$coefficients[3])))))
}

plot.df <- data.frame(z = zmesh, eggs =  eggs_z(zmesh) )
ggplot(data = plot.df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced",
       subtitle = "Gizzard Shad") + 
  scale_x_continuous(limits = c(0,U.shad), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,700), breaks = seq(0,700,100))+
  geom_point(data = egg.size.data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_classic()+
  theme(text = element_text(size=16),
        aspect.ratio = .7)

### SURVIVAL AGE0
## import data from Michaletz paper
Michaletz.data <- read.csv('./Michaletz2009.csv',header=TRUE,sep=",")
names(Michaletz.data) <- c('density', 'survival')
##setting negative values to 0

surv.den.exp <- nls(survival/100~s0*exp(-b*density), 
                    data = Michaletz.data,
                    start = list(s0=0.3, b=0.01))

surv.density <- function(d) # probability of survival of age-0 fish dependent on density d
{
  return(coef(surv.den.exp)[1]*exp(-coef(surv.den.exp)[2]*d))
}

dmesh <- seq( from = min(Michaletz.data$density), 
              to =max(Michaletz.data$density), 
              length.out = length(Michaletz.data$density) )

plot.df <- data.frame(x = dmesh , prob = surv.density(dmesh) )
ggplot(data = plot.df,
       aes( x = x, y = prob))+
  geom_line(color = "blue", size = 1)+
  labs(x = "density (age-0 per 1000 m^3)",
       y = "probability of survival",
       title = "Survival Probability for Age-0",
       subtitle = "Gizzard Shad") + 
  scale_x_continuous(limits = c(0,max(Michaletz.data$density)), 
                     breaks = seq(0,max(Michaletz.data$density),200))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.1))+
  geom_point(data = Michaletz.data, 
             aes(x = density, y = survival/100))+
  theme_classic()+
theme(text = element_text(size=16),
      aspect.ratio = .7)
#dev.off()


######################
## Section 3: Model Set-up
######################

################# Growth function
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation L(t) = L_inf(1-e^r(t-t0))
# to find L(t+1) = L_inf*(1-e^(-r)) + e^(-r)*L(t)

G_z1z <- function(z1, z, m.par)
{
  mu <- m.par["grow.max"]*(1-exp(-m.par["grow.rate"])) + exp(-m.par["grow.rate"]) * z           # mean size next year
  sig <- m.par["grow.sd"]                                    # sd about mean
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)             # pdf that you are size z1 given you were size z
  return(p.den.grow)
}

################# Survival function, sigmoidal
s_z <- function(z, m.par)
{
  m.par["surv.min"] + (m.par["surv.max"]-m.par["surv.min"])/(1+exp(m.par["surv.beta"]*(log(z)-log(m.par["surv.alpha"]))))
  }

#################Reproduction

# survival age-0
surv_age0 <- function(n,z, m.par)
{
  #distribution of VIABLE age0 from current population n
  age0_dist <- m.par["egg_viable"]*(eggs_z(z)*1000)*m.par["prob_spawn"]*n
  age0_density <- 10**(-3)*(sum(age0_dist)*delta.z)
  return(surv.density(age0_density))
  #return(m.par['sage0_int']*exp(-m.par['sage0_decay']*age0_density)/100)
}

## Recruit size pdf

c_0z1 <- function(z1, m.par)
{
  mu <- m.par["recruit.mean"]
  sig <- m.par["recruit.sd"]
  p.den.recruit <- dnorm(z1, mean = mu, sd = sig)
  p.den.recruit <- p.den.recruit/(sum(p.den.recruit)*delta.z)
  return(p.den.recruit)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Define the survival kernel
P_z1z <- function(z1, z, m.par) {
  G_matrix <- matrix(0,N,N)
  for (x in 1:N){
    G_matrix[,x] <- G_z1z(z, rep(z[x], times = N), m.par)
    G_matrix[,x] <- G_matrix[,x]/(sum(G_matrix[,x])*delta.z)
  }
#  return(s_z(z, m.par) * G_matrix)
  return(G_matrix %*% diag(s_z(z, m.par)))
  }

## Define the fecundity kernel
F_z1z <- function (z1, z,n, m.par) {
 age1_dist <- m.par["prob_spawn"]*(eggs_z(z)*1000)*m.par["egg_viable"]*surv_age0(n,z, m.par) 
 # returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z) 
  return( outer(c_0z1(z1,m.par), age1_dist))
 #  return( p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par))
}

####### Survival Function

#### Model Analysis. No comparison to data
N <- 10 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z = (U.shad - L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####

Tf <- 40 # number of years
n <- matrix(0,length(zmesh),Tf)
n0.total <- 5000
n[,1] <-dnorm(zmesh, mean = (U.shad + L.shad)/2, sd = 50)  
n[,1] <- (n[,1]/sum(n[,1]))*n0.total/delta.z

# Note: sum(n[,1])*delta.z = n0.total

# Dynamical system
for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}

# Population size vs time
n.total <- rep(0,Tf)
for (i in 1:Tf){
  n.total[i] <- sum(n[,i])*delta.z
}

#plt.df <- data.frame(time.years=seq(Tf), number = n.total)
#ggplot(data=plt.df, aes(x=time.years, y=number)) +
#  geom_path() + 
#  theme_bw()

# survival0 vs time
surv_t <- rep(0, times=(Tf-10))
for (i in 1:(Tf-10-1)){
  surv_t[i] <- surv_age0(n=n[,i], z=zmesh, m.par)
}

### FIGURE 2(a)
plt.df <- data.frame(time.years=seq(Tf-10), prob = surv_t)
ggplot(data=plt.df, 
       aes(x=time.years, y=prob)) +
  geom_line(color = "blue", size =1)+
  labs(x = "time (in years)",
       y = "probability of survival",
       title = "Survival Probability of Age-0",
       subtitle = "Gizzard Shad") + 
  scale_x_continuous(limits = c(0,Tf-10), breaks = seq(0,Tf-10,5),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,0.008), breaks = seq(0,0.008,.001),
                     expand = c(0,0))+
  theme_classic()


# lambda vs time
lambda <- rep(0,Tf)
for (i in 1:(Tf-1)){
  lambda[i] <- n.total[i+1] / n.total[i]
}

### FIGURE 2(b)
plt.df <- data.frame(time.years=seq(Tf), ratio = lambda)
ggplot(data=plt.df, 
       aes(x=time.years, y=ratio)) +
  geom_line(color = "blue", size =1)+
  labs(x = "time (in years)",
       y = "lambda",
       title = "Relative Growth Rate - Gizzard Shad") + 
  scale_x_continuous(limits = c(0,Tf-10), breaks = seq(0,Tf-10,5),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0.4,2.6), breaks = seq(0.4,2.6,0.2),
                     expand = c(0,0))+
  theme_classic()


##### Graph n over time
#### Model Analysis. No comparison to data
N <- 50 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z = (U.shad - L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####


Tf <- 40 # number of years
n <- matrix(0,length(zmesh),Tf)
n0.total <- 5000
n[,1] <-dnorm(zmesh, mean = (U.shad + L.shad)/2, sd = 50)  
n[,1] <- (n[,1]/sum(n[,1]))*n0.total/delta.z

# Note: sum(n[,1])*delta.z = n0.total

# Dynamical system
for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}

# Population size vs time
n.total <- rep(0,Tf)
for (i in 1:Tf){
  n.total[i] <- sum(n[,i])*delta.z
}

plot.df <- data.frame(z = zmesh, density = n[,1:5] )
ggplot(data = plot.df,
       aes( x = z, y = density.1, color = "t=0"))+
  geom_line(size =1)+
  geom_line(data = plot.df,
            aes( x = z, y = density.2, color="t=1"))+
  geom_line(data = plot.df,
            aes( x = z, y = density.3, color="t=2"))+
  geom_line(data = plot.df,
            aes( x = z, y = density.4, color="t=3"))+
#  geom_line(data = plot.df,
#            aes( x = z, y = density.5, color="purple"))+
  labs(color = "Year",
       x = "length (in mm)",
       y = "population density",
       title = "Population of Gizzard Shad") + 
  scale_x_continuous(limits = c(0,U.shad), breaks = seq(0,500,100), 
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10), 
                     expand = c(0,0))+
  scale_color_manual(breaks = c("t=0", "t=1", "t=2", "t=3"),
                     values = c("t=0" = "blue", "t=1" = "green", 
                                "t=2" = "orange", "t=3" = "red")) +
  theme_classic()




######################
#### LTMR Data from LG
######################

##### FOR COMPARISON WITH DATA
N <- 11 # number of size classes
# must use 11 for data comparison
L.shad <- 60   # lower size limit in mm
U.shad <- 390    # upper size limit in mm - we want this to be larger than L-infty
# As we test the model, we will need some equal spaced meshpoints. This can be moved down later.
delta.z <- (U.shad-L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####

# initial length distribution from LTMR data
lg_data <- read.csv('./LGdata.csv',header=TRUE,sep=",")
lg2002_data <- lg_data[lg_data$year==2002,]
lg2005_data <- lg_data[lg_data$year==2005,]

Tf <- 5 # number of years

## Use counts
n <- matrix(0,length(zmesh),Tf)
n[,1] <- lg2002_data$count/delta.z

plt.df <- data.frame(length=zmesh, number = n[,1])
ggplot(data=plt.df, aes(x=length, y=number)) +
  geom_bar(stat="identity") + 
  theme_bw()

# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,0],m.par)
round(K*delta.z,2)


#max(abs(Re(eigen(K*delta.z)$values)))

#### Dynamical system
##Define dynamical system

for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}

#### Graphs
## 2002 to 2005 data
## NEEDED?

#plt.df <- data.frame( LGdata = c(rep('2002',N),rep('2005',N)), 
#                      length = rep(zmesh,2),
#                      frequency = c(lg2002_data$count/sum(lg2002_data$count),lg2005_data$count/sum(lg2005_data$count)))
#ggplot(plt.df, aes(x=length, y=frequency, fill = LGdata))+
#  geom_bar(stat="identity", position = position_dodge()) + theme_bw()

###Figure 3
## 2005 predicted vs data
plt.df <- data.frame( Legend = c(rep('predicted',N),rep('2005',N)), 
                      length = rep(zmesh,2),
                      frequency = c(n[,4]/sum(n[,4]),lg2005_data$count/sum(lg2005_data$count)))
ggplot(plt.df, aes(x=length, y=frequency, fill = Legend))+
  geom_bar(stat="identity", position = position_dodge()) + 
  theme_bw()


########## LTMR Data - La Grange
ltmr.data <- read.csv('ltrm_fish_data2.csv',header=TRUE,sep=",")
lg.data <- ltmr.data[ltmr.data$site == "LaGrange",]

#### Model Analysis. No comparison to data
N <- 50 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z = (U.shad - L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####

Tf <- length(lg.data$year) # number of years
n <- matrix(0,length(zmesh),Tf)
n0.total <- lg.data$electrofish[1]
n[,1] <-dnorm(zmesh, mean = (U.shad + L.shad)/2, sd = 50)  
n[,1] <- (n[,1]/sum(n[,1]))*n0.total/delta.z

### Dynamical System
for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}
### n totals
# Population size vs time
n.total <- rep(0,Tf)
for (i in 1:length(lg.data$year)){
  n.total[i] <- sum(n[,i])*delta.z
}

plot.df <- data.frame(year=seq(first(lg.data$year),last(lg.data$year)), 
                      number = n.total)
ggplot(data = lg.data,
       aes( x = year, y = electrofish))+
  geom_line(color = "red") +
  geom_line(data = plot.df,
            aes( x = year, y = number))+
  labs(x = "time (years)",
       y = "total density",
       title = "Density vs La Grange Observations - Gizzard Shad") + 
  scale_x_continuous(limits = c(first(lg.data$year),last(lg.data$year)), 
                     breaks = seq(first(lg.data$year),last(lg.data$year),5))+
  scale_y_continuous(limits = c(0,50000), breaks = seq(0,50000,10000))+
#  geom_point(data = Michaletz.data, 
#             aes(x = density, y = survival/100))+
  theme_classic()


########## LTMR Data - Pool 26
ltmr.data <- read.csv('ltrm_fish_data2.csv',header=TRUE,sep=",")
pool26.data <- ltmr.data[ltmr.data$site == "Pool26",]

#### Model Analysis. No comparison to data
N <- 50 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z = (U.shad - L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####

Tf <- length(pool26.data$year) # number of years
n <- matrix(0,length(zmesh),Tf)
n0.total <- pool26.data$electrofish[1]
n[,1] <-dnorm(zmesh, mean = (U.shad + L.shad)/2, sd = 50)  
n[,1] <- (n[,1]/sum(n[,1]))*n0.total/delta.z

### Dynamical System
for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}
### n totals
# Population size vs time
n.total <- rep(0,Tf)
for (i in 1:length(pool26.data$year)){
  n.total[i] <- sum(n[,i])*delta.z
}

plot.df <- data.frame(year=seq(first(pool26.data$year),last(pool26.data$year)), 
                      number = n.total)
ggplot(data = pool26.data,
       aes( x = year, y = electrofish))+
  geom_line(color = "red") +
  geom_line(data = plot.df,
            aes( x = year, y = number))+
  labs(x = "time (years)",
       y = "total density",
       title = "Density vs La Grange Observations - Gizzard Shad") + 
  scale_x_continuous(limits = c(first(pool26.data$year),last(pool26.data$year)), 
                     breaks = seq(first(pool26.data$year),last(pool26.data$year),5))+
  scale_y_continuous(limits = c(0,14000), breaks = seq(0,14000,2000))+
  #  geom_point(data = Michaletz.data, 
  #             aes(x = density, y = survival/100))+
  theme_classic()





#########################################
#### Using data from Life History and Ecology of the Gizzard Shad, Dorosoma 
#### cepedianum (Le Sueur) with Reference to Elephant Butte Lake

N <- 29 # number of size classes
L.shad <- 70   # lower size limit in mm
U.shad <- 360    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z <- (U.shad-L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####

####Now rerun functions - some depend on delta.z

nm_data <- read.csv('./NMdata.csv',header=TRUE,sep=",")
nm1965_data <- nm_data[nm_data$year==1965,]

Tf <- 30 # number of years

## Use counts
n <- matrix(0,length(zmesh),Tf)
n[,1] <- nm1965_data$count/delta.z

plt.df <- data.frame(length=zmesh, number = n[,1])
ggplot(data=plt.df, aes(x=length, y=number)) +
  geom_bar(stat="identity") + 
  theme_bw()

# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,0],m.par)
round(K*delta.z,2)

#library(plot.matrix)
#plot(K*delta.z)



###### Paper graphs

#### Model Analysis. No comparison to data
N <- 5 # number of size classes
L.shad <- 0.01   # lower size limit in mm
U.shad <- 500.0    # upper size limit in mm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta.z = (U.shad - L.shad)/N
zmesh <-  L.shad + ((1:N) - 1/2) * (U.shad-L.shad)/N    # midpoint of N intervals
####



Tf <- 40 # number of years
n <- matrix(0,length(zmesh),Tf)
n0.total <- 5000
n[,1] <-dnorm(zmesh, mean = (U.shad + L.shad)/2, sd = 50)  
n[,1] <- (n[,1]/sum(n[,1]))*n0.total/delta.z

# Note: sum(n[,1])*delta.z = n0.total

plt.df <- data.frame(length=zmesh, number = n[,1])
ggplot(data=plt.df, aes(x=length, y=number)) +
  geom_bar(stat="identity") + 
  theme_bw()

K <- P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,1],m.par)
round(K*delta.z,2)

for (i in 1:(Tf-1)) {
  K_iter <- (P_z1z(zmesh,zmesh,m.par) + F_z1z(zmesh,zmesh, n[,i],m.par))*delta.z
  n[,i+1] <- K_iter %*% n[,i]
  #   plot(n[,i+1])
}

# Population size vs time
n.total <- rep(0,Tf)
for (i in 1:(Tf-1)){
  n.total[i] <- sum(n[,i])*delta.z
}

plt.df <- data.frame(time.years=seq(Tf), number = n.total)
ggplot(data=plt.df, aes(x=time.years, y=number)) +
  geom_point() + 
  theme_bw()

