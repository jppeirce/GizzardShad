import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## IBM for Gizzard Shad (Fall 2019)

# Gizzard shad (gs) at LaGrange Ranch(LG) before carp invasion 1983-1999
# A dictionary of parameter values
mpar = { # Growth and survival
        # A size- and age-structured model to estimate 
        # fish recruitment, growth, mortality, and gear selectivity
        'grow_max': 394.30,
        'grow_rate': 0.60,
        #
        'grow_sd': 10.0,
        'surv_min': 0.002,
        'surv_max': 0.70,
        'surv_alpha': 120.0,
        'surv_beta': -5.0,
        # Reproduction
        'sage0_int': 26.87216326509913/100,
        'sage0_decay': 0.003058740051047354,
        'recruit_mean': 112.0,
        'recruit_sd': 40.0,
        'egg_viable': 0.002, 
        'repr_int': -73.9542,
        'repr_z': 0.43019586,
        'egg_slope': 15000/7}

N = 20 # number of size classes
L_shad = 0.01   # lower size limit in cm
U_shad = 500.0    # upper size limit in cm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta_z = (U_shad - L_shad)/N
zmesh =  L_shad + (np.arange(N) + 1/2) * (U_shad-L_shad)/N    # midpoint of N intervals
####

################# Growth function ####################
######################################################
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation L(t) = L_inf(1-e^(-r(t-t0)))
# to find L(t+1) = L_inf*(1-e^(-r)) + e^(-r)*L(t)

def normPDF(z,mean,sd):
    return ( 1/(sd*np.sqrt(2*np.pi))*np.exp(-(z-mean)**2/(2*sd**2)) )

def G_z1z(z1,z,mpar):
    mu = mpar['grow_max']*(1-np.exp(-mpar['grow_rate']))+ \
            np.exp(-mpar['grow_rate'])*z
    sig = mpar['grow_sd']
#   pdf of new size z1 for current size = z
    p_den_grow = normPDF(z1, mean = mu, sd = sig)
    return p_den_grow

# Survival
def s_z(z, mpar):
    return mpar['surv_min'] + (mpar['surv_max']-mpar['surv_min'])/(1+ \
         np.exp(mpar['surv_beta']*(np.log(z)-np.log(mpar['surv_alpha']))))
     
############### Density-dependent reproduction ###########
##########################################################    
    
# From R code:
#bodola.data<-data.frame( # from Bodola(1955)
#  z = c(rep(112,13), 225, 236, 282, 285, 292, 293, 305, 322, 328, 343, 348, 363, 355),
#  Repr = c(rep(0,13), rep(1,13)),
#  Eggs = c(rep(0,13), 22405, 96560, 543912, 524580, 211378, 258345, 356713, 406174, 367670, 260509, 267216, 350283, 215331)
#)
# 
#mod.Repr <- glm(Repr ~ z, family = binomial, data = bodola.data)
##mod.Eggs <- glm(Eggs ~ z, family = poisson, data = bodola.data)

#m.par["repr.(Intercept)"] = coef(mod.Repr)[1]
#m.par["repr.z"] = coef(mod.Repr)[2]
    
# Probability of females spawning
def p_bz(z, mpar):
    # linear predictor
    linear_p = mpar['repr_int'] + mpar['repr_z'] * z  
    # logistic transformation to probability with max 90%
    return (1/(1+np.exp(-linear_p)))*0.9  


plt.plot(zmesh,p_bz(zmesh,mpar))
plt.xlabel('length)') 
plt.ylabel('probability of reproduction') 
plt.title('Reproduction Probabality - Gizzard shad') 
ax = plt.gca();
ax.set_yticks(np.arange(0, 1, .1));
plt.savefig('paper/figures/repro_prob.png', bbox_inches='tight')
plt.show()
#plt.bar(x = zmesh, height = p_bz(zmesh,mpar), width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
#plt.show() 

# number of adults in each bin that reproduce    
#plt.bar(x = zmesh, height = p_bz(zmesh,mpar)*n[:,0], width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
#plt.show()

# Using the batch fecundity vs length (figure 1a) from Jons & Miranda (1997)
# linear model
def eggs_z(z, mpar):
    Eggs = np.zeros(len(z))
    Eggs[z>140] = mpar['egg_slope']*(z[z>140]-140)
    return Eggs

plt.plot(zmesh,eggs_z(zmesh,mpar))
plt.xlabel('length') 
plt.ylabel('eggs (hundred thousand)') 
plt.title('Eggs Produced - Gizzard shad') 
ax = plt.gca();
ax.set_yticks(np.arange(0, 800000, 100000));
ax.set_yticklabels(np.arange(0, 9, 1));
#plt.savefig('paper/figures/eggs.png', bbox_inches='tight')
plt.show()

#plt.bar(x = zmesh, height = eggs_z(zmesh,mpar), width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
#plt.show()

# number of eggs produced in each bin    
#plt.bar(x = zmesh, height = eggs_z(zmesh,mpar)*p_bz(zmesh,mpar)*n[:,0], width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
#plt.show()

#def density(n, mpar): #density of VIABLE age0 from current population n
#    age0_dist = mpar['egg_viable']*eggs_z(zmesh,mpar)*p_bz(zmesh,mpar)*n
#    return 10**(-3)*(age0_dist.sum()*delta_z)

def surv_age0(z, n, mpar): #Michaletz (2010)
    #distribution of VIABLE age0 from current population n
    age0_dist = mpar['egg_viable']*eggs_z(z,mpar)*p_bz(z,mpar)*n
    age0_density = 10**(-3)*(age0_dist.sum()*delta_z)
    return mpar['sage0_int']*np.exp(-mpar['sage0_decay']*age0_density)

plt.plot(np.arange(0,1000),mpar['sage0_int']*np.exp(-mpar['sage0_decay']*np.arange(0,1000)))
plt.xlabel('density (age-0 per 10 m^3)') 
plt.ylabel('survival probability') 
plt.title('Survival Probabality for Age-0 - Gizzard shad') 
ax = plt.gca();
ax.set_yticks(np.arange(0, .4, .1));
plt.savefig('paper/figures/surv_age0.png', bbox_inches='tight')
plt.show()



#    return mpar['sage0_int']*np.exp(-mpar['sage0_decay']*density(n,mpar))/100

#dmesh = np.linspace(start=1, stop=5001, num=100)
#plt.plot(dmesh,mpar['sage0_int']*np.exp(-mpar['sage0_decay']*dmesh)/100)

# Recruits PDF
def C_1z1(z1, mpar):
    mu = mpar['recruit_mean']
    sig = mpar['recruit_sd']
#    p_den_recruit = norm.pdf(z1, loc = mu, scale = sig)
    p_den_recruit = normPDF(z1, mean = mu, sd = sig)
    p_den_recruit = p_den_recruit/(p_den_recruit.sum()*delta_z)
    return p_den_recruit


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the growth and survival kernel
def P_z1z(z1, z, mpar):
#    return s_z(z, mpar) * G_z1z(z1, z, mpar)
# Define growth matrix
    G_matrix = np.zeros((N,N))
    for x in range(len(G_matrix)):
        G_matrix[:,x] = G_z1z(zmesh, np.repeat(zmesh[x],N), mpar)
# To ensure that G(z',z)*delta_z is a probability each column has sum 1
        G_matrix[:,x] = G_matrix[:,x]/(G_matrix[:,x].sum()*delta_z)
# and now times survival
    return s_z(z, mpar) * G_matrix


# OLD
#P_matrix = np.zeros((N,N))
#for x in range(len(P_matrix)):
#    P_matrix[:,x] = P_z1z(zmesh, np.repeat(zmesh[x],N), mpar)

## Define the fecundity kernel
def F_z1z(z1, z, n, mpar): 
    age1_dist =  p_bz(z, mpar) * eggs_z(z, mpar) * mpar['egg_viable'] *  \
    surv_age0(zmesh,n,mpar)
    # returns fecuntity kernel (as a matrix). Recruits= F.dot(n*delta_z) 
    return np.outer(C_1z1(z1, mpar), age1_dist)

## OLD
#    return p_bz(z, mpar) * eggs_z(z, mpar) * mpar['egg_viable'] *  \
#    surv_age0(n,mpar)*c_0z1(z1, mpar)

##Define dynamical system
    
# initial length distribution
Tf = 30 # number of years
n = np.zeros((len(zmesh),Tf))
n0_total = 5000
n[:,0]= normPDF(zmesh, mean = (U_shad + L_shad)/2, sd = 50)
n[:,0]= (n[:,0]/n[:,0].sum())*n0_total/delta_z

plt.bar(x = zmesh, height = n[:,0], width = 0.7*delta_z)
plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
plt.show() 

#F_matrix = np.zeros((N,N))

#for x in range(len(F_matrix)):
#    F_matrix[:,x] = F_z1z(zmesh, np.repeat(zmesh[x],N), n[:,0], mpar)


# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,0],mpar)
np.round(K*delta_z,2)

fig = plt.figure()
plt.imshow(K*delta_z, origin = 'upper');
plt.colorbar()
plt.xlabel('current size bin') 
plt.ylabel('next size bin') 
plt.title('Iteration Matrix - Gizzard shad') 
ax = plt.gca();
ax.set_xticks(np.arange(0, N, 1));
ax.set_yticks(np.arange(0, N, 1));
ax.set_xticklabels(np.arange(1, N+1, 1));
ax.set_yticklabels(np.arange(1, N+1, 1));
plt.show()

fig.savefig('figure1.png')

abs(np.real(np.linalg.eigvals(K*delta_z))).max()

# Dynamical system
for x in range(1,Tf):
    K_iter = (P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,x-1],mpar))*delta_z
    n[:,x] = np.dot(K_iter,n[:,x-1])

# Population size vs time
n_total = np.zeros(Tf)
for x in range(0,Tf):
    n_total[x] = n[:,x].sum()*delta_z
plt.plot(np.arange(Tf),n_total)    
    
# surv_age0 vs time
surv_t = np.zeros(Tf)
for x in range(0,Tf):
    surv_t[x] = surv_age0(zmesh, n[:,x],mpar)
    
plt.plot(np.arange(Tf),surv_t);
plt.xlabel('time') 
plt.ylabel('s0') 
plt.title('Survival Probability of Age-0 Gizzard Shad')    
plt.savefig('paper/figures/age0time.png', bbox_inches='tight')
plt.show()

# Now graph
f, ax = plt.subplots(1);
for t in range(0,5):
    ax.plot(zmesh,n[:,t])
 #    plt.bar(x = zmesh, height = n[:,t], width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0))  
plt.xlabel('length') 
plt.ylabel('population density') 
plt.title('Population of Gizzard Shad') 
plt.savefig('paper/figures/popdensity.png', bbox_inches='tight')
plt.show()


# Relative growth vs time
lamb = np.zeros(Tf-1)
for x in range(0,Tf-1):
    lamb[x] = n[:,x+1].sum()/n[:,x].sum()
    
plt.plot(np.arange(Tf-1), lamb);
plt.xlabel('time (years)') 
plt.ylabel('lambda') 
plt.title('Relative Growth Rate - Gizzard Shad') 
plt.savefig('paper/figures/lambda.png', bbox_inches='tight')
plt.show()

############
# Population size vs time
n_total = np.zeros(Tf)
for x in range(0,Tf):
    n_total[x] = n[:,x].sum()*delta_z
plt.plot(np.arange(Tf),n_total)
plt.plot(np.arange(5,30)+1990,n_total[5:30])

############
## Versus LTMR data
####La Grange
# initial length distribution
Tf = 30 # number of years
n = np.zeros((len(zmesh),Tf))
n0_total = 100000
n[:,0]= normPDF(zmesh, mean = (U_shad + L_shad)/2, sd = 50)
n[:,0]= (n[:,0]/n[:,0].sum())*n0_total/delta_z

# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,0],mpar)
np.round(K*delta_z,2)

# Dynamical system
for x in range(1,Tf):
    K_iter = (P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,x-1],mpar))*delta_z
    n[:,x] = np.dot(K_iter,n[:,x-1])

n_total = np.zeros(Tf)
for x in range(0,Tf):
    n_total[x] = n[:,x].sum()*delta_z

ltmr_data = pd.read_csv("./ltrm_fish_data2.csv")
lg_data = ltmr_data.loc[(ltmr_data["site"]=="LaGrange")]
#plt.plot(lg_data['year'],lg_data['catch'])
plt.plot(np.arange(2,27)+1993,n_total[2:27])
plt.plot(lg_data['year'],lg_data['electrofish'],'-ro')
plt.xlabel('time (years)') 
plt.ylabel('total density') 
plt.title('Density vs La Grange Observations - Gizzard Shad') 
plt.savefig('paper/figures/lagrange.png', bbox_inches='tight')
plt.show()



# initial length distribution
Tf = 30 # number of years
n = np.zeros((len(zmesh),Tf))
n0_total = 5000
n[:,0]= normPDF(zmesh, mean = (U_shad + L_shad)/2, sd = 50)
n[:,0]= (n[:,0]/n[:,0].sum())*n0_total/delta_z


# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,0],mpar)
np.round(K*delta_z,2)

# Dynamical system
for x in range(1,Tf):
    K_iter = (P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,x-1],mpar))*delta_z
    n[:,x] = np.dot(K_iter,n[:,x-1])

n_total = np.zeros(Tf)
for x in range(0,Tf):
    n_total[x] = n[:,x].sum()*delta_z

pool26_data = ltmr_data.loc[(ltmr_data["site"]=="Pool26")]
#plt.plot(pool26_data['year'],pool26_data['catch'])
plt.plot(np.arange(5,30)+1990,n_total[5:30]);
plt.plot(pool26_data['year'],pool26_data['electrofish'], '-ro')
plt.xlabel('time (years)') 
plt.ylabel('total density') 
plt.title('Density vs Pool 26 Observations - Gizzard Shad') 
plt.savefig('paper/figures/pool26.png', bbox_inches='tight')
plt.show()
