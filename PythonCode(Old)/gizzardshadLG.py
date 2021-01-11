import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## IBM for Gizzard Shad (Spring 2020)

mpar = { # Growth and survival
        # A size- and age-structured model to estimate 
        # fish recruitment, growth, mortality, and gear selectivity
        'grow_max': 394.30,
       # 'grow_max': 300,
        'grow_rate': 0.60,
        'grow_sd': 10.0,
        'surv_min': 0.002,
        'surv_max': 0.70,
        'surv_alpha': 120.0,
        'surv_beta': -5.0,
        # Reproduction
        'sage0_int': 26.87216326509913,
        'sage0_decay': 0.003058740051047354,
        'recruit_mean': 112.0,
        'recruit_sd': 40.0,
        'egg_viable': (0.002), 
        'repr_int': -73.9542,
        'repr_z': 0.43019586,
        'egg_slope': 15000/7}

#### Using data from Life History and Ecology of the Gizzard Shad, Dorosoma 
#### cepedianum (Le Sueur) with Reference to Elephant Butte Lake

N = 11  # number of size classes
L_shad = 60   # lower size limit in cm
U_shad = 390    # upper size limit in cm - we want this to be larger than L-infty

####
# As we test the model, we will need some equal spaced meshpoints. 
# This can be moved down later.
delta_z = (U_shad - L_shad)/N
zmesh =  L_shad + (np.arange(N) + 1/2) * (U_shad-L_shad)/N    # midpoint of N intervals
####

################# Growth function ####################
######################################################
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation L(t) = L_inf(1-e^r(t-t0))
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
    
# Probability of females spawning
def p_bz(z, mpar):
    # linear predictor
    linear_p = mpar['repr_int'] + mpar['repr_z'] * z  
    # logistic transformation to probability with max 90%
    return (1/(1+np.exp(-linear_p)))*0.9  

# Using the batch fecundity vs length (figure 1a) from Jons & Miranda (1997)
# linear model
def eggs_z(z, mpar):
    Eggs = np.zeros(len(z))
    Eggs[z>140] = mpar['egg_slope']*(z[z>140]-140)
    return Eggs

def surv_age0(n, mpar):
    #distribution of VIABLE age0 from current population n
    age0_dist = mpar['egg_viable']*eggs_z(zmesh,mpar)*p_bz(zmesh,mpar)*n
    age0_density = 10**(-3)*(age0_dist.sum()*delta_z)
    return mpar['sage0_int']*np.exp(-mpar['sage0_decay']*age0_density)/100

# Recruits PDF
def c_0z1(z1, mpar):
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
# Define growth matrix
    G_matrix = np.zeros((N,N))
    for x in range(len(G_matrix)):
        G_matrix[:,x] = G_z1z(zmesh, np.repeat(zmesh[x],N), mpar)
# To ensure that G(z',z)*delta_z is a probability each column has sum 1
        G_matrix[:,x] = G_matrix[:,x]/(G_matrix[:,x].sum()*delta_z)
# and now times survival
    return s_z(z, mpar) * G_matrix


## Define the fecundity kernel
def F_z1z(z1, z, n, mpar): 
    age1_dist =  p_bz(z, mpar) * eggs_z(z, mpar) * mpar['egg_viable'] *  \
    surv_age0(n,mpar)
    # returns fecuntity kernel (as a matrix). Recruits= F.dot(n*delta_z) 
    return np.outer(c_0z1(z1, mpar), age1_dist)

# initial length distribution from NW data
lg_data = pd.read_csv("./LGdata.csv")
lg2002_data = lg_data.loc[(lg_data["year"]==2002)]
lg2005_data = lg_data.loc[(lg_data["year"]==2005)]

Tf = 30 # number of years

## Use counts
n = np.zeros((len(zmesh),Tf))
n[:,0] = lg2002_data['count']/delta_z


plt.bar(x = zmesh, height = n[:,0], width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0)) 
plt.show() 


# INITIAL Gizzard Shad Matrix K
K = P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,0],mpar)
np.round(K*delta_z,2)

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


abs(np.real(np.linalg.eigvals(K*delta_z))).max()

#### Dynamical system
##Define dynamical system

for x in range(1,Tf):
    K_iter = (P_z1z(zmesh,zmesh,mpar) + F_z1z(zmesh,zmesh, n[:,x-1],mpar))*delta_z
    n[:,x] = np.dot(K_iter,n[:,x-1])

# Now graph
f, ax = plt.subplots(1)
for t in range(0,3):
 #   ax.plot(zmesh,n[:,t])
     plt.bar(x = zmesh, height = n[:,t], width = 0.7*delta_z)
#plt.xticks(ticks = np.round(zmesh,0), labels = np.round(zmesh,0))  
plt.xlabel('length') 
plt.ylabel('population size') 
plt.title('Population of Gizzard Shad') 

### Plot simulation vs NM data
lg2005_data = lg_data.loc[(lg_data["year"]==2005)]

ax = plt.subplot(111)
ax.bar(x = zmesh-10, height = lg2002_data['count']/lg2002_data['count'].sum(), \
       width = 10, color='b', align='center', label='2002')
ax.bar(x = zmesh, height = lg2005_data['count']/lg2005_data['count'].sum(), \
       width = 10, color='g', align='center',  label='2005')
ax.legend()
plt.xlabel('length')
plt.ylabel('frequency')
plt.title('Gizzard Shad LTMR Data - La Grange')
plt.savefig('paper/figures/LGdata.png', bbox_inches='tight')
plt.show()

ax = plt.subplot(111)
ax.bar(x = zmesh-10, height = n[:,0]/n[:,0].sum(), width = 10, color='b',\
       align='center',label='2002')
ax.bar(x = zmesh, height = n[:,3]/n[:,3].sum(), width = 10, color='g', \
       align='center', label='2005')
ax.legend()
plt.xlabel('length')
plt.ylabel('frequency')
plt.title('Gizzard Shad Simulation - La Grange')
plt.savefig('paper/figures/LGsim.png', bbox_inches='tight')
plt.show()

ax = plt.subplot(111)
ax.bar(x = zmesh-10, height = lg2005_data['count']/lg2005_data['count'].sum(),\
       width = 10, color='r', align='center', label='LTMR Data')
ax.bar(x = zmesh, height = n[:,3]/n[:,3].sum(), width = 10, color='b', \
       align='center', label='Simulation')
ax.legend()
plt.xlabel('length')
plt.ylabel('frequency')
plt.title('Gizzard Shad 2005 (with initial = 2002) - La Grange')
plt.savefig('paper/figures/LGdata_sim.png', bbox_inches='tight')
plt.show()

#
#
#
#
ax = plt.subplot(111)
ax.bar(x = zmesh-10, height = lg2002_data['count'], \
       width = 10, color='b', align='center', label='2002')
ax.bar(x = zmesh, height = lg2005_data['count'], \
       width = 10, color='g', align='center',  label='2005')
ax.legend()
plt.title('Gizzard Shad LTMR Data - La Grange')
#plt.savefig('paper/figures/LGdata.png', bbox_inches='tight')
plt.show()




