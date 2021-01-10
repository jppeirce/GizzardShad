import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# import data from Michaletz paper
dens_surv = pd.read_csv('survival_age0_data.csv', names=['density','survival'])
#setting negative values to 0
dens_surv[dens_surv<0]=0 

# scatter plot of data
dens_surv.plot(x='density',y='survival',kind='scatter', color='red')


# define exponential model
def expmodel(x, a, b):
    return a * np.exp(- b * x)

xdata = dens_surv.density.to_numpy()
ydata = dens_surv.survival.to_numpy()
popt, pcov = curve_fit(expmodel, xdata, ydata, bounds=(0, [100., 1.]))

xmesh = np.linspace(start = dens_surv.min(0).density, stop = dens_surv.max(0).density, \
         num = len(dens_surv.density)  )

dens_surv.plot(x='density',y='survival', kind='scatter', \
               color='red');
plt.plot(xmesh, expmodel(xmesh, *popt), 'blue')
plt.ylabel('probability of survival')
plt.title('Surival Age-0 from Michaletz(2010)')

del dens_surv # Done with df