import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("ltrm_fish_data.csv")

#create a new column of year caught
data['sdate']=pd.to_datetime(data['sdate'])
data['year'] =  data['sdate'].dt.year

# extract only la grange data
lg_data = data.loc[(data["site"]=="LA GRANGE") | (data["site"]=="LAGRANGE TWZ") |
         (data["site"]=="LA GRANGE TW") | (data["site"]=="LA GRANGETWZ")]

catch = lg_data.groupby(["year"])["catch"].sum()
year_catch = np.stack(( np.sort(lg_data.year.unique()), np.array(catch)), axis=0)
plt.plot(year_catch[0], year_catch[1],'o')
plt.plot(year_catch[0], year_catch[1])

ltmr_data = pd.read_csv("ltrm_fish_data2.csv")
lg_data = ltmr_data.loc[(ltmr_data["site"]=="LaGrange")]
plt.plot(lg_data['year'],lg_data['catch'])
plt.plot(lg_data['year'],lg_data['electrofish'])

pool26_data = ltmr_data.loc[(ltmr_data["site"]=="Pool26")]
plt.plot(pool26_data['year'],pool26_data['catch'])
plt.plot(pool26_data['year'],pool26_data['electrofish'])