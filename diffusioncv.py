#import
import numpy as np 
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from typing import Tuple, Optional, Iterable
from scipy.optimize import curve_fit

    
c0 = 5
compound = 'meept'
solvent = 'acn'
salt = 'tbapf6'
scanrates = [0.02,0.025,0.04,0.05,0.08, 0.1, 0.16, 0.2, 0.25, 0.32]
exp = [1]
bounds = [-0.18,0.25]


def CVdiffusion(c0,compound, solvent, salt, scanrate, runs, bounds):

    for run in runs:
        path = f'{compound}_{solvent}_{salt}_{scanrate}_{run}.txt'
        data = pd.read_csv(path, delim_whitespace=True)

        if compound == 'azobenzene':
            time = data.iloc[:, 0]
            current = data.iloc[:, 1]
            potential = data.iloc[:, 2]

        else:
            time = data.iloc[:, 2]
            current = data.iloc[:, 1]
            potential = data.iloc[:, 0]
        
        scan = data.iloc[:, 3]
        data = np.zeros((int(1/3*len(scan)), 3))
        for i in range(0, int(len(scan)*2/3)):
            if scan[i] != 1:
                data[i - int(2/3*len(scan))] = [time[i], (current[i]+current[i + int(1/3*len(scan))])/2, (potential[i]+potential[i +int(1/3*len(scan))])/2]
        
        ymax = max(data[:,1])
        ydata = []
        xdata = []
        for j in range (len(data[:,2])):
            ydata.append(data[j,1])
            xdata.append(data[j,2])
            if data[j,1] == ymax:
                xmax = data[j,2]
                break
        
        
        ydata_fit = []
        xdata_fit = []
        for j in range (0, len(data[:,2])):
            if bounds[0]< data[j, 2] < bounds[1]:
                ydata_fit.append(data[j,1])
                xdata_fit.append(data[j,2])
        
        data_f = np.column_stack((xdata_fit, ydata_fit))
        data_fi = data_f[np.lexsort((-data_f[:, 0], data_f[:, 1]))]
        data_fit = data_fi[len(data_fi) // 2:]
        
        def yfunc(x, a, b):
            return a*x+b
        
        popt, pcov = curve_fit(yfunc, data_fit[:, 0], data_fit[:, 1])
        
    
        i = ymax - yfunc(xmax, popt[0], popt[1])
        return i

ipsv = []
for scan in scanrates:
    ipsv.append(CVdiffusion(c0, compound, solvent, salt, scan, exp, bounds))
    

mu_root = np.sqrt(scanrates)
def yfunc(x, a, b):
            return a*x+b
popt, pcov = curve_fit(yfunc, mu_root, ipsv)
plt.scatter(mu_root,ipsv)
plt.plot(mu_root, yfunc(mu_root, popt[0], popt[1]))
plt.show()

F = 96485
R = 8.3145
T = 298.15
A = 7e-6
#c0 = float(input(f'starting concentration analyte: '))
D0 = (popt[0] / (0.446*F*A*c0))**2 * R*T / F
print(D0)
    