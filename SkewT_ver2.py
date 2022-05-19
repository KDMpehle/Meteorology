# -*- coding: utf-8 -*-
"""
Prototype functions for a skewT diagram in matplotlib

@author: Khaya Mpehle
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Thermo_functions import *

#from Thermo_functions import *
#define some default inputs
T0_DA = [n* 10 for n in range(-10,8)] #P0 temperatures for the dry adiabats
T0_SA = [n*5 for n in range(-4,9)] #P0 temperatures of the saturated adiabats
rs_init = [0.1, 0.4, 1, 2, 3, 5, 8, 12, 20] 
P0 = 1050 # in hPa, the lowest pressure height
Pf = 100 # in hPa, the highest pressure height 
rotation = 35
n = int(np.floor((P0-Pf)/50.))+1 # number of points in the dry adiabats
#Define the transform for the skewT variable
def skewT(T, P, rotation = rotation, P0 = P0): 
    '''
    A function to generate the X-axis variable of the skewT
    plot.
    -----
    input:
    T: The temperature. a 1-d array.
    P: The pressure. a 1-d array. Shape must matche the temperature's
    rotation: The rotation that skews the temperature variable. defaults
    at 1050.
    P0: The reference pressure taken as the Y = 0 on the diagram. Given as 
    hPa. 
    -----
    output: 
    X: the skewed temperature variable. 1-d array of equal shape to
    the temperature and pressure arrays.
    '''
    return T + rotation * np.log(P0 / P)

### This is just a test ###
def blank_skewT(T0_DA = T0_DA, rs = rs_init,T0_SA = T0_SA, P0 = P0, Pf = Pf, n = n):
    '''
    Plot a barebones plot with dry adiabats
    -----
    input:
    T0: an array of initial temperatures 
    P0: lowest pressure level. Given in hPa
    Pf: highest pressure level. Given in hPa
    n: number of points in the dry adiabats.
    -----
    output:
    An axis object showing the adiabats and isohumes
    '''
    #adiabats = dry_adiabats(T0, P0, Pf, n) # get the adiabats as an array
    P = np.linspace(P0, Pf, n) # i will reverse and log these later. 
    rs_gg  = np.array(rs) / 1000 # rs in g/g
    fig, ax = plt.subplots(figsize = (10,10)) #initialise an axes object.
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_yticks(P)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.grid(True)
    ax.margins(y = 0) # Zero margin between the diagram's xaxis and tick labels
    ax.set_xlim(-35,45)
    for index, label in enumerate(ax.yaxis.get_ticklabels()):
        if index % 2 == 0: # this hides the pressure tick labels ending in 50 eg 1050, 950, 750 etcl   
            label.set_visible(False)
    for index, line in enumerate(ax.yaxis.get_major_ticks()):
        if index % 2 == 0 or index == n-1: # this hides the tick line of labels ending in 50 and at the top of the plot.
            line.tick1On = False # a way to remove the axis tick while retaining the grid line
    #X = skewT(adiabats, P)
    for k in range(len(T0_DA)): # loop over each adiabat, plot them
        adiabat = (T0_DA[k] + 273.15) * (P / P[0])**0.286 -273.15 # convert to Kelvin for Poisson eqn then convert to celsius
        X_adiabat = skewT(adiabat, P)
        X_isotherm = skewT(T0_DA[k], P)
        ax.plot(X_adiabat, P, color = 'g', alpha = 0.5 )
        if abs(T0_DA[k]) < 1E-3:
            ax.plot(X_isotherm, P, color = 'k', alpha = 0.5)
        else:
            ax.plot(X_isotherm, P, color = 'orange', alpha = 0.5)
            ax.annotate('{}'.format(T0_DA[k]), (X_isotherm[-1], P[-1]), color = 'orange')
    for k in range(len(rs)):
        #P_trunc = P[ P >=400]
        isohumes = isohume(P, rs_gg[k]) # get the temperature isohumes
        X_isohume = skewT(isohumes, P)
        ax.plot(X_isohume, P, '--', color = 'b', alpha = 0.5)
        ax.annotate('{}'.format(rs[k]), (X_isohume[-int(n/2)],P[-int(n/2)]), color = 'b')
    n_sat = int(np.floor((P0-Pf)/2.))+1 # number of points for the saturated adiabats
    P_sat = np.linspace(P0, Pf, n_sat) # get the pressure intervals integrating saturated adiabat
    for k in range(len(T0_SA)): # loop over and plot each saturated adiabat
        saturated_adiabat = moist_adiabat(T0_SA[k], P_sat)
        X_moist_adiabats = skewT(saturated_adiabat[0::25], P)
        ax.plot(X_moist_adiabats, P, '--', color='g', alpha = 0.5)
        ax.annotate('{}'.format(T0_SA[k]), (X_moist_adiabats[-2], P[-2]), color = 'g')
    return ax

if __name__ == "__main__":  # an execution guard if importing to other scripts
    #a sample trace to plot
    P = np.array([900, 850, 800, 700, 600, 500])
    T = np.array([15, 11.8, 9.2, 2.6, 2, -5.3])
    Td = np.array([8.5,3.8,7.2, -2.4, -36.5, -50.3])
    T_lifted,P_new = lift_trace(T,Td,P, 150)
    Tw = wetbulb_trace(T, Td, P) #get wet bulb trace
    #print(P_new)
    t = skewT(T, P)
    td = skewT(Td, P)
    tw = skewT(Tw, P)
    t_lifted = skewT(T_lifted,P_new)
    ax = blank_skewT()
    ax.plot(t, P, color = 'r')
    ax.plot(td, P, color = 'r')
    ax.plot(tw, P, color = 'g')
    #ax.plot(t_lifted,P_new,'-o',color = 'magenta')
    plt.show()
        
    #20 degree saturated potential temperature adiabat test
    #n_sat = int(np.floor((P0-Pf)/2.))+1 # number of points for the saturated adiabats
    #P_sat = np.linspace(P0, Pf, n_sat) # get the standard pressure intervals
    #P = np.linspace(P0,Pf,n)
    #print(P_sat[0::50])
    #sat_adiabat = moist_adiabat(20, P_sat) # 20 degrees at 1000 hPa
    #print(len(sat_adiabat[0::50]))
    #X = skewT(sat_adiabat[0::50], P)
    #print(sat_adiabat[0::50],X)

   
    
