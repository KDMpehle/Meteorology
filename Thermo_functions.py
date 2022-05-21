# -*- coding: utf-8 -*-
"""
A script with functions useful in
atmospheric thermodynamics

@author: Khaya Mpehle
"""
import numpy as np
import scipy.optimize
# Define some parameters
eps = 0.622 #g/g ratio of dry air specifc gas const to moist gas const.
Lv = 2.5E6 #Latent heat of vaporisation
Cp = 1003.5 #Specific heat capcity of dry air at constant pressure in Jkg-1K-1
Rd = 287.058 # The specific gas constant of dry air in Jkg-1K-1 
g = 9.8 # gravitational acceleration in m/s/s
rhol = 1000 # the rough density of liquid water in kg/m^3

# Calculate the Saturation vapour pressure
def vapour_pres_sat(T):
    '''
    A function for calculating the saturation vapour pressure
    input: 
    T: Temperature in degrees celcius. As a scalar or array of values
    -----
    output:
    es: Saturation vapour pressure in hPa/mb. Returns an array or scalar 
    depending on the temperature input.
    
    notes: The formuala used is August-Roche-Magnus formula for the saturation
    vapour pressure. See "Improved Magnus' form approximation of Saturation
    Vapour Pressure" released by NOAA.
    '''
    #Here you'd need to put in error checking. 
    return 6.1094 * np.exp(17.625 * T / (T + 243.04))
        
def vapour_pres(Td):
    '''
    Function for calculating the vapour pressure given the dew point.
    
    input:
    Td: dew point in degrees celcius. As a scalar or array
    output:
    e: The vapour pressure at the temperature associated with Td. Returns object
    with the same shape as Td
    '''
    #Here put in the error checking.
    return 6.1094 * np.exp(17.625 * Td / (Td + 243.04))
def mix_ratio_sat(P, T):
    '''
    Function to calculate the saturation mixing at a particular pressure and
    temperature.
    -----
    input:
    T: Temperature in degrees celcius. As a scalar or array of values
    P: Associated pressure. As a scalar or array. Must have the same shape
    and data type as the temperature, with pressure and temp of the same level
    at equivalent array index.
    -----
    output:
    The saturation mixing ratio. Shape is the same as the shape of temp and 
    pressure arrays.
    '''
    #Here you need to put in error checking
    eps = 0.622 # g/g ratio of dry specific gas constant to mois gas const
    return eps * vapour_pres_sat(T) / (P - vapour_pres_sat(T))

def mix_ratio(P,Td):
    '''
    Function to calculate the mixing at a particular pressure and
    dew point.
    -----
    input:
    Td: Dew point in degrees celcius. As a scalar or array of values
    P: Associated pressure. As a scalar or array. Must have same shape and data
    type as the temperature, with pressure and temp of the same level in 
    corresponding array positions.
    -----
    output:
    The mixing ratio in [WHAT UNITS?]
    '''
    eps = 0.622
    return eps * vapour_pres(Td) / ( P - vapour_pres(Td))

def isohume(P,r):
    '''
    Function to calculate the temperature-pressure profile when the mixing 
    ratio rs is kept constant. 
    -----
    input:
    P: The pressure in hPa. As a scalar or array of values
    rs: The mixing ratio desired. Given in g/g (?)
    -----
    output:
    Temperature in degreess celsius. Same shape as the inputted pressure. 
    '''
    eps = 0.622 
    a = 6.1094 # constants for the August-Roche-Magnus formula
    c = 17.625
    d = 243.04
    # The actual temperature function
    arg = P * r / (r + eps) / 6.1094 
    return 243.04 * np.log(arg) / (17.625 - np.log(arg))

def dry_adiabats(T0, P0 = 1050., Pf = 100., n = 50) :
    '''
    Script to calculate the dry adiabats
    -----
    input:
    T0: Initially temperature given as 1-d array
    temp_type: Defaults to 'celcius' and converts temperatures to Kelvin.
    P0: The surface pressure. positive scalar. Defaluts to 1050 hPa
    Pf: The highest pressure chosen. Defaults 
    -----
    output: 
    adiabats: an array with shape (n,len(T0)) holding the unsaturated 
    adiabat curves.
    '''
    P = np.linspace(P0, Pf, n)
    #print(len(T0))
    T= np.zeros((m, n)) # initialise array to hold the adiabats
    for k in range(len(T0)):
        T[k,:] = T0[k] * (P[:] / P[0])**3.496  #the 3.496 is Cp/R = 1003.5 / 287.06
    return T 

def SALR(P, T,temp_type = 'celsius'):
    '''
    Function to calculate the saturated adiabatic lapse rate (SALR) at a
    particular pressure and temperature.
    -----
    input:
    T: Temperature. Scalar or array
    P: Pressure in hPa/mb. Scalar or array
    temp_type: defaults to 'celsius' and converts temepratures to Kelvin. If 
    not 'Kelvin' or 'celsius', returns an error.
    -----
    output: The saturated adiabatic lapse rate as C/hPa. same shape as input
    temperature and pressure
    '''
    if temp_type == 'celsius':
        r = mix_ratio_sat(P,T) # get the mixing ratio, with temp given in centigrade
        T = T + 273.15 # convert temperature to Kelvin
    elif temp_type == 'kelvin' or temp_type == 'Kelvin':
        r = mix_ratio_sat(P, T - 273.15) # get the mixing ratio
    else:
        print('error: temp_type must be given as "celsius" or "Kelvin"')
        return
    numerator = Rd * T / Cp + r * Lv / Cp
    denomenator = P * (1 + Lv**2 * r * eps / (Cp * Rd * T**2) )
    #print(Rd / Cp, Lv/Cp, Lv**2 * eps/ (Cp * Rd))
    return numerator / denomenator
def Tlcl_eqn(T,Td,P):
    '''
    A function that sets up an equation relating Tl, the temperature at the LCL
    for a parcel starting at initial point (T0,Td0,p0). The function can be
    used in other scripts to fix T and Td, and then solve for Tl iteratively
    -----
    input:
    T: temperature in celcius. Must be given as a scalar
    Td: dew point in celcius. Must be given as scalar
    P: Associated pressure. Must be given as a scalar
    -----
    output:
    F: a function that should be equal to zero theoretically F(Tl) = 0, 
    and should be solved iteratively to find the temperature at the LCL.
    
    Note: The formula is taken from 
    "The computation of the equivalent potential temperature" (Bolton, 1980)
    '''
    e = vapour_pres(Td) # Get the vapor pressure
    r = mix_ratio(P, Td) # YOU HAVE NOT YET DEFINED r(P, Td)
    T = T + 273.15 #convert the temperature to kelvin
    def F(Tl):
        return (np.log(6.112 / e) + np.log(T / Tl) / (0.2854 * 1 - 2.8E-4 * r)
                + 17.67 * (Tl - 273.15) / (Tl - 29.65))
    return F, e, r, T 

def find_Tlcl(T,Td,P):
    '''
    Finds the Lcl temperature and pressure using Bolton's formula. Written
    to be adapted to a root finding approach if ever need( don't think so,
    Bolton's formula already so accurate)
    -----
    input:
    T: Temperature in celcius. Must be given as a scalar
    Td: Dew point in celcius. Must be given as a scalar
    P: Associated pressure in hPa/mb. Must be given as a scalar
    -----
    output:
    [Tl,Pl]: The temperature and pressure at the lcl. two component array
    '''
    F, e, r, T = Tlcl_eqn(T,Td,P) # Get eqn, e, r and kelvin temp
    Tl = 55 + 2840 / (3.5 * np.log(T) -np.log(e) -4.805) # initial guess
    #print(Tl)
    #Tl = scipy.optimize.fsolve(F, Tl0)[0] # find the lcl temperature
    #Now use Poisson's formula to find the pressure level
    Pl = P * (Tl/T)**3.496  # the 3.496 is Cp/R = 1003.5 / 287.06
    Tl = Tl - 273.15
    return np.array([Tl,Pl])

def moist_adiabat(Ti,P_vals):
    '''
    Solves for the moist saturated adiabat using the moist lapse rate
    dT/dp = SALR(T, Td, P) starting from the level (Ti, Pi) using a 
    simple forward Euler numerical method
    Script assumes lapse rate is C/hPa
    -----
    input:
    Ti: The initial temperature in degrees celcius. Given as scalar
    P_vals: Array of pressure values over which to calculate the saturated adiabats.
    n: The number of points on the saturated adiabat curve
    -----
    output:
    saturated adibatic. a 2 x n array holding the temperature values in the 
    first row (in degrees celcius) and the pressure values in the second row
    (in hPa).
    '''
    #P_vals = np.linspace(Pi,Pf,n) # get the array of P_Vals
    T_vals = np.zeros(len(P_vals))
    deltaP = P_vals[1] - P_vals[0] # get the pressure increment.
    T_vals[0] = Ti # populate first element of the temperature and pressure arrays
    #P_vals[0] = Pi
    for k in range(len(P_vals) - 1):
        #print(SALR(P_vals[k],T_vals[k]))
        T_vals[k+1] = T_vals[k] + SALR(P_vals[k],T_vals[k]) * deltaP
    return T_vals


def lift_trace(T, Td, P, DeltaP):
    '''
    To lift a temperature trace by some pressure increment deltaP.
    The script will lift at the dry adiabatic lapse rate until the
    lifting condensation level is reached (if it is), and then lift
    along a saturated adiabat
    -----
    input:
    T: Array of temperatures in sounding data. given as a one dimensional
    array (1,n) [adapt for scalars?]
    Td: Array of dew points in sounding data. 1D array (1,n).
    P: Array of pressures in the sounding data. 1D array (1,n).
    deltaP: The amount we wish to lift by. enter as hPa/mb.
    enter as positive number(lift by 100 mb means DeltaP = 100)
    -----
    output: 
    T_lifted: Array of the lifted the temperatures
    P_lifted: Corresponding array of pressures.
    '''
    #list to hold initial conditions for the saturated lifting
    P_new = P - DeltaP # get the new 
    T_new = np.zeros(len(T))
    for k in range(len(T)): 
        Tlcl, Plcl = find_Tlcl(T[k],Td[k],P[k]) # get the lcl quantities
        if P_new[k] >= Plcl:
            T_new[k] = (T[k] + 273.15) * (P_new[k] / P[k])**0.286 - 273.15 #Dry adiabat lift to new level.
        else: # if the parcel saturates before lifting is over, finish on Saturated adiabat
            P_interval = np.linspace(Plcl, P_new[k],100) # figure out the interval size later
            T_satcurve = moist_adiabat(Tlcl,P_interval)
            T_new[k] = T_satcurve[-1] # only interested final temp after
    return T_new, P_new # return the new temperatures and 

def wetbulb_trace(T, Td, P, dP = 2.):
    '''
    A function to find the wet bulb temperature using a trace. I do not know
    of a simple formula, iterative solver, etc for finding wet bulb
    temperatures, so I simply lift to the LCL then descend down a saturated adiabat
    back to the original pressure.
    -----
    input:
    T: Array of temperatures in sounding data. given as a one dimensional
    array (1,n) [adapt for scalars?]
    Td: Array of dew points in sounding data. 1D array (1,n).
    P: Array of pressures in the sounding data. 1D array (1,n).
    dP: The resolution of the saturated adiabats, defaults to
    2 hPa
    -----
    output:
    Tw: The wetbulb potential temperature trace. equal shape to the 
    temperature
    '''
    Tw = np.zeros(len(T)) # initialise the wetbulb trace array.
    #print(Tw)
    for k in range(len(T)):
        Tlcl, Plcl = find_Tlcl(T[k], Td[k], P[k])
        n_steps = int(np.floor((P[k]-Plcl) / dP)) + 1 # number of steps for saturated adiabats
        print(n_steps)
        P_vals = np.linspace(Plcl, P[k], n_steps) # array of pressure values over which to calculate the lcl
        T_satcurve = moist_adiabat(Tlcl, P_vals)
        Tw[k] = T_satcurve[-1] # The final temperature gives the wetbul temperature. 
    return Tw
def precip_water(T, Td, P):
    '''
    Calculates the precipitable water through a trace. Numerical method
    utilised is the trapezoidal rule as implemented in numpy.
    -----
    input:
    T: Array of temperatures in sounding data. given as a one dimensional
    array (1,n).
    Td: Array of dew points in sounding data. 1D array (1,n).
    P: Array of pressures in the sounding data. 1D array (1,n).
    -----
    output:
    Pw: The precipitable water in mm. A sclar
    '''
    r = mix_ratio(P, Td) # calculate the mixing ratio of each data point
    return -np.trapz(r, x = P) / (rhol * g) # negative sign is to get positive result since we integrate from bottom to top of the atm.

            