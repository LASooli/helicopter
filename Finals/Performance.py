'''
Title: Performance of Helicopter
Description: This code aims to develop the equations of perfomrance of a helicopter
using the estimated parameters
'''
# Libraries
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
def glauert(vih, Vh, W, CDfus, Sfus, error):
    """
    :param vih:
    :param Vh:
    :param W:
    :param CDfus:
    :param Sfus:
    :param error:
    :return:
    """
    Vbar = Vh / vih
    D = 0.5 * rho * (Vh**2) * Sfus * CDfus
    alphad = np.arcsin(D/W)
    V1 = Vbar * np.cos(alphad)
    V2 = Vbar * np.sin(alphad)
    vibar0 = 1
    ite = 100
    for i in range(ite):
        print(i)
        vibar = np.sqrt(1/(V1 **2 + (V2 + vibar0)**2 ))
        if np.abs((vibar - vibar0) / vibar) <= error:
            return vibar
        vibar0 = vibar
    return vibar


# --- Helicopter Parameters
# lb2kg = 0.4539
g = 9.81
Wemp = 1051 # [kg]
Wfuel = 405 #[kg]
W = (Wemp + Wfuel)*g # [kg]
blade_r = 5.345 # [m]
omega = 213 / blade_r # [rad/s]
Sfus = 38.585
CDfus = 0.02
n_blades = 3
blade_chord = 0.3 # [m]
# --- Induced velocity
T = W   # [N]
rho = 1.225 # [kg/m^3]
Adisc = np.pi * blade_r ** 2 # [m^2]
Vhmax = 287 / 3.6  # [m/s]
vh = np.linspace(0,Vhmax,int(Vhmax/0.1) + 1) # [m/s]
vi_hov = (T/(2*rho*Adisc)) ** 0.5
error = 10**(-3)
vi_bar = []
for Vh in vh:
    vi_bar.append(glauert(vi_hov, Vh, W, CDfus, Sfus, error))

Vcruise = 245 / 3.6 #  [m/s]
vi_forh = glauert(vi_hov, Vcruise, W, CDfus, Sfus, error) * vi_hov



# --- Ideal power
Pi = W * vi_hov

# --- Hover power
# - ACT
FM = 0.65
Pact_hover = Pi / FM

# -- BEM
# - Induced power
k = 1.1
P_ind = k * T * vi_hov

# - Profile drag
# Rotor solidity
sigma = n_blades * blade_chord / (np.pi * blade_r)

# Profile drag power

P_p = (sigma * Cdp / 8) * rho * ((omega * blade_r) ** 3) * np.pi * blade_r **2

# Power hover
Pbem_hover = P_ind + P_p









