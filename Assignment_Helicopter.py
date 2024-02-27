'''
Title: Performance of Helicopter
Description: This code aims to develop the equations of perfomrance of a helicopter
using the estimated parameters
'''
# Libraries
import numpy as np
import matplotlib.pyplot as plt

# --- Helicopter Parameters
# lb2kg = 0.4539
Wemp = 1051 # [kg]
Wfuel = 405 #[kg]
W = Wemp + Wfuel # [kg]
blade_r = 5.345 # [m]
omega = 213 / blade_r # [rad/s]
Cdp = 0.02


# TODO:
#   Análise no xfoil ou xflr5 do aerofólio da helice em duas seções
#   E fazer a média do Cd
n_blades = 3
blade_chord = 0.3 # [m]
# --- Induced velocity
g = 9.81
T = W*g   # [N]
rho = 1.225 # [kg/m^3]
Adisc = np.pi * blade_r ** 2 # [m^2]
vh = np.linspace(0,600,600) # [m/s]

vi_hov = (T/(2*rho*Adisc)) ** 0.5
vi_for = (-(vh ** 2 ) / 2 + (((vh ** 2) / 2) ** 2 + vi_hov**4) ** 0.5) ** 0.5

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










def plot_vi():
    f1 = plt.figure()
    plt.xlabel('Forward velocity [m/s]')
    plt.ylabel('Induced velocity [m/s]')
    # plt.title("")
    plt.yticks(np.arange(0,45+0.1, 5))
    plt.ylim(0,45)
    plt.grid(which="major", linewidth=1)
    plt.grid(which="minor", linewidth=0.2)
    plt.plot(vh, vi_for, label="Ex")
    # plt.legend(shadow=True, loc="upper center", fontsize=15)
    plt.show()
plot_vi()