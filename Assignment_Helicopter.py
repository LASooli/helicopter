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
m_emp = 1051 # [kg]
m_fuel = 405 # [kg]
g = 9.81 # [m/s^2]
M = m_emp +  m_fuel # [kg]
W = M * g
blade_r = 5.345 # [m]
vtip = 213
omega = vtip / blade_r # [rad/s]
n_blades = 3
blade_chord = 0.3 # [m]
# --- Induced velocity
T = W  # [N]
rho = 1.225 # [kg/m^3]
Adisc = np.pi * (blade_r ** 2) # [m^2]
vh = np.linspace(0,100,100) # [m/s]

vi_hov = (T/(2*rho*Adisc)) ** 0.5
vi_for = (-(vh ** 2 ) / 2 + (((vh ** 2) / 2) ** 2 + vi_hov**4) ** 0.5) ** 0.5

# --- Ideal power
Pi = W * vi_hov

# --- Hover power
# - ACT
FM = 0.7
Pact_hover = Pi / FM

# -- BEM
# - Induced power
k = 1.1
P_ind = k * T * vi_hov

# - Profile drag
# Rotor solidity
sigma = n_blades * blade_chord / (np.pi * blade_r)

# Profile drag power
lambdai = vi_hov / vtip
CT = 2 * lambdai **2  # Thrus coefficient
CDp = (1/FM - k) /((sigma/8) / (CT * np.sqrt(CT/2)))

P_p = (sigma * CDp / 8) * rho * ((vtip) ** 3) * np.pi * blade_r **2

# Power hover
Pbem_hover = P_ind + P_p



def plot_vi():
    f1 = plt.figure()
    plt.xlabel('Forward velocity [m/s]')
    plt.ylabel('Induced velocity [m/s]')
    # plt.title("")
    plt.yticks(np.arange(0,10+0.1, 1))
    plt.ylim(0,10)
    plt.grid(which="major", linewidth=1)
    plt.grid(which="minor", linewidth=0.2)
    plt.minorticks_on()
    plt.plot(vh, vi_for, label="Ex")
    # plt.legend(shadow=True, loc="upper center", fontsize=15)
    plt.show()
plot_vi()
