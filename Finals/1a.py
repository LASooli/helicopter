'''
Title: Performance of Helicopter
Description: This code aims to display the power distributproduce the total power of our given helicopter
This is used to find the speed for best endurance and best range.
Authors: Matteus and Nathan
'''
# Libraries
import numpy as np
import matplotlib.pyplot as plt


# --- Helicopter Parameters
Wemp = 1051                                             # [kg]
Wfuel = 405                                             # [kg]
W = Wemp + Wfuel                                        # [kg]
blade_r = 5.345                                         # [m]
omega = 213 / blade_r                                   # [rad/s]
Cdp = 0.0053
n_blades = 3
blade_chord = 0.3                                       # [m]
g = 9.81
T = W*g                                                 # [N]
rho = 1.225                                             # [kg/m^3]
Adisc = np.pi * blade_r ** 2                            # [m^2]
sigma = n_blades * blade_chord / (np.pi * blade_r)      # Solidity


# -- Tail rotor values
Rtr = 1.86/2                                            # [m]
chord_tr = 0.185                                        # [m]
sigma_tr = 0.127                                        # Solidity
n_blades_tr = 2
omega_tr = 199 / Rtr                                    # [rad/s]
ktr = 1.3
ltr = 5.78                                              # [m]
Adisc_tr = np.pi * Rtr ** 2                             # [m^2]


def mu_(V, omega, R):
    return V / (omega*R)

# -- Tail Rotor Thrust
def Ttr_(Phov, omega_tr, ltr):
    return Phov / (omega_tr * ltr)

# -- Tail Rotor Power
def Ptail(ktr, Ttr, vitr, Rtr, sigma_tr, Cdp, omega_tr, mu_tr):
    return 1.1 * ktr * Ttr * vitr + (sigma_tr*Cdp)/8*rho*(omega_tr*Rtr)**3*np.pi*Rtr**2*(1+4.65*mu_tr**2)


max_velocity = 80                       # 287 kph = 80 m/s
num_iterations = 200                    # high num_iterations for accurate best range speed
vh = np.linspace(0, max_velocity, num_iterations) 

# --- Induced velocity
vi_hov = (T/(2*rho*Adisc)) ** 0.5
vi_for = (-(vh ** 2 ) / 2 + (((vh ** 2) / 2) ** 2 + vi_hov**4) ** 0.5) ** 0.5

# --- Ideal power
Pi = W * vi_hov

# --- Hover power
# -- ACT
FM = 0.65
Pact_hover = Pi / FM

# -- BEMs
k = 1.1
P_ind = k * T * vi_hov                                              # Induced power

# - Profile drag Power
P_p = (sigma * Cdp / 8) * rho * ((omega * blade_r) ** 3) * np.pi * blade_r **2

Pbem_hover = P_ind + P_p                                            # Power in hover

# --- Total Power in Forward Flight

k_for = 2 
eq_flt_plt_area = 0.5                                               # Graph lec notes p 52

meu = vh / (omega * blade_r)                                        # Advance ratio

Ppro = P_p * (1 + 4.65 * meu ** 2)                                  # Profile Power

Pi = k * T * vi_for                                                 # Induced Power

Ppar = eq_flt_plt_area * 0.5 * rho * vh**3                          # Parasite Fuselage Power

P_for = Ppro + Pi + Ppar                                            # Total Forward Power

# -- Tail Rotor
Ttr = Ttr_(P_for, omega, ltr)                                       # Tail Rotor Thrust
vi_tr = (Ttr/(2*rho*Adisc_tr)) ** 0.5
mu_tr = mu_(vh, omega_tr, Rtr)
Ptr = Ptail(ktr, Ttr, vi_tr, Rtr, sigma_tr, Cdp, omega_tr, mu_tr)   # Tail Rotor Power
# add tail rotor power to total power
P_for += Ptr

# take min value of total power curve
min_y_index = np.argmin(P_for)
v_min_power = vh[min_y_index]
print(f"Speed for best endurance = {v_min_power:.2f} m/s")      # Speed for best endurance

# Calculate tangent on power curve which intersects near (0,0)
gradients = np.diff(P_for) / np.diff(vh)
cs = []
for i, gradient in enumerate(gradients):
    # Get gradient of line at each point
    c = P_for[i] - ( gradient * vh[i])
    cs.append(abs(c))
v_best_range = vh[cs.index(np.min(cs))]
print(f"Speed for best range = {v_best_range:.2f} m/s")         # Speed for best range



plt.xlabel('Forward velocity [m/s]')
plt.ylabel('Power [kW]')
plt.plot(vh, P_for/1000, label='Total') #/1000 to get into kW
plt.plot(vh, Ppro/1000, label='Profile')
plt.plot(vh, Pi/1000, label='Induced')
plt.plot(vh, Ppar/1000, label='Parasite')
plt.plot(vh, Ptr/1000, label='Tail Rotor Power')
plt.legend()
plt.show()