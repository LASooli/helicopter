'''
    Title: Helicopter Rotor Dynamics
    Description: Present the behaviour of blade
    in motion using dynamics relations
'''
import numpy as np
# --- Helicopter parameters

Vh = 20  # Forward velocity [m/s]
q = 20  # Pitch rate [deg/s]
p = 10  # Roll pitch rate [deg/s]

thatablade = 6 * np.pi / 180 # [rad] Collective pitch angle
thetarotor = 2 * np.pi / 180 # [rad] Angle for forward flight
phirotor = 1 * np.pi / 180   # [rad] Angle for lateral motion

# Helicopter didn't change its angle, only the rotor
class helicopter:
    m_emp = 1051  # [kg]
    m_fuel = 405  # [kg]
    g = 9.81  # [m/s^2]
    M = m_emp + m_fuel  # [kg]
    W = M * g
    R = 5.345  # [m]
    vtip = 213
    omega = vtip / blade_r  # [rad/s]
    A = np.pi * R**2
    vih = (W/(2*rho*A)) ** 0.5
    n = 3
    c = 0.3  # [m]


# - 1 - Blade falpping (one complete blade revolution)
# - General variableS

R = helicopter.R
omega = helicopter.omega   #  [rad/s]
vi = (-(Vh ** 2 ) / 2 + (((Vh ** 2) / 2) ** 2 + helicopter.vih**4) ** 0.5) ** 0.5       #  [m/s]
c = helicopter.c       #  [m]
Ib = 995    #  [kg.m^2]
Clapha = 2 * np.pi
gamma =  (Clapha * rho * c * R**4) / Ib   #  [-]
alphac = 0  #  [rad]
theta = 0   #  [rad]
psi = 0     #  [rad]
mu = Vh * np.cos(alphac) / (omega * R)  #  [-]
lambdac = Vh * np.sin(alphac) / (omega * R)  #  [-]
lambdai = vi / (omega * R)  #  [-]


#  Aero damping term
Caerod = (1 + (4 / 3) * mu * sin(psi)) * omega

#  Centrifugal Spring + Aero spring
Kacs = (1 + gamma * (mu * np.cos(psi) / 6 + (mu**2) * np.sin(2*psi))) * (omega**2)
# - Aerodynamic effects
#  Coriollis effect
cori = -2 * q * omega * np.sin(psi)

# - Roll and pitch
c0 = (theta * ((1 + mu)**2) - 4/3 * (lambdai + lambdac)) #  constant - constant
c1 = mu * (8/3 * theta - 2 * (lambdac + lambdai)) #  Constant - sin(psi)
c2 = (q/omega) #  Constant - cos(psi)
c3 = (q/(3 * omega)) #  Constant - sin(2psi)
c4 = (theta * (mu**2))
Ma = (gamma * (omega**2)/8) * (c0 + c1 * np.sin(psi) + c2 * np.cos(psi) + c3 * np.sin(2 * psi) + c4 * np.cos(2 * psi))


#  Coning angle
a0 = (gamma / 8) * (theta * ((1 + mu)**2) - 4/3 * (lambdai + lambdac))

#  Longitudinal disc tilt
a1 = (-16/gamma * q/omega + 8/3 * mu * theta - 2 * mu * (lambdac + lambdai)) / (1 - 1/(2 * mu**2))

#  Lateral disc tilt
b1 = (-q/omega + 4/3 * mu * a0) / (1 + 1/(2 * mu**2))

