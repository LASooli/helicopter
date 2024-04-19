import matplotlib
import numpy as np
matplotlib.use('TkAgg')
from numpy import pi, cos, sin, arctan
import matplotlib.pyplot as plt
from Trim_function import getTrimAngles


# ---- PHUGOID ANALYSIS
# -- Helicopter parameters

class helicopter:
    rho = 1.225                    # [kg/m^3]
    m_emp = 1051                  # [kg]
    m_fuel = 405                  # [kg]
    g = 9.81                      # [m/s^2]
    M = m_emp + m_fuel            # [kg]
    W = M * g                     # [N]
    Sfus = 38.585                 # [m^2]
    CDfus = 0.02                  # [-]
    R = 5.345                     # [m]
    vtip = 213                    # [m/s]
    omega = vtip / R              # [rad/s]
    A = pi * R**2              # [m^2
    vih = np.sqrt((W/(2*rho*A)))  # [m/s]
    n = 3
    c = 0.3                       # [m]
    Iy = 5313
    cla = 5.7
    Ib = 995 / 3   #  [kg.m^2] - Each blade
    gamma = (cla * rho * c * R**4) / Ib
    h = 1.2

# -- Basic parameters
omega = helicopter.omega
cla = helicopter.cla
CD = helicopter.CDfus
Sfus = helicopter.Sfus
rho = helicopter.rho
g = helicopter.g
W = helicopter.W
n = helicopter.n
R = helicopter.R
c = helicopter.c
sigma = n * c / (np.pi * R)
Iy = helicopter.Iy
m = helicopter.M
h = helicopter.h
gamma = helicopter.gamma
A = helicopter.A

# -- Trim condition
V = 46.3  # [m/s]
[thetac, theta0, lambdai] = getTrimAngles(V)

a1 = float(thetac)
thetac = float(thetac)
theta0 = float(theta0)
lambdai = float(lambdai)

# -- Initial condition
D = 0.5 * CD * Sfus * V**2 * rho
T = np.sqrt((W**2 + D**2))
thetaf = np.arctan(-D/W)
u = V * cos(thetaf)
u0 = u
w = V * sin(thetaf)
w0 = w
thetaf0 = thetaf
q = 0
Ct = T/(rho * ((omega*R)**2) * A)

alphac = thetac - arctan(w/u)
mu = V/(omega*R) * cos(alphac)
lambdac = V * sin(alphac) / (omega*R)
thetac = float(a1)


# --- Loop
U = []
Fx = []
Fm = []
Q = []
pertubation = 0.001
iter = 20
dt = 0.00001
t = np.linspace(0, dt*iter, iter)
glau = True
for i in range(iter):
    if i == 2:
        q = q + pertubation
    print(q)
    V = np.sqrt(w**2 + u**2)
    D = 0.5 * CD * Sfus * (V ** 2) * rho
    alphac = thetac - arctan(w / u)
    thetaf = np.arctan(-D / W)
    mu = V / (omega * R) * cos(alphac)
    lambdac = V * sin(alphac) / (omega * R)
    if glau == True:
        lambdai0 = lambdai
        error = 10**(-16)
        iter = 10000
        for i in range(iter):
            a1 = (8 / 3 * mu * theta0 - 2 * mu * (lambdac + lambdai0) - 16 / gamma * q / omega) / (1 - 0.5 * mu ** 2)
            V1 = (V/(omega*R) * cos(alphac-a1))
            V2 = (V/(omega*R) * sin(alphac-a1) + lambdai0)
            Ctglau = 2*lambdai0*np.sqrt(V1**2 + V2**2)
            Ctbem = 1/4 * cla * sigma * (2/3*theta0*(1+3/2*mu**2)-(lambdac+lambdai0))
            if (np.abs((Ctglau - Ctbem)/Ctbem)) < error:
                lambdai = lambdai0
                Ct = Ctbem
                a1 = (8/3 * mu*theta0 - 2*mu*(lambdac+lambdai)-16/gamma * q/omega)/(1-0.5*mu**2)
                break
            else:
                Ct = Ctbem
                lambdai0 = (Ct/2) / np.sqrt(V1**2 + V2**2)
    T = Ct * (rho * (omega*R) ** 2 * A)
    udot = -g * sin(thetaf) - D/m * u/V + T/m * sin(thetac - a1) -q*w
    wdot = g*cos(thetaf) - D/m * w/V - T/m * cos(thetac - a1) + q*u
    qdot = -T/Iy * h * sin(thetac - a1)
    thetafdot = q
    # X = T * sin(thetac - a1) - D * u / V
    # M = -T * h * sin(thetac - a1)
    u += udot*dt
    w += wdot*dt
    q += qdot*dt
    thetaf += thetafdot*dt
    DX = m*(udot+w*q)
    DM = Iy*qdot
    U.append(u)
    Fx.append(DX)
    Fm.append(DM)
    Q.append(q)

U = np.asarray(U)
Fx = np.asarray(Fx)
Fm = np.asarray(Fm)
# -- When pertubation is in u:
Xum = (Fx[2]-Fx[1])/pertubation
xum = Xum/m

Muiy = (Fm[2]-Fm[1])/pertubation
muiy= Muiy/Iy

# -- When pertubation is in q:
Mqiy = (Fm[2]-Fm[1])
mqiy = Mqiy/pertubation /Iy


# --- Final values computed
Xum = -13.864
Muiy = 16.243
Mqiy = -1370
xum = Xum/m
muiy = Muiy/Iy
mqiy = Mqiy/Iy
g = 9.81
# -- Eigenvalues from phugoid mode
matrix = [[xum, 0, -g], [muiy, mqiy, 0], [0, 1, 0]]
eig, v = np.linalg.eig(matrix)
omegan = abs(eig[1])/(np.sqrt(2))
zeta = np.real(eig[1])/omegan