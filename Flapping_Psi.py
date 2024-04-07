'''
    Title: Helicopter Rotor Dynamics
    Description: Present the behaviour of blade
    in motion using dynamics relations
'''
import numpy as np
from scipy.integrate import odeint
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
        vibar = np.sqrt(1/(V1 **2 + (V2 + vibar0)**2 ))
        if np.abs((vibar - vibar0) / vibar) <= error:
            return vibar
    return vibar


def solve_beta(psi):
    def beta_derivatives(x, psi, b, bdot, c):
        bd = x[1]
        bdd = -(bdot(psi) * x[1] + b(psi) * x[0]) + c(psi)
        return [bd, bdd]

    global cbdot, cbeta, const
    theta = theta0 - thetalat * np.cos(psi) - thetalong * np.sin(psi)
    cbdot = g * (1 + (4 / 3) * mu * np.sin(psi))
    cbeta = (1 + g * mu * (4 / 3 * np.cos(psi) + mu * np.sin(2 * psi)))
    c0 = g * (theta * ((1 + mu) ** 2) - 4 / 3 * lambdar + 2 / 3 * mu * pbar)  # constant - "constant"
    c1 = g * (mu * (8 / 3 * theta - 2 * lambdar) + pbar) - 2 * qbar  # Constant - sin(psi)
    c2 = g * qbar + 2 * pbar  # Constant - cos(psi)
    c3 = g * (qbar * mu * 2 / 3)  # Constant - sin(2psi)
    c4 = g * (-theta * (mu ** 2) - 2 / 3 * mu * pbar)  # Constant - cos(2psi)
    const = (c0 + c1 * np.sin(psi) + c2 * np.cos(psi) + c3 * np.sin(2 * psi) + c4 * np.cos(2 * psi))
    bdot = cbdot



    def cbdot2(psi):
        cbdot = g * (1 + (4 / 3) * mu * np.sin(psi))

        return cbdot

    def cbeta2(psi):
        cbeta = (1 + g * mu * (4 / 3 * np.cos(psi) + mu * np.sin(2 * psi)))

        return cbeta

    def const2(psi):
        theta = theta0 - thetalat * np.cos(psi) - thetalong * np.sin(psi)
        c0 = g * (theta * ((1 + mu**2)) - 4 / 3 * lambdar + 2 / 3 * mu * pbar)  # constant - "constant"
        c1 = g * (mu * (8 / 3 * theta - 2 * lambdar) + pbar) - 2 * qbar  # Constant - sin(psi)
        c2 = g * qbar + 2 * pbar  # Constant - cos(psi)
        c3 = g * (qbar * mu * 2 / 3)  # Constant - sin(2psi)
        c4 = g * (-theta * (mu ** 2) - 2 / 3 * mu * pbar)  # Constant - cos(2psi)
        const = (c0 + c1 * np.sin(psi) + c2 * np.cos(psi) + c3 * np.sin(2 * psi) + c4 * np.cos(2 * psi))
        return const
    n = len(psi)
    beta0 = [0,0]
    bet, betdot = odeint(beta_derivatives, beta0, psi,  args=(cbeta2,cbdot2, const2)).T
    return np.asarray(bet), np.asarray(betdot)



class helicopter:
    rho = 1.25                    # [kg/m^3]
    m_emp = 1051                  # [kg]
    m_fuel = 405                  # [kg]
    g = 9.81                      # [m/s^2]
    M = m_emp + m_fuel            # [kg]
    W = M * g                     # [N]
    Sfus = 38.584                 # [m^2]
    CDfus =0.002                  # [-]
    R = 5.345                     # [m]
    vtip = 213                    # [m/s]
    omega = vtip / R              # [rad/s]
    A = np.pi * R**2              # [m^2
    vih = np.sqrt((W/(2*rho*A)))  # [m/s]
    n = 3
    c = 0.3                       # [m]

# - 1 - Blade falpping (one complete blade revolution)
# - General variableS
Vh = 20  # Forward velocity [m/s]
rho = helicopter.rho
R = helicopter.R
omega = helicopter.omega   #  [rad/s]
vih = helicopter.vih
W = helicopter.W
Sfus = helicopter.Sfus
CDfus = helicopter.CDfus

# ---- Solving induced velocity
error = 10**(-3)
vibar = glauert(vih, Vh, W, CDfus, Sfus, error)
vi = vibar * vih

c = helicopter.c       #  [m]
Ib = 995 / 3   #  [kg.m^2] - Each blade
Clapha = 5.7    # NACA0012
gamma =  (Clapha * rho * c * R**4) / Ib   #  [-]

# --- Solving for forward flight with lateral and longitudianl cyclic
deg2rad = np.pi/180
q = 20 * deg2rad  # Pitch rate [deg/s]
p = 10 * deg2rad  # Roll pitch rate [deg/s]
alphasp = np.tan(vi/Vh)
theta0 = 6 * deg2rad                         # [rad] Collective pitch angle
thetalong = 2 * deg2rad                      # [rad] Longitudinal cyclic
thetalat = 1 * deg2rad                       # [rad] Lateral cyclic
alphac = thetalong -  alphasp                # [rad] AoA from control plane
mu = Vh * np.cos(alphac) / (omega * R)       # [-]
lambdac = Vh * np.sin(alphac) / (omega * R)  # [-]
lambdai = vi / (omega * R)                   # [-]
lambdar = lambdai + lambdac                  # [-]
pbar = p / omega                             # [-] Dimensionless pitch rate
qbar = q / omega                             # [-] Dimensionless roll rate
g = (gamma) / 8  # Variable     # [(rad/s)^2] - Constant parameter

dpsi = 1
psi = np.linspace(0, 360*5, int(360*5/dpsi) + 1) * deg2rad
# - Roll and pitch
beta, llll = solve_beta(psi)
fig = plt.figure()
fs = 15
plt.xlabel("Ψ [deg]", fontsize=fs)
plt.ylabel("β [deg]", fontsize=fs)
plt.ylim(0,6)
# plt.title("Beta angle in advancing blade")
plt.plot(psi / deg2rad, beta / deg2rad, color="blue")
plt.grid()

# -- alpha variation ---> beta dot = 0, beta ---> static
dpsi = 1
psi = np.linspace(0, 360, int(360/dpsi) + 1) * deg2rad
n = 5
psi2 = np.linspace(0, 360 * n, int(360*n/dpsi) + 1) * deg2rad


betas, betasdot = solve_beta(psi2)
betas_hom = betas[360*(n-1):360*n+1]
fig = plt.figure()
fs = 15
plt.xlabel("Ψ [deg]", fontsize=fs)
plt.ylabel("β [deg]", fontsize=fs)
plt.ylim(0,6)
# plt.title("Beta angle in advancing blade")
plt.plot(psi / deg2rad, betas_hom / deg2rad, color="blue")
plt.grid()

betasdot_hom = betasdot[360*(n-1):360*n+1]
alpha = []
dr = 0.05
r = 1
xr = np.linspace(0.1, r, int(r/dr))
x1 = xr  #  0-1 - r
theta = theta0 - thetalong * np.cos(psi) - thetalat * np.sin(psi)
for x in x1:
    x = float(x)
    den = x + mu * np.sin(psi)
    # print(den)
    alpha1 = theta - (lambdai + betasdot_hom * x + betas_hom * np.cos(psi) * mu + lambdac + x * (qbar * np.cos(psi) + pbar * np.sin(psi))) / den
    alpha.append(alpha1)
alpha = np.asarray(alpha)

# -- Test circle


cos = np.cos(psi)
sin = np.sin(psi)


X = xr[:, None] * sin
Y = xr[:, None] * -cos
# fig = plt.figure(figsize=(4,4))
#
# for xp, yp, pt in zip(X, Y, alpha):
#     for xp2, yp2, pt2 in zip(xp, yp, pt):
#         plt.text(xp2, yp2, s=str(round(pt2 * 180/np.pi)), fontsize=14)
# plt.scatter(X, Y)
# plt.show()
alphan = alpha / deg2rad

mini = np.round(alphan).min()-1
maxi = np.round(alphan).max()+1
dt = 1

lev_max = 8
lev_min =-2
levels = np.linspace(lev_min, lev_max, int(((lev_max-lev_min)/dt))+1)
fig = plt.figure(figsize=(5,5))
fs = 15
j = plt.contourf(X, Y, np.round(alphan), levels=levels, cmap = "jet")
pk = plt.contour(X, Y, np.round(alphan),levels=levels, colors= "black")
plt.clabel(pk, inline=1, fontsize=15)
plt.colorbar(j)
plt.xticks(np.linspace(-1,1,5))
plt.yticks(np.linspace(-1,1,9))
plt.xlabel("x/R [-]", fontsize=fs)
plt.ylabel("y/R [-]", fontsize=fs)
plt.axhline(0, color="k", linestyle="--")
plt.axvline(0,color="k", linestyle="--")

plt.show()


# ------------

theta = theta0 - thetalong * np.cos(psi) - thetalat * np.sin(psi)
a0 = g * theta0 *(1+mu**2) - g * lambdar * 4/3 + g * 2/3 * mu * pbar - g*thetalong * 4/3 * mu
a1 = (8/3 * theta0 * mu - 2 * qbar/g - 2* mu * lambdar + pbar - thetalong*(1+1.5*mu**2)) / (1- (mu**2) / 2)
b1 = (4/3 * a0 * mu - 2 * pbar / g - qbar + thetalat * (1+0.5*mu**2)) / (1 + (mu**2) / 2)
betateste = (a0 - a1 * np.cos(psi) - b1 * np.sin(psi)) / deg2rad
print("a0: ", a0/deg2rad, "a1: ", a1/deg2rad, "b1: ", b1/deg2rad)
fig = plt.figure()
fs = 15
plt.xlabel("Ψ [deg]", fontsize=fs)
plt.ylabel("β [deg]", fontsize=fs)
# plt.title("Beta angle in advancing blade")
plt.plot(psi/deg2rad, betateste, label="$\\beta$ - Fourier Series", color="r")
plt.plot(psi / deg2rad, betas_hom / deg2rad, label="$\\beta$ - Exact Solution", color="blue")
plt.legend(shadow=True, fontsize=fs)
plt.grid()
plt.show()