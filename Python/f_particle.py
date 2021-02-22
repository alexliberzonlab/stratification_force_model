# f_particle calculates the trajectory of an inertial particle traversing
# a stratified fluid layer. The z-axis is defined downwards, and
# layer 1 (the top layer) is for z<zu and layer 2 (the bottom layer) is for
# z>zl. 
#
# SYNTAX [t, zp, V] = f_particle(z0, tend, rhop, d, g, 
#                                zu, zl, rho1, rho2, nu1, nu2, fitfun, 
#                                options)
#
# INPUT ARGUMENTS
# ---------------
#        z0 : initial particle position        [m]
#      tend : simulation time                  [s]
#      rhop : particle density                 [kg/m3]
#         d : particle diameter                [m]
#         g : gravitational acceleration       [m/s2]
#        zu : start of density interface       [m]
#        zl : end of density interface         [m]
#      rho1 : density of top layer             [kg/m3]
#      rho2 : density of bottom layer          [kg/m3]
#       nu1 : kinematic viscosity of top layer [m2/s]
#       nu2 : kinematic viscosity of bottom layer [m2/s]
#    fitfun : a function of with arguments (z, zu, zl, f1, f2) to create
#             continuous density and viscosity profiles
#   options : ODE option function
#
# OUTPUT ARGUMENTS
# ----------------
#         t : time interval     [s]
#        zp : particle position [m]
#         V : particle velocity [m/s]

from settling_velocity import *
import numpy as np
from scipy.integrate import solve_ivp # initival value problem


tzl = 0



def f_particle(z0, tend, rhop, d, g, zu, zl, rho1, rho2, nu1, nu2, lam, options):

    # Vp(d) exists in settling_velocity
    h  = zl - zu                                     # interface thickness  [m]
    N  = (2*g*(rho2 - rho1)/h/(rho1 + rho2))**0.5    # buoyancy frequency   [1/s]

    # Construct density and viscosity functions
    def rho(z): return rho2 - 0.5*(rho2-rho1)*(1-np.tanh((z-0.5*(zl + zu))/(lam*h)))
    def nu(z): return  nu2 + 0.5*(nu1-nu2)*(1-np.tanh((z-0.5*(zl + zu))/(lam*h)))

    # Auxiliary functions
    # Cd(Re) exists
    def trec_d2nu2(Re): return 13./Re                                # recovery time 
    # def Vc0_f(Fr,Re,Vp): return 0.13*Fr**0.75 * Vp + 0.*Re           # caudal volume
    # see Marco e-mail from Feb. 21, 2021
    def Vc0_f(Fr,Re,Vp): return (1 - (1.85 / Fr )) * Vp + 0.*Re

    # calculate settling velocity
    sol = settling_velocity(rhop,rho1,g,d,nu1)
    V1 = sol.x
    sol = settling_velocity(rhop,rho2,g,d,nu2)
    V2 = sol.x

    # Determine Fr1, Vc0 and trec
    Fr1  = np.abs(V1) / (N * d)
    Re1  = np.abs(V1) * d / (nu1)
    trec = trec_d2nu2(V2*d/nu2)*(d**2/nu2)
    Vc0  = np.max(Vc0_f(Fr1,Re1,Vp(d)),0)

    # Solve equation of motion
    y0 = [z0, V1[0]]
    
    # odefun = @(t, y) particle_ode(t, y, rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0, trec)
    # [t, y] = ode15s(odefun, [0 tend], y0, options)
    sol = solve_ivp(particle_ode, [0, tend], y0, args=(rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0, trec))
    
    t = sol.t
    zp = sol.y[0]
    V = sol.y[1]
        
    return t, zp, V


def particle_ode(t, y, rho1, rho, nu , Cd, g, rhop, d, zu, zl, Vc0, trec):
    
    global tzl
                             
    zp   = y[0]                # particle position      [m]
    V    = y[1]                # particle velocity      [m/s] 
    Cam = 0.5                  # added mass coefficient [-]
    # global tzl                  # time at which the particle reaches zl

    # determine Vc
    if (zp <= zl):
        Vc = Vc0
        tzl = t
    else:
        Vc = Vc0 * np.exp(-(t-tzl)/trec)

    # stratification force
    FS = (rho(zp) - rho1)*Vc*g

    # force balance on the particle 
    dzpdt = V
    dVdt  = ( (rhop-rho(zp)) * Vp(d) * g \
            - 0.5 *Cd(V*d/nu(zp))* rho(zp) * Ap(d) * np.abs(V) * V \
            - FS \
            ) / (rhop*Vp(d) + Cam*rho(zp)*Vp(d))

    dydt = [dzpdt, dVdt]

    return dydt
