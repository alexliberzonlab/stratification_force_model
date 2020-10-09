from scipy.optimize import least_squares
import numpy as np


def Cd(Re): return 0.4 + 24.0/Re + 6.0/(1+np.sqrt(Re))  
def Vp(d): return np.pi*d**3/6
def Ap(d): return np.pi*d**2/4

def settling_velocity(rhop=1.05, rhof=1.0 , g=9.81, d=100e-6, nu=1e-6):
    """ Settling velocity 
    
    Inputs:
        rhop : particle density, (kg/m^3)
        rhof : fluid density, (kg/m^3)
        g : gravity acceleration, default = 9.81 ( m/s^2 )
        d : particle diameter, (m)
        nu : kinematic viscosity, (m^2/s), water is 1e-6
    
    
    Outputs:
        ws : 
    
    """
    
    def f(wp): return 1e12*((rhop-rhof)*Vp(d)*g - 0.5*Cd(wp*d/nu)*rhof*Ap(d)*np.abs(wp)*wp)


    # Calculate the settling speed from the force balance in the steady-state conditions
    sol = least_squares(f, .5, bounds=[1e-11,1],ftol=1e-15,xtol=1e-8)

    return sol