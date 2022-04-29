import numpy as np
import stream
import vorticity
import energy
import velocity
import temperature
import tester


#Input Parameters
m       = 10
n       = 10
dX      = 1/m
dY      = 1/n
Pr      = 1
Ra      = 10^3
phi     = 90
r       = 0.83

#Boundary Conditions
bound_u_vel     = [0,0,0,0]
bound_v_vel     = [0,0,0,0]
bound_strm      = [0,0,0,0]
bound_temp      = [0,0,0,0]
bound_vort      = [0,0,0,0]


def discretize(m,n):
    return(np.zeros((m,n)))

test        = discretize(m,n)
u_vel       = discretize(m,n)
v_vel       = discretize(m,n)
strm        = discretize(m,n)
temp        = discretize(m,n)
vort        = discretize(m,n)

#Calculated Values
div         = ((2/(dX*dX))+(2/(dY*dY)))


#vort = vorticity.vort(vort,strm,temp,m,n,dX,dY,div,Pr,Ra,phi,r)


temp = energy.energy_init(temp,m,n,dX)
print(temp)

temp = energy.energy_bound(temp,m,n)
print(temp)

