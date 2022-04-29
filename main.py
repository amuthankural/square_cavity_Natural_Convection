import numpy as np
import stream
import vorticity
import energy


#Input Parameters
m       = 10
n       = 10
dX      = 1/m
dY      = 1/n
Pr      = 1
Ra      = 1000
phi     = 90
r       = 0.83


#Calculated Values
div         = ((2/(dX*dX))+(2/(dY*dY)))


def discretize(m,n):
    return(np.zeros((m,n)))


u_vel       = discretize(m,n)
v_vel       = discretize(m,n)
vort_o      = discretize(m,n)
strm_o      = discretize(m,n)
temp_o      = discretize(m,n)
vort_calc   = discretize(m,n)
strm_calc   = discretize(m,n)
temp_calc   = discretize(m,n)

vort_o      = vorticity.vort_init(vort_o,m,n)
strm_o      = stream.strm_init(strm_o,m,n)
temp_o      = energy.energy_init(temp_o,m,n,dX)


vort_o      = vorticity.vort_bound(vort_o,strm_o,m,n,dX,dY)
temp_o      = energy.energy_bound(temp_o,m,n)


vort_calc   = vorticity.vort(vort_o,strm_o,temp_o,m,n,dX,dY,div,Pr,Ra,phi)
vort_calc   = vorticity.vort_ur(vort_o,vort_calc,m,n,r)

strm_calc   = stream.strm(strm_o,vort_calc,m,n,dX,dY,div)
strm_calc   = stream.strm_ur(strm_o,strm_calc,m,n,r)

vort_calc   = vorticity.vort_bound(vort_calc,strm_calc,m,n,dX,dY)
vort_calc   = vorticity.vort_ur(vort_o,vort_calc,m,n,r)

temp_calc   = energy.energy(temp_o,strm_calc,m,n,dX,dY,div)
temp_calc   = energy.energy_ur(temp_o,temp_calc,m,n,r)


temp_calc   = energy.energy_bound(temp_calc,m,n)
temp_calc   = energy.energy_bound_ur(temp_o,temp_calc,m,n,r)


print('Vorticity:\n',vort_calc,'\nStream_Function:\n',strm_calc,'\nTemperature:\n',temp_calc)