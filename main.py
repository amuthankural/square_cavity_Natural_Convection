import numpy as np
import stream
import vorticity
import energy


#Input Parameters
m           = 100
n           = 100
dX          = 1/m
dY          = 1/n
Pr          = 1
Ra          = 1000
phi         = 90
r           = 0.83
convergence = 0.0001


#Calculated Values
div         = ((2/(dX*dX))+(2/(dY*dY)))


def discretize(m,n):
    return(np.zeros((m,n)))


def prog():
    i = 1
    vort_o      = discretize(m,n)
    strm_o      = discretize(m,n)
    temp_o      = discretize(m,n)
    vort_calc   = discretize(m,n)
    strm_calc   = discretize(m,n)
    temp_calc   = discretize(m,n)
    status      = False

    residue      = []
    vort_residue = []
    strm_residue = []
    temp_residue = []
    while (status != True):
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


        vort_residue.append(vorticity.converge(vort_calc,strm_calc,temp_calc,m,n,dX,dY,div,Pr,Ra,phi))
        residue.append(vort_residue)
        strm_residue.append(stream.converge(strm_calc,vort_calc,m,n,dX,dY,div))
        residue.append(strm_residue)
        temp_residue.append(energy.converge(temp_calc,strm_calc,m,n,dX,dY,div))
        residue.append(temp_residue)
        rel_residue = np.std(residue)

        if rel_residue <= convergence:
            status = True

        if i%10 == 0:
            print("\nIteration no: ",i,"\n")
        
        i += 1

    return(i,vort_calc,strm_calc,temp_calc,vort_residue,strm_residue,temp_residue,rel_residue)


output = prog()


print('No. of Iterations: ',output[0])
print('\nVorticity:\n',output[1],'\nStream_Function:\n',output[2],'\nTemperature:\n',output[3])
print('\nVorticity Residue:\n',output[4])
print('\nStream function Residue:\n',output[4])
print('\nTemperature Residue:\n',output[5])
print('\nRelative residue:\n',output[6])