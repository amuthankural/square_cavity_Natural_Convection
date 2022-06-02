import numpy as np
import math
import logging
from pathlib import Path

#Output folder setup
output = Path('./output/log').expanduser()
output.mkdir(parents=True, exist_ok=True)


vort_log = logging.getLogger(__name__)
vort_log.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('./output/log/vorticity.log',mode='w')

file_handler.setFormatter(formatter)
vort_log.addHandler(file_handler)

def vort_init(domain,m,n):
    for j in range(n):
        domain[0][j]      = 0
        domain[m-1][j]    = 0
    for i in range(m):
        domain[i][0]      = 0
        domain[i][n-1]    = 0
    return(domain)


def vort_bound(vort,strm,m,n,dX,dY,iter):
    vort_log.debug("Iteration No.: {}".format(iter))
    vort_log.debug("Vorticity at recalculation entry: \n{}".format(vort))
    for i in range(m):
        vort[i][0]      = (-2*strm[i][1])/(dY*dY)
        vort[i][n-1]    = (-2*strm[i][n-2])/(dY*dY)
    for j in range(n):
        vort[0][j]      = (-2*strm[1][j])/(dX*dX)
        vort[m-1][j]    = (-2*strm[m-2][j])/(dX*dX)
    vort_log.debug("Vorticity at recalculation exit: \n{}".format(vort))
    vort_log.debug("____________________________________________________")
    return vort


def vort_bound_ur(vort_o,vort_calc,m,n,r):
    vort_n      = np.copy(vort_calc)
    for i in range(m):
        vort_n[i][0]      = vort_o[i][0] + r*(vort_calc[i][0] - vort_o[i][0])
        vort_n[i][n-1]    = vort_o[i][n-1] + r*(vort_calc[i][n-1] - vort_o[i][n-1])
    for j in range(n):
        vort_n[0][j]      = vort_o[0][j] + r*(vort_calc[0][j] - vort_o[0][j])
        vort_n[m-1][j]    = vort_o[m-1][j] + r*(vort_calc[m-1][j] - vort_o[m-1][j])
    return(vort_n)


def vort(vort_o,strm,temp,m,n,dX,dY,div,Pr,Ra,phi):
    vort_calc   = np.copy(vort_o)
    mul         = (-1/(4*dX*dY*Pr))

    for i in range(1,m-1):
        for j in range(1,n-1):
            strm_i_diff = (strm[i+1,j]-strm[i-1,j])
            strm_j_diff = (strm[i,j+1]-strm[i,j-1])
            temp_i_diff = ((temp[i+1,j]-temp[i-1,j])*math.sin(phi))/(2*dX)
            temp_j_diff = ((temp[i,j+1]-temp[i,j-1])*math.cos(phi))/(2*dY)
            vort_i_diff = (vort_o[i+1,j]-vort_o[i-1,j])
            vort_j_diff = (vort_o[i,j+1]-vort_o[i,j-1])
            vort_i_sum  = (vort_o[i+1,j]+vort_o[i-1,j])/(dX*dX)
            vort_j_sum  = (vort_o[i,j+1]+vort_o[i,j-1])/(dY*dY)
            vort_calc[i,j] = (( (mul*((strm_j_diff*vort_i_diff)-(strm_i_diff*vort_j_diff))) + vort_i_sum + vort_j_sum + Ra*(temp_i_diff - temp_j_diff) ))/(div)
    return(vort_calc)



def vort_ur(vort_o,vort_calc,m,n,r):
    vort_n      = np.copy(vort_calc)
    for i in range(1,m-1):
        for j in range (1,n-1):
            vort_n[i][j] = vort_o[i][j] + r*(vort_calc[i][j] - vort_o[i][j])
    return(vort_n)


def converge(vort_o,strm,temp,m,n,dX,dY,div,Pr,Ra,phi,iter):
    vort_residue   = np.zeros((m,n))
    mul         = (-1/(4*dX*dY*Pr))
    vort_max    = int(np.amax(np.abs(vort_o)))
    if vort_max == 0:
        vort_max = 1
    vort_log.debug("Vorticity max value: \n{}".format(vort_max))
    for i in range(1,m-1):
        for j in range (1,n-1):
            strm_i_diff = (strm[i+1,j]-strm[i-1,j])
            strm_j_diff = (strm[i,j+1]-strm[i,j-1])
            temp_i_diff = ((temp[i+1,j]-temp[i-1,j])*math.sin(phi))/(2*dX)
            temp_j_diff = ((temp[i,j+1]-temp[i,j-1])*math.cos(phi))/(2*dY)
            vort_i_diff = (vort_o[i+1,j]-vort_o[i-1,j])
            vort_j_diff = (vort_o[i,j+1]-vort_o[i,j-1])
            vort_i_sum  = (vort_o[i+1,j]+vort_o[i-1,j])/(dX*dX)
            vort_j_sum  = (vort_o[i,j+1]+vort_o[i,j-1])/(dY*dY)
            vort_residue[i,j] = ((( (mul*((strm_j_diff*vort_i_diff)-(strm_i_diff*vort_j_diff))) + vort_i_sum + vort_j_sum + Ra*(temp_i_diff - temp_j_diff) )/(div))\
                - vort_o[i,j])/vort_max
    vort_log.debug("Vorticity residue domain: \n{}".format(vort_residue))
    return np.std(vort_residue)