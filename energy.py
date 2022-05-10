import numpy as np
import logging
from pathlib import Path

#Output folder setup
output = Path('./output/log').expanduser()
output.mkdir(parents=True, exist_ok=True)

en_log = logging.getLogger(__name__)
en_log.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('./output/log/energy.log',mode='w')

file_handler.setFormatter(formatter)
en_log.addHandler(file_handler)


def energy_init(temp,m,n,dX):
    for i in range(1,m-1):
        for j in range (1,n-1):
            temp[j, i] = 1 - j*dX
    return(temp)


def energy_bound(temp,m,n,iter):
    en_log.debug("Iteration No.: {}--------C".format(iter))
    en_log.debug("Temperature at bouondary calculation entry: \n{}".format(temp))
    for j in range(n):
        temp[0,j]      = 1
        temp[m-1,j]    = 0
    for i in range(m):
        temp[i,0]      = (4/3)*( temp[i,2]- (temp[i, 3]/4) )
        temp[i,n-1]    = (4/3)*( temp[i,j-1] - (temp[i, j-2]/4) )
    en_log.debug("Temperature at boundary calculation exit: \n{}".format(temp))
    en_log.debug("____________________________________________________")
    return(temp)


def energy_bound_ur(temp_o,temp_calc,m,n,r,iter):
    temp_n      = np.copy(temp_calc)
    en_log.debug("Iteration No.: {}--------D".format(iter))
    en_log.debug("Temperature at boundary UR calculation entry: \n{}".format(temp_calc))
    for i in range(m):
        temp_n[i,0]      = temp_o[i,0] + r*(temp_calc[i,0] - temp_o[i,0])
        temp_n[i,n-1]    = temp_o[i,n-1] + r*(temp_calc[i,n-1] - temp_o[i,n-1])
    for j in range(n):
        temp_n[0,j]      = temp_o[0,j] + r*(temp_calc[0,j] - temp_o[0,j])
        temp_n[m-1,j]    = temp_o[m-1,j] + r*(temp_calc[m-1,j] - temp_o[m-1,j])
    en_log.debug("Temperature at boundary UR calculation exit: \n{}".format(temp_n))
    en_log.debug("____________________________________________________")
    return(temp_n)


def energy(temp_o,strm,m,n,dX,dY,div,iter):
    en_log.debug("Iteration No.: {}--------A".format(iter))
    en_log.debug("Temperature at calculation entry: \n{}".format(temp_o))
    temp_calc   = np.copy(temp_o)
    mul         = (-1/(4*dX*dY))
    for i in range(1,m-1):
        for j in range (1,n-1):          
            strm_i_diff = (strm[i,j+1]-strm[i,j-1])
            strm_j_diff = (strm[i,j+1]-strm[i,j-1])
            temp_i_diff = (temp_o[i+1,j]-temp_o[i-1,j])
            temp_j_diff = (temp_o[i,j+1]-temp_o[i,j-1])
            temp_i_sum = (temp_o[i+1,j]+temp_o[i-1,j])
            temp_j_sum = (temp_o[i,j+1]+temp_o[i,j-1])
            temp_calc[i,j] = ( ( (mul*((strm_j_diff*temp_i_diff)-(strm_i_diff*temp_j_diff))) + temp_i_sum + temp_j_sum )/div )

    en_log.debug("Temperature at calculation exit: \n{}".format(temp_calc))
    en_log.debug("____________________________________________________")
    return temp_calc


def energy_ur(temp_o,temp_calc,m,n,r,iter):
    en_log.debug("Iteration No.: {}--------B".format(iter))
    en_log.debug("Temperature at UR calculation entry: \n{}".format(temp_calc))
    temp_n      = np.copy(temp_calc)
    for i in range(1,m-1):
        for j in range (1,n-1):
            temp_n[i,j] = temp_o[i,j] + r*(temp_calc[i,j] - temp_o[i,j])
    en_log.debug("Temperature at UR calculation exit: \n{}".format(temp_n))
    en_log.debug("____________________________________________________")
    return(temp_n)


def converge(temp_o,strm,m,n,dX,dY,div):
    temp_residue = np.zeros((m,n))
    mul         = (-1/(4*dX*dY))
    for i in range(1,m-1):
        for j in range (1,n-1):         
            strm_i_diff = (strm[i+1,j]-strm[i-1,j])
            strm_j_diff = (strm[i,j+1]-strm[i,j-1])
            temp_i_diff = (temp_o[i+1,j]-temp_o[i-1,j])
            temp_j_diff = (temp_o[i,j+1]-temp_o[i,j-1])
            temp_i_sum = (temp_o[i+1,j]+temp_o[i-1,j])/(dX*dX)
            temp_j_sum = (temp_o[i,j+1]+temp_o[i,j-1])/(dY*dY)
            temp_residue[i,j] = ( ( (mul*((strm_j_diff*temp_i_diff)-(strm_i_diff*temp_j_diff))) + temp_i_sum + temp_j_sum )/div )- temp_o[i,j]
    return np.std(temp_residue)