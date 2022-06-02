import numpy as np
import logging
from pathlib import Path

#Output folder setup
output = Path('./output/log').expanduser()
output.mkdir(parents=True, exist_ok=True)


strm_log = logging.getLogger(__name__)
strm_log.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('./output/log/stream.log',mode='w')

file_handler.setFormatter(formatter)
strm_log.addHandler(file_handler)


def strm_init(domain,m,n):
    for j in range(n):
        domain[0][j]      = 0
        domain[m-1][j]    = 0
    for i in range(m):
        domain[i][0]      = 0
        domain[i][n-1]    = 0  
    return(domain)


def strm(strm_o,vort,m,n,dX,dY,div,iter):
    strm_log.debug("Iteration No.: {}--------A".format(iter))
    strm_log.debug("Temperature at bouondary calculation entry: \n{}".format(strm_o))
    strm_calc   = np.copy(strm_o)
    for i in range(1,m-1):
        for j in range(1,n-1):
            strm_i_sum  = (strm_o[i+1,j]+strm_o[i-1,j])/(dX*dX)
            strm_j_sum  = (strm_o[i,j+1]+strm_o[i,j-1])/(dY*dY)
            strm_calc[i][j]   = ( strm_i_sum + strm_j_sum + vort[i,j] ) / (div)
    strm_log.debug("Vorticity at recalculation exit: \n{}".format(strm_calc))
    strm_log.debug("____________________________________________________")
    return(strm_calc)


def strm_ur(strm_o,strm_calc,m,n,r):
    strm_n      = np.copy(strm_calc)
    for i in range(1,m-1):
        for j in range(1,n-1):
            strm_n[i][j] = strm_o[i][j] + r*(strm_calc[i][j] - strm_o[i][j])
    return(strm_n)


def converge(strm_o,vort,m,n,dX,dY,div,iter):
    strm_residue    = np.zeros((m,n))
    strm_max        = int(np.amax(np.abs(strm_o)))
    if strm_max     == 0:
        strm_max    = 1
    strm_log.debug("Iteration level:{}\nStream function max value: \n{}".format(iter,strm_max))
    for i in range(1,m-1):
        for j in range(1,n-1):
            strm_i_sum  = (strm_o[i+1,j]+strm_o[i-1,j])/(dX*dX)
            strm_j_sum  = (strm_o[i,j+1]+strm_o[i,j-1])/(dY*dY)            
            strm_residue[i][j]   = ((( strm_i_sum + strm_j_sum + vort[i,j] ) / (div)) - strm_o[i][j])/strm_max
    strm_log.debug("Stream function residue domain: \n{}".format(strm_residue))
    return np.std(strm_residue)