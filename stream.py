import numpy as np


def strm_init(domain,m,n):
    for j in range(n):
        domain[0][j]      = 0
        domain[m-1][j]    = 0
    for i in range(m):
        domain[i][0]      = 0
        domain[i][n-1]    = 0
    return(domain)


def strm(strm_o,vort,m,n,dX,dY,div):
    strm_calc   = np.copy(strm_o)
    for i in range(m-1):
        for j in range(n-1):
            strm_calc[i][j]   = (   (strm_o[i+1][j]+strm_o[i-1][j])/(dX*dX) +\
                (strm_o[i][j+1]+strm_o[i][j-1])/(dY*dY) +\
                    vort[i][j]  )/(div)
    return(strm_calc)


def strm_ur(strm_o,strm_calc,m,n,r):
    strm_n      = np.copy(strm_o)
    for i in range(m-1):
        for j in range(n-1):
            strm_n[i][j] = strm_o[i][j] + r*(strm_calc[i][j] - strm_o[i][j])
    return(strm_n)