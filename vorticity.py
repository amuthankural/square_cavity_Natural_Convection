import numpy as np
import math

def vort_init(domain,m,n):
    for j in range(n):
        domain[0][j]      = 0
        domain[m-1][j]    = 0
    for i in range(m):
        domain[i][0]      = 0
        domain[i][n-1]    = 0
    return(domain)


def vort_bound_ur(vort_o,vort_calc,m,n,r):
    vort_n      = np.copy(vort_o)
    for i in range(m):
        vort_n[i][0]      = vort_o[i][0] + r*(vort_calc[i][0] - vort_o[i][0])
        vort_n[i][n-1]    = vort_o[i][n-1] + r*(vort_calc[i][n-1] - vort_o[i][n-1])
    for j in range(n):
        vort_n[0][j]      = vort_o[0][j] + r*(vort_calc[0][j] - vort_o[0][j])
        vort_n[m-1][j]    = vort_o[m-1][j] + r*(vort_calc[m-1][j] - vort_o[m-1][j])
    return(vort_n)


def vort_bound(vort,strm,m,n,dX,dY):
    for i in range(m):
        vort[i][0]      = (-2*strm[i][2])/(dY*dY)
        vort[i][n-1]    = (-2*strm[i][n-2])/(dY*dY)
    for j in range(n):
        vort[0][j]      = (-2*strm[2][j])/(dX*dX)
        vort[m-1][j]    = (-2*strm[m-2][j])/(dX*dX)
        return vort


def vort(vort_o,strm,temp,m,n,dX,dY,div,Pr,Ra,phi):
    vort_calc   = np.copy(vort_o)
    for i in range(1,m-1):
        for j in range (1,n-1):
            vort_calc[i][j] = ( (-1/(4*dX*dY*Pr)) * ((strm[i][j+1]-strm[i][j-1])*(vort_o[i+1][j]-vort_o[i-1][j])\
                -(strm[i+1][j]-strm[i-1][j])*(vort_o[i][j+1]-vort_o[i][j-1]) )\
                    + ((vort_o[i+1][j]+vort_o[i-1][j])/(dX*dX))\
                        + ((vort_o[i][j+1]+vort_o[i][j-1])/(dY*dY))\
                            - Ra*( (temp[i+1][j]-temp[i-1][j])*(math.sin(phi)/(2*dX))\
                                - (temp[i][j+1]-temp[i][j-1])*(math.cos(phi)/(2*dY)) ))/(div) 
    return(vort_calc)



def vort_ur(vort_o,vort_calc,m,n,r):
    vort_n      = np.copy(vort_o)
    for i in range(1,m-1):
        for j in range (1,n-1):
            vort_n[i][j] = vort_o[i][j] + r*(vort_calc[i][j] - vort_o[i][j])
    return(vort_n)