import numpy as np


def energy_init(temp,m,n,dX):
    for i in range(1,m-1):
        for j in range (1,n-1):
            temp[j, i] = 1 - j*dX
    return(temp)


def energy_bound(temp,m,n):
    for j in range(n):
        temp[0][j]      = 1
        temp[m-1][j]    = 0
    for i in range(m):
        temp[i][0]      = (4/3)*( temp[i,2]- (temp[i, 3]/4) )
        temp[i][n-1]    = (4/3)*( temp[i,j-1] - (temp[i, j-2]/4) )
    return(temp)


def energy_bound_ur(temp_o,temp_calc,m,n,r):
    temp_n      = np.copy(temp_o)   
    for i in range(m):
        temp_n[i][0]      = temp_o[i][0] + r*(temp_calc[i][0] - temp_o[i][0])
        temp_n[i][n-1]    = temp_o[i][n-1] + r*(temp_calc[i][n-1] - temp_o[i][n-1])
    for j in range(n):
        temp_n[0][j]      = temp_o[0][j] + r*(temp_calc[0][j] - temp_o[0][j])
        temp_n[m-1][j]    = temp_o[m-1][j] + r*(temp_calc[m-1][j] - temp_o[m-1][j])
    return(temp_n)


def energy(temp_o,strm,m,n,dX,dY,div):
    temp_calc = temp_o.copy()
    for i in range(1,m-1):
        for j in range (1,n-1):
            temp_o[i,j] = ( ( (-1/(4*dX*dY)) * \
                ( (strm[i,j+1]-strm[i,j-1])*(temp_o[i+1,j]-temp_o[i-1,j])-\
                    (strm[i,j+1]-strm[i,j-1])*(temp_o[i+1,j]-temp_o[i-1,j])))+\
                        ((temp_o[i+1,j]+temp_o[i-1,j])/(dX*dX)) +\
                            ((temp_o[i,j+1]+temp_o[i,j-1])/(dY*dY)))/(div)
    return temp_calc


def energy_ur(temp_o,temp_calc,m,n,r):
    temp_n      = np.copy(temp_o)
    for i in range(1,m-1):
        for j in range (1,n-1):
            temp_n[i][j] = temp_o[i][j] + r*(temp_calc[i][j] - temp_o[i][j])
    return(temp_n)