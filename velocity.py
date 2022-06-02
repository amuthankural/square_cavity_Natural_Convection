def velocity_init(domain,m,n):
    pass

def velocity_bound(domain,m,n):
    for j in range(n):
        domain[0][j]      = 41
        domain[m-1][j]    = 23

    for i in range(m):
        domain[i][0]      = 12
        domain[i][n-1]    = 34

    return(domain)


def u_velocity(strm,u,dY,m,n):    
    for i in range(1,m-1):
        for j in range(1,n-1):
            u[i][j] = (strm[i][j+1]-strm[i][j-1])/(2*dY)
    return u

def v_velocity(strm,v,dX,m,n):    
    for i in range(1,m-1):
        for j in range(1,n-1):
            v[i][j] = -(strm[i+1][j]-strm[i-1][j])/(2*dX)
    return v