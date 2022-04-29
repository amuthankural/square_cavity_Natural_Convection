def pressure_init(domain,m,n):
    pass

def pressure_bound(domain,m,n):
    for j in range(n):
        domain[0][j]      = 41
        domain[m-1][j]    = 23

    for i in range(m):
        domain[i][0]      = 12
        domain[i][n-1]    = 34

    return(domain)