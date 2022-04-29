import numpy as np




arr = np.zeros((2,4)).T


def tester(domain,m,n):

    for j in range(n):
        domain[0, j]      = 41
        domain[m-1, j]    = 23

    for i in range(m):
        domain[i, 0]      = 12
        domain[i, n-1]    = 34


    return domain