import os
from re import L
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt

def plotter(N,domain,m,n,title):

    xlist = np.linspace(0.0, 1.0, m)
    ylist = np.linspace(0.0, 1.0, n)

    X, Y = np.meshgrid(xlist, ylist)
    Z = domain

    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z, N, cmap="viridis")
    fig.colorbar(cp) # Add a colorbar to a plot

    ax.set_title(title)
    ax.set_xlabel('x (unit)')
    ax.set_ylabel('y (unit)')


    title = title+'.png'
    fname = os.path.join('./output/img/',title)


    plt.savefig(fname, bbox_inches='tight')

def vel_plot(domain, axis, m, n, title):
    if m%2 != 0 and n%2 != 0:
        m += 1
        n += 1

    if axis == 'x':
        l=int(m/2)
        xlist = np.linspace(0.0, 1.0, m-1)
        ylist = domain[l,0:n-1]
    elif axis == 'y':
        b=int(n/2)
        print("b:{}\nm:{}\nn:{}\nDomain:\n{}".format(b,m,m,domain))
        xlist = domain[0:m-1,b]
        ylist = np.linspace(0.0, 1.0, n-1)
        print("xlist:\n{}\nylist:\n{}".format(xlist,ylist))

    fig,ax=plt.subplots(1,1)
    sc = plt.scatter(xlist, ylist, c ="blue")

    ax.set_title(title)
    ax.set_xlabel('x (unit)')
    ax.set_ylabel('y (unit)')

    title = title+'.png'
    fname = os.path.join('./output/img/',title)
    plt.savefig(fname, bbox_inches='tight')

    