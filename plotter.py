import os
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