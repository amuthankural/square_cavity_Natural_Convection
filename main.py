import sys
import logging
import json
import numpy as np
from pathlib import Path


import stream
import vorticity
import energy
import nusselt as nu
import plotter
import velocity


np.set_printoptions(precision=2)
np.set_printoptions(threshold=sys.maxsize)


#Output folder setup
data = Path('./output/data').expanduser()
data.mkdir(parents=True, exist_ok=True)
img = Path('./output/img').expanduser()
img.mkdir(parents=True, exist_ok=True)


#File logger setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('./output/Natural_convection.log',mode='w')
file_handler.setFormatter(formatter)    
logger.addHandler(file_handler)


#Domain discretization function
def discretize(m,n):
    domain = np.zeros((m,n))
    return domain



def prog():
    with open('./input.json','r') as f:
        input           = json.load(f)


    #'While' loop flags
    i           = 1                     #Iteration number
    status      = False                 #Convergence status


    #Input Parameters
    m           = input["m"]            # Divisions in X - axis
    n           = input["n"]            # Divisions in Y - axis
    Pr          = input["Pr"]           # Prandtl number
    Ra          = input["Ra"]           # Rayleigh number
    phi         = input["phi"]          # Angle of orientation
    r           = input["r"]            # Relaxation factor for inner nodes
    rb          = input["rb"]           # Relaxation factor for boundary nodes
    convergence = input["convergence"]  # Convergence criteria
    N           = input["resolution"]   # Contour Plot resolution


    print("Data imported")
    logger.info('Data imported from JSON file\n---------------------------------------------------------------------')
    logger.info('\nDomain data:\nX_Div:\t\t{}\nY_Div:\t\t{}\nPrandtl:\t{}\nRayleigh:\t{}\nAngle:\t\t{}\nRelaxation:\t{}\nConvergence:\t{}\n'\
        .format(m,n,Pr,Ra,phi,r,convergence))
    logger.info('\n---------------------------------------------------------------------')


    #Calculated Values
    dX          = 1/m
    dY          = 1/n
    dX2         = dX*dX
    dY2         = dY*dY
    div         = ((2/dX2)+(2/dY2))



    #Domain definition
    vort_o      = discretize(m,n)
    strm_o      = discretize(m,n)
    temp_o      = discretize(m,n)
    vort_calc   = discretize(m,n)
    strm_calc   = discretize(m,n)
    temp_calc   = discretize(m,n)
    nslt        = discretize(m,n)
    u_vel       = discretize(m,n)
    v_vel       = discretize(m,n)


    #Residue definition
    residue      = []
    vort_residue = []
    strm_residue = []
    temp_residue = []


    #STEP - I: Domain initial condition
    vort_o      = vorticity.vort_init(vort_o,m,n)
    strm_o      = stream.strm_init(strm_o,m,n)
    temp_o      = energy.energy_init(temp_o,m,n,dX)
    logger.debug('Vorticity - Initial Domain: \n{}'.format(vort_o))
    logger.debug('Stream function - Initial Domain: \n{}'.format(strm_o))
    logger.debug('Temperature - Initial Domain: \n{}'.format(temp_o))


    #FDM main loop
    while (status != True):
        logger.debug("\n\n__________________________________________________________________________________________")
        logger.debug("\nIteration no: \t{}".format(i))


        #STEP - II: Applying boundary conditions:
        vort_o      = vorticity.vort_bound(vort_o,strm_o,m,n,dX,dY,i)
        logger.debug('Vorticity - First Bounded Domain: \n{}'.format(vort_o))
        strm_o      = stream.strm_init(strm_o,m,n)
        logger.debug('Stream function - First Bounded Domain: \n{}'.format(strm_o))
        temp_o      = energy.energy_bound(temp_o,m,n,i)
        logger.debug('Temperature - First Bounded Domain: \n{}'.format(temp_o))

        
        #STEP - III: Vorticity domain is calculated:
        vort_calc   = vorticity.vort(vort_o,strm_o,temp_o,m,n,dX,dY,div,Pr,Ra,phi)
        logger.debug('Vorticity - Calculated Domain: \n{}'.format(vort_calc))
        #STEP - IV: Under relaxation for vorticity domain:
        vort_calc   = vorticity.vort_ur(vort_o,vort_calc,m,n,r)
        logger.debug('Vorticity - Calculated UR Domain: \n{}'.format(vort_calc))


        #STEP - V: Stream function domain is calculated:
        strm_calc   = stream.strm(strm_o,vort_calc,m,n,dX,dY,div,i)
        logger.debug('Stream function - Calculated Domain: \n{}'.format(strm_o))
        #STEP - VI: Under relaxation for stream function
        strm_calc   = stream.strm_ur(strm_o,strm_calc,m,n,r)
        logger.debug('Stream function - Calculated UR Domain: \n{}'.format(strm_o))


        #STEP - VII: Vorticity boundary condition is calculated for new stream function:
        vort_calc   = vorticity.vort_bound(vort_calc,strm_calc,m,n,dX,dY,i)
        logger.debug('Vorticity - Second Bounded Domain: \n{}'.format(vort_o))
        #STEP - VIII: Under relaxation for vorticity domains boundary:
        vort_calc   = vorticity.vort_ur(vort_o,vort_calc,m,n,rb)
        logger.debug('Vorticity - Second Bounded UR Domain: \n{}'.format(vort_o))


        #STEP - IX: Energy domain is calculated:
        temp_calc   = energy.energy(temp_o,strm_calc,m,n,dX,dY,div,i)
        logger.debug('Temperature - Calculated Domain: \n{}'.format(temp_o))
        #STEP - X: Under relaxation for energy domain:
        temp_calc   = energy.energy_ur(temp_o,temp_calc,m,n,r,i)
        logger.debug('Temperature - Calculated UR Domain: \n{}'.format(temp_o))


        #STEP - XI: Energy boundary condition is applied:
        temp_calc   = energy.energy_bound(temp_calc,m,n,i)
        logger.debug('Temperature - Second Bounded Domain: \n{}'.format(temp_o))
        """
        #STEP - XI(b): Under relaxation for temperature boundary:
        temp_calc   = energy.energy_bound_ur(temp_o,temp_calc,m,n,r,i)
        logger.debug('Temperature - Second Bounded UR Domain: \n{}'.format(temp_o))"""


        #Residuals for Vorticity, Stream function & Temperature are calculated
        vort_residue.append(vorticity.converge(vort_calc,strm_calc,temp_calc,m,n,dX,dY,div,Pr,Ra,phi,i))
        residue.append(vort_residue)
        strm_residue.append(stream.converge(strm_calc,vort_calc,m,n,dX,dY,div,i))
        residue.append(strm_residue)
        temp_residue.append(energy.converge(temp_calc,strm_calc,m,n,dX,dY,div,i))
        residue.append(temp_residue)


        #Relative residue and max among the above three are calculated:
        rel_residue = np.std([vort_residue[-1],strm_residue[-1],temp_residue[-1]])
        max_residue = max([vort_residue[-1],strm_residue[-1],temp_residue[-1]])


        #Check for logging data:
        if i%10 == 0:
            logger.info("\nIteration no: \t\t\t{}".format(i))
            logger.info('Vorticity residual: \t\t{}'.format(vort_residue[-1]))
            logger.info('Stream function residual: \t{}'.format(strm_residue[-1]))
            logger.info('Temperature residual: \t\t{}'.format(temp_residue[-1]))
            logger.debug("Relative residue: \t\t{}".format(rel_residue))
            logger.debug("Max residue: \t\t\t{}".format(max_residue))
            logger.info("Status: \t\t\t{}".format(status))
            print("Iteration no.:\t{}".format(i))
        

        #COnvergence condition check:
        if  max_residue <= convergence or rel_residue <= convergence or i > 1000:
            status = True


        vort_o  = vort_calc
        strm_o  = strm_calc
        temp_o  = temp_calc
        i      += 1


    #Nusselt number is calculated along the domain boundary:
    nusselt = nu.nusselt(nslt,temp_calc,m,n,dX)

    #Velocity components calculated:
    u_vel   = velocity.u_velocity(strm_calc,u_vel,dY,m,n)
    v_vel   = velocity.v_velocity(strm_calc,v_vel,dX,m,n)


    #Numpy matrix is transposed to obtain correctly alligned output data
    vort_plot = vort_calc.T
    strm_plot = strm_calc.T
    temp_plot = temp_calc.T
    nslt_plot = nslt.T
    u_vel_plot= u_vel.T
    v_vel_plot= -(v_vel.T)



    #Domain files are saved as numpy dataframes:    
    np.save(data/'Vorticity.npy',vort_calc)
    np.save(data/'StreamFunction.npy',strm_calc)
    np.save(data/'Temperature.npy',temp_calc)
    np.save(data/'Nusselt.npy',nusselt[0])
    np.save(data/'u_velocity.npy',u_vel)
    np.save(data/'v_velocity.npy',v_vel)
    logger.info('Data saved\n---------------------------------------------------------------------')
    print("Data saved")


    #Plotting the domain data into a contour plot:
    plotter.plotter(N,strm_plot,m,n,"Contour of Stream function")
    plotter.plotter(N,vort_plot,m,n,"Contour of Vorticity")
    plotter.plotter(N,temp_plot,m,n,"Contour of Temperature")
    plotter.plotter(N,nslt_plot,m,n,"Nusselt number on hot and cold side")
    plotter.vel_plot(u_vel_plot,'y',m,n,"U-velocity")
    plotter.vel_plot(v_vel_plot,'x',m,n,"V-velocity")
    logger.info('Images saved\n---------------------------------------------------------------------')
    print("Images saved")
    print("Log file is saved as Natural_convection.log")

    return(i,vort_calc,strm_calc,temp_calc,vort_residue,strm_residue,temp_residue,rel_residue,nusselt,u_vel,v_vel)


output = prog()



#Logging of final values for domains and residuals
logger.info('No. of Iterations: \t\t\t{}'.format(output[0]))
logger.info('\nRelative residue: \t\t\t{}'.format(output[7]))
logger.info('\nStream function Residue: \t\t{}'.format(output[5]))
logger.info('\nVorticity Residue: \t\t\t{}'.format(output[4]))
logger.info('\nTemperature Residue: \t\t\t{}'.format(output[6]))
logger.info('\nVorticity: \n{}\nStream_Function:\n{}\nTemperature:\n{}\nNusselt No:{}\nu_Velocity:\n{}\nv_velocity:\n{}'.format(output[1],output[2],output[3],output[8],output[9],output[10]))
print('No. of Iterations: ',output[0])
print('\nRelative residue:\n',output[7])
print('\nVorticity:\n',output[1],'\nStream_Function:\n',output[2],'\nTemperature:\n',output[3])
print('Nusselt domain:\n',output[8])
