import numpy as np
import logging
from pathlib import Path

#Output folder setup
output = Path('./output/log').expanduser()
output.mkdir(parents=True, exist_ok=True)

Nu_log = logging.getLogger(__name__)
Nu_log.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('./output/log/nusselt.log',mode='w')

file_handler.setFormatter(formatter)
Nu_log.addHandler(file_handler)


def nusselt(Nu,temp_calc,m,n,dX):
    NuH_mid = 0
    NuC_mid = 0
    for j in range(n):
        Nu[0,j]     = ( 1 - temp_calc[1,j])/ dX
        Nu[m-1,j]   = (temp_calc[m-2,j])/ dX
    
    Nu_log.debug('Nusselt number in domain:\n{}'.format(Nu))


    NuH_end         = (Nu[0,0]+Nu[0,n-1])/2
    NuC_end         = (Nu[m-1,0]+Nu[m-1,n-1])/2
    Nu_log.debug('Nusselt number on end nodes of hot side:\t{}'.format(NuH_end))
    Nu_log.debug('Nusselt number on end nodes of cold side:\t{}'.format(NuC_end))

    for j in range(1,n-1):
        NuH_mid     += (Nu[0,j])
        NuC_mid     += (Nu[m-1,j])
    
    NuH              = (NuH_end + NuH_mid)*dX
    NuC              = (NuC_end + NuC_mid)*dX

    Nu_log.debug('Nusselt number on hot side:\t{}'.format(NuH))
    Nu_log.debug('Nusselt number on cold side:\t{}'.format(NuC))
    
    return(Nu,[NuH,NuC])
    