# ----------------------------------------------------------------------
#
# FUNCTION CFL
#
#    PURPOSE: Calculate the timestep to guarantee numerical stability
#
#    INPUT ARGUMENTS: those necessary to calculate the sound speed, namely
#           rho, momz, pres
#    COMMON BLOCKS: Note that some necessary input is passed via common
#                       blocks (like gam or the grid zz)
#    OUTPUT:  the delta t 
# 
# ----------------------------------------------------------------------

import numpy as np
import grid_common as gr
import params_common as par

def cfl(rho,momx,momz,pres):
    cs = np.sqrt(par.gam*pres/rho)

    vchar1 = np.abs(momx/rho + cs)
    vchar2 = np.abs(momx/rho - cs)
    vchar3 = np.abs(momz/rho + cs)
    vchar4 = np.abs(momz/rho - cs) 
    
    maxi=np.maximum(vchar1,vchar2)
    maxi=np.maximum(maxi,vchar3)
    maxi=np.maximum(maxi,vchar4)    

    dt = par.fcfl*np.min(gr.dx,gr.dz) / np.max(maxi)

    return dt

