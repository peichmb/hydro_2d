# Calculate fluxes

from dens2vpt import *
import numpy as np

def fluxes(rho,momx,momz,energ):

    vx,vz,pres=dens2vpt(rho,momx,momz,energ)

    # Density flux
    flm_x=momx
    flm_z=momz
    
    # Momentum flux
    flcx_x=momx*momx/rho + pres
    flcx_z=momx*momz/rho
    
    flcz_x=flcx_z
    flcz_z=momz*momz/rho + pres
    
    # Energy flux
    fle_x=(energ+pres)*vx
    fle_z=(energ+pres)*vz
    
    return (flm_x,flm_z), (flcx_x,flcx_z), (flcz_x,flcz_z), (fle_x,fle_z)
