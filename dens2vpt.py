# Calculate primitive variables from densities:

import params_common as par

def dens2vpt(rho,momx,momz,energ):

    vx = momx/rho
    vz = momz/rho
    pres = (par.gam-1.)*(energ-(momx*momx+momz*momz)/2./rho)
    
    return vx,vz,pres
