# ----------------------------------------------------------------------
# ROUTINE UPDATE 
#
#    PURPOSE:  Calculate variables at timestep n+1/2 (half) or n+1 (full).
#
#    INPUT ARGUMENTS:  - the densities at timestep n, namely rho, momz, energ
#                      - the fluxes mflz, momflzz, energflz
#                      - the request of half or full timestep (/half, /full)
#
#    COMMON BLOCKS: Note that some necessary input is passed via common
#                       blocks (like the grid parameters and array zz, etc)
#
#    OUTPUT:  the densities at timestep n+1/2 (half) or at timestep n (full),
#                       namely rhon, momzn, energn
# ----------------------------------------------------------------------



from deriv import *
import grid_common as gr
import numpy as np



def update(dt,rho,momx,momz,energ,flm,flcx,flcz,fle,mode):

    

    # ----------------------
    # DERIVATIVES OF FLUXES 
    # ----------------------
    
    flm_x, flm_z = flm
    flcx_x, flcx_z = flcx
    flcz_x, flcz_z = flcz
    fle_x, fle_z = fle
    
    dx_flm_x = ddx(flm_x)
    dz_flm_z = ddz(flm_z)    
    dx_flcx_x = ddx(flcx_x)
    dz_flcx_z = ddz(flcx_z)
    dx_flcz_x = ddx(flcz_x)
    dz_flcz_z = ddz(flcz_z)
    dx_fle_x = ddx(fle_x)
    dz_fle_z = ddz(fle_z)    
        
    # -----------------------
    # UPDATE TO HALF TIMESTEP
    # -----------------------

    if mode=='half':                       
    
        # Promedios de las densidades:
        
        rhoprom=np.empty([gr.npx-1,gr.npz-1],dtype=float)
        momxprom=np.empty([gr.npx-1,gr.npz-1],dtype=float)
        momzprom=np.empty([gr.npx-1,gr.npz-1],dtype=float)
        energprom=np.empty([gr.npx-1,gr.npz-1],dtype=float)
        
        for i in range(gr.npx-1):
            for j in range(gr.npz-1):
                rhoprom[i,j]=(rho[i,j]+rho[i+1,j]+rho[i,j+1]+rho[i+1,j+1])/4.
                momxprom[i,j]=(momx[i,j]+momx[i+1,j]+momx[i,j+1]+momx[i+1,j+1])/4.
                momzprom[i,j]=(momz[i,j]+momz[i+1,j]+momz[i,j+1]+momz[i+1,j+1])/4.
                energprom[i,j]=(energ[i,j]+energ[i+1,j]+energ[i,j+1]+energ[i+1,j+1])/4.                   
        
        #print rhoprom #@@@@@@@@@@@@@@@@@@@@@                        
        
        # Derivadas de los flujos en los puntos (i-1/2,j-1/2)
        dx_flm_x=0.5*(dx_flm_x[:,1:]+dx_flm_x[:,:-1]) # Tienen que tener (npz-1)^2 componentes
        dz_flm_z=0.5*(dz_flm_z[1:,:]+dz_flm_z[:-1,:])        
        
        dx_flcx_x=0.5*(dx_flcx_x[:,1:]+dx_flcx_x[:,:-1]) # Tienen que tener (npz-1)^2 componentes
        dz_flcx_z=0.5*(dz_flcx_z[1:,:]+dz_flcx_z[:-1,:])  

        dx_flcz_x=0.5*(dx_flcz_x[:,1:]+dx_flcz_x[:,:-1]) # Tienen que tener (npz-1)^2 componentes
        dz_flcz_z=0.5*(dz_flcz_z[1:,:]+dz_flcz_z[:-1,:])  
                
        dx_fle_x=0.5*(dx_fle_x[:,1:]+dx_fle_x[:,:-1]) # Tienen que tener (npz-1)^2 componentes
        dz_fle_z=0.5*(dz_fle_z[1:,:]+dz_fle_z[:-1,:])  
                                      
        # Densidades en los puntos (i-1/2,j-1/2)                              
        rhon = rhoprom - dt/2.*(dx_flm_x+dz_flm_z) # Tienen que tener (npz-1)^2 componentes
        momxn = momxprom - dt/2.*(dx_flcx_x+dz_flcx_z) 
        momzn = momzprom - dt/2.*(dx_flcz_x+dz_flcz_z)         
        energn = energprom - dt/2.*(dx_fle_x+dz_fle_z)
        
        #print 'HALF: ',np.min(energn) 
        
    # -----------------------
    # UPDATE TO FULL TIMESTEP
    # -----------------------

    elif mode=='full':
    
        rhon=np.empty([gr.npx,gr.npz],dtype=float)
        momxn=np.empty([gr.npx,gr.npz],dtype=float)
        momzn=np.empty([gr.npx,gr.npz],dtype=float)
        energn=np.empty([gr.npx,gr.npz],dtype=float)
        
        # Derivadas en los puntos (i,j)
        dx_flm_x=0.5*(dx_flm_x[:,1:]+dx_flm_x[:,:-1]) # Tienen que tener (npz-2)^2 componentes
        dz_flm_z=0.5*(dz_flm_z[1:,:]+dz_flm_z[:-1,:])        
        
        dx_flcx_x=0.5*(dx_flcx_x[:,1:]+dx_flcx_x[:,:-1]) # Tienen que tener (npz-2)^2 componentes
        dz_flcx_z=0.5*(dz_flcx_z[1:,:]+dz_flcx_z[:-1,:])  

        dx_flcz_x=0.5*(dx_flcz_x[:,1:]+dx_flcz_x[:,:-1]) # Tienen que tener (npz-2)^2 componentes
        dz_flcz_z=0.5*(dz_flcz_z[1:,:]+dz_flcz_z[:-1,:])  
                
        dx_fle_x=0.5*(dx_fle_x[:,1:]+dx_fle_x[:,:-1]) # Tienen que tener (npz-2)^2 componentes
        dz_fle_z=0.5*(dz_fle_z[1:,:]+dz_fle_z[:-1,:])  
        
        
        rhon[1:-1,1:-1]=rho[1:-1,1:-1]-dt*(dx_flm_x+dz_flm_z)
        momxn[1:-1,1:-1]=momx[1:-1,1:-1]-dt*(dx_flcx_x+dz_flcx_z)
        momzn[1:-1,1:-1]=momz[1:-1,1:-1]-dt*(dx_flcz_x+dz_flcz_z)
        energn[1:-1,1:-1]=energ[1:-1,1:-1]-dt*(dx_fle_x+dz_fle_z)                        
        
    return rhon,momxn,momzn,energn

