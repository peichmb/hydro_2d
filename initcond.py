# coding=latin-1

# ----------------------------------------------------------------------
# ROUTINE INITCOND
#
#    PURPOSE:  Calculate the arrays of density, velocity and pressure 
#                 at time t=0. 
#    INPUT ARGUMENTS: 
#              condition:   char variable with the choice of initial condition
#                 shape :   char variable with some subsidiary choice
#    MODULES: Note that some necessary input is passed via modules
#                (like the grid settings, problem parameters...)
#    OUTPUT:  the arrays vx0, vz0, rho0, pres0
# ----------------------------------------------------------------------

import numpy as np
import grid_common as gr
import params_common as par

def initcond():
    
    # SOUND WAVES
    if par.inittype == 'sound wave':
        
        # SINE WAVE PROPAGATING IN x DIRECTION
        if par.shape == 'sine-x':      
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.empty([gr.npx,gr.npz],dtype=float)
            vz0=np.zeros([gr.npx,gr.npz],dtype=float)
            h = np.cos(2.*np.pi*((gr.xx-gr.x0)/(gr.xf-gr.x0)))         
            for i in range(gr.npz):
                pres0[:,i]=par.pres00*(1.+par.amp*par.gam*h)
                rho0[:,i]=par.rho00*(1.+par.amp*h)
                vx0[:,i]=par.vx00+par.cs00*par.amp*h

        # SINE WAVE PROPAGATING IN z DIRECTION
        if par.shape == 'sine-z':      
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.zeros([gr.npx,gr.npz],dtype=float)
            vz0=np.empty([gr.npx,gr.npz],dtype=float)
            h = np.cos(2.*np.pi*((gr.zz-gr.z0)/(gr.zf-gr.z0)))
            for i in range(gr.npx):
                pres0[i,:]=par.pres00*(1.+par.amp*par.gam*h)
                rho0[i,:]=par.rho00*(1.+par.amp*h)
                vz0[i,:]=par.vz00+par.cs00*par.amp*h

        # SINE WAVE PROPAGATING IN A GENERAL DIRECTION
        if par.shape == 'sine':      
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.zeros([gr.npx,gr.npz],dtype=float)
            vz0=np.zeros([gr.npx,gr.npz],dtype=float)
            h=np.empty([gr.npx,gr.npz],dtype=float)           
            for i in range(gr.npx):
                for j in range(gr.npz):
                    h[i,j]=np.sin( par.kx*(gr.xx[i]-gr.x0)+par.kz*(gr.zz[j]-gr.z0) )             
            pres0=par.pres00*(1.+par.amp*par.gam*h)
            rho0=par.rho00*(1.+par.amp*h)
            if par.kz != 0. and par.kx != 0.:
                alpha=np.arctan(par.kz/par.kx)
                vx0=par.vx00+par.cs00*np.cos(alpha)*par.amp*h
                vz0=par.vz00+par.cs00*np.sin(alpha)*par.amp*h
            elif par.kz == 0.:
                vx0=par.vx00+par.cs00*par.amp*h
            elif par.kx == 0.:
                vz0=par.vz00+par.cs00*par.amp*h
        
        # NICE RANDOM WAVE
        if par.shape == 'nice-xz':      
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.empty([gr.npx,gr.npz],dtype=float)
            vz0=np.empty([gr.npx,gr.npz],dtype=float)
            h = np.cos(2.*np.pi*((gr.xx-gr.x0)/(gr.xf-gr.x0)))
            for i in range(gr.npz):
                pres0[:,i]=par.pres00*(1.+par.amp*par.gam*h)
                rho0[:,i]=par.rho00*(1.+par.amp*h)
                vx0[:,i]=par.vx00+par.cs00/np.sqrt(2)*par.amp*h
                vz0[i,:]=par.vz00+par.cs00/np.sqrt(2)*par.amp*h  
                
    # SOUND PULSES, WAVE PACKETS                    
    if par.inittype == 'pulse':
    
        # GAUSSIAN PACKET #1
        if par.shape == 'gaussian1':
            W=gr.zlen/12.
            h=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.zeros([gr.npx,gr.npz],dtype=float)
            vz0=np.zeros([gr.npx,gr.npz],dtype=float)            
            for i in range(gr.npx):
                for j in range(gr.npz):
                    h[i,j]=np.exp(-((gr.xx[i]-gr.xmid)**2+(gr.zz[j]-gr.zmid)**2)/W**2)* \
                           np.sin(par.kx*(gr.xx[i]-gr.x0)+par.kz*(gr.zz[j]-gr.z0))
            if par.kz != 0. and par.kx != 0.:
                alpha=np.arctan(par.kz/par.kx)
                vx0=par.vx00+par.cs00*np.cos(alpha)*par.amp*h
                vz0=par.vz00+par.cs00*np.sin(alpha)*par.amp*h
            elif par.kz == 0.:
                vx0=par.vx00+par.cs00*par.amp*h
            elif par.kx == 0.:
                vz0=par.vz00+par.cs00*par.amp*h
            pres0=par.pres00*(1.+par.amp*par.gam*h)
            rho0=par.rho00*(1.+par.amp*h)
            
        # GAUSSIAN PACKET #2
        #   Plane wave modulated by gaussian in x, non uniform sound speed
        if par.shape == 'gaussian2':
            par.x0=gr.xmid/2.
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.empty([gr.npx,gr.npz],dtype=float)
            vz0=np.zeros([gr.npx,gr.npz],dtype=float)
            W=0.1
            h = np.cos(par.kx*((gr.xx-gr.x0)/(gr.xf-gr.x0)))*np.exp(-(gr.xx-par.x0)**2/W**2)        
            for i in range(gr.npz):
                pres0[:,i]=par.pres00*(1.+par.amp*par.gam*h)
                rho0[:,i]=par.rho00*(1.+par.amp*h)
                vx0[:,i]=par.vx00+par.cs00*par.amp*h
        if par.shape == 'gaussian3':
            par.xpack0=gr.xmid/2.
            pres0=np.empty([gr.npx,gr.npz],dtype=float)
            rho0=np.empty([gr.npx,gr.npz],dtype=float)
            vx0=np.empty([gr.npx,gr.npz],dtype=float)
            vz0=np.zeros([gr.npx,gr.npz],dtype=float)
            W=0.1
            h = np.cos(par.kx*((gr.xx-gr.x0)/(gr.xf-gr.x0)))*np.exp(-(gr.xx-par.xpack0)**2/W**2)
            for i in range(gr.npz):
                pres0[:,i]=par.pres00*(1.+par.amp*par.gam*h)
                rho0[:,i]=par.rho00*np.exp((gr.xx-par.xpack0)/3.)+par.rho00*par.amp*h
                vx0[:,i]=par.vx00+par.cs00*par.amp*h
            
    # OTHER INITIAL CONDITIONS YET TO IMPLEMENT        
    if par.shape == 'other':
        exit('routine initcond: no other condition is implemented yet')

    return vx0,vz0,rho0,pres0
