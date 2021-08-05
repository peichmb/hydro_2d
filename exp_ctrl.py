# This routine controls when the experiment being carried out must stop, and gives appropriate output.

import numpy as np
import params_common as par
from deriv import *
import grid_common as gr
import __main__

# As each experiment is different the module takes the values
#   it needs for each experiment, instead of being written
#   on the module's symbol table from other routines such as read_input.
first=True

# Variables for the 1d Ray Tracing Experiment:
k0=par.kx
x0=1.0
phi=(par.xpack0-gr.x0)/(gr.xf-gr.x0)*k0
f=0

def cs(x):
    return np.sqrt(par.gam*par.p00/(par.rho00*np.exp((x-par.xpack0)/3.)))
    
def ddx_cs(x):
    return -1./6./(par.gam*par.p00*par.rho00*np.exp((x-par.xpack0)/3.))

def exp_ctrl(filename):    
    global first
    stop=False
    if filename == 'exps/ray1d.dat':
        global k0,x0,phi,f
        if first:
            first=False
            f=open('ray1dout.dat','w')
        else:          
            # Calculo magnitudes trazado de rayos 
            k0=k0-k0*ddx_cs(x0)*__main__.dt
            x0=x0+cs(x0)*__main__.dt
            phi=(x0-gr.x0)/(gr.xf-gr.x0)*k0
            
        f.write(str(__main__.time)+' '+str(x0)+' '+str(k0)+'\n')
        #print(str(__main__.time)+' '+str(x0)+' '+str(k0)) +' '+str(phi)
        
        
        
            
    
    return stop       
        
