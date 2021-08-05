# Calculate derivatives

import grid_common as gr
import numpy as np

def ddx(var):
    lenx=len(var[:,0])
    lenz=len(var[0,:])
    dx_var = np.empty([lenx-1,lenz],dtype=float)
    for i in range(lenz):
        dx_var[:,i]=(var[1:,i]-var[:-1,i])/gr.dx
    return dx_var
        
def ddz(var):
    lenx=len(var[:,0])
    lenz=len(var[0,:])
    dz_var = np.empty([lenx,lenz-1],dtype=float)
    for i in range(lenx):
        dz_var[i,:]=(var[i,1:]-var[i,:-1])/gr.dz
    return dz_var        
        
