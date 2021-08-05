# coding=latin-1

# 2-D CODE FOR THE SIMULATION OF A COMPRESSIBLE FLUID

# ---------------------------------------------------------
# THIS IS THE MAIN PROGRAM OF THE CODE. 
#   -> It constitutes the decision and task-distribution center. 
#   -> It carries out a limited number of assignments itself.
# ---------------------------------------------------------

# --------------------------------------
# PROGRAM INITIALIZATION
# --------------------------------------

# CHECK IF AN INPUT FILE HAS BEEN SPECIFIED
from sys import argv        # argv contains the name of the program and the command line arguments
if len(argv)==1:            # No file specified...
    exit('Please, specify input file...')
else:
    filename=argv[1]        # The filename identifies the experiment being carried out throughout the program.
                            #   (See, for instance, exp_ctrl.py).

# READ INPUT FILE AND SET UP PROBLEM PARAMETERS
from read_input import *
read_input(filename)

# Import modules containing common variables, parameters and routines

import numpy as np              # Mathematical package to perform IDL-style array operations
from routines import *          # Module containing references to all the routines called from __main__
import grid_common as gr        # Module containing the numerical mesh and related parameters        
import params_common as par     # Module containing problem parameters, such as equilibrium density, etc...

# PRINT OUT INFORMATION ABOUT THE EXPERIMENT
print
print 'Experiment '+filename+' on course...'
print 'Taking '+str(gr.nintx)+'x'+str(gr.nintz)+' cells with '+str(gr.x0)+' < x < '+str(gr.xf)+' and '+str(gr.z0)+' < z < '+str(gr.zf)+'.'

# After this call to read_input(), the following variables are
# available in the __main__ module's symbol table:
#   itmax               -> Maximum number of iterations
#   timefinal           -> Time at which simulation ends
#   plottimeinterv      -> ]
#   plotstinterv        -> ] Cadence for output


# --------------------------------------
# PROGRAM INITIALIZATION - END ---------
# --------------------------------------


# -------------------------
#  INITIAL CONDITIONS 
# -------------------------

vx0,vz0,rho0,pres0=initcond()

momx0 = rho0*vx0                                            # Initial x-momentum density
momz0 = rho0*vz0                                            # Initial z-momentum density
energ0 = pres0/(par.gam-1.) + rho0*(vx0*vx0+vz0*vz0)/2.     # Initial energy density

# -------------------------
#  INITIAL CONDITIONS - end
# -------------------------
 
# --------------------------------------------------
# INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP
# --------------------------------------------------

vx=np.copy(vx0)
vz=np.copy(vz0)
pres=np.copy(pres0)

rho=np.copy(rho0)
momx=np.copy(momx0)
momz=np.copy(momz0)
energ=np.copy(energ0)

rhon=np.copy(rho0)
momxn=np.copy(momx0)
momzn=np.copy(momz0)
energn=np.copy(energ0)

it=0
time=0.
fstop=False

# --------------------------------------------------------
# INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP - end
# --------------------------------------------------------

# PLOT INITIAL CONDITIONS:
if plottimeinterv !=0:
    plotrout(vx0,vz0,rho0,pres0,time,it)
    raw_input('press enter...')




# ----------------------------
# BIG LOOP BEGINS
# ----------------------------

while it<itmax and time<timefinal and not fstop: # fstop === Force stop

#   CALCULATE FLUXES FROM DENSITIES
    flm,flcx,flcz,fle=fluxes(rho,momx,momz,energ)
    
#   TIMESTEP
    dt=cfl(rho,momx,momz,pres)
    
#   UPDATE TO HALF TIMESTEP  
    rhon,momxn,momzn,energn=update(dt,rho,momx,momz,energ,flm,flcx,flcz,fle,'half')

#   CALCULATE FLUXES FROM DENSITIES IN NEW GRID POINTS
    flm,flcx,flcz,fle=fluxes(rhon,momxn,momzn,energn)
    
#   UPDATE TO FULL TIMESTEP    
    rhon,momxn,momzn,energn=update(dt,rho,momx,momz,energ,flm,flcx,flcz,fle,'full') 

#   BOUNDARY CONDITIONS
    rhon=bcs(rhon,'periodic')
    momxn=bcs(momxn,'periodic')
    momzn=bcs(momzn,'periodic')    
    energn=bcs(energn,'periodic')

#   CALCULATE PRIMITIVE VARIABLES FROM THE DENSITIES
    vxn,vzn,presn=dens2vpt(rhon,momxn,momzn,energn)
    
#   EXCHANGE NEW AND OLD VARIABLES
    rho=rhon
    momx=momxn
    momz=momzn
    energ=energn
    vx=vxn
    vz=vzn
    pres=presn

#   UPDATE TIME AND NUMBER OF ITERATIONS
    it+=1
    time+=dt

#   STORE RESULTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#   CHECK IF EXPERIMENT IS COMPLETE
    fstop=exp_ctrl(filename) 

#   PLOT RESULTS 
    if plottimeinterv != 0 and (it%plottimeinterv==0 or it==itmax or time>timefinal or fstop):
        plotrout(vx,vz,rho,pres,time,it) 
         
# ---------------
# BIG LOOP - end
# ---------------

print('program finished.')
print('     it='+str(it))
print('   time='+str(time))
print
raw_input('press enter...')
print

