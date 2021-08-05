# ----------------------------------------
#
# MODULE: grid_common
#   Contains grid variables and parameters
#
# ----------------------------------------

import numpy as np

nintx=100
x0=-5.
xf=5.
npx=nintx+2
xmid=(x0+xf)/2.
xlen=(xf-x0)
dx=xlen/nintx
xx=np.linspace(x0-dx/2.,xf+dx/2.,npx)

nintz=100
z0=-5.
zf=5.
npz=nintz+2
zmid=(z0+zf)/2.
zlen=(zf-z0)
dz=zlen/nintz
zz=np.linspace(z0-dz/2.,zf+dz/2.,npz)

