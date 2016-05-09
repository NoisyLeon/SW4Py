import numpy as np


vs=np.ones(10)*3.5
vp=np.ones(10)*4.9
rho=np.ones(10)*2.7
vpr=np.append(vs,vp);
vpr=np.append(vpr,rho);
vpr=vpr.reshape(3,10)
vpr=vpr.T
vpr=vpr.reshape(30)
