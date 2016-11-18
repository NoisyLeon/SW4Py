import vmodel
import numpy as np
rmodel=vmodel.rModel(attenuation =1)
rmodel.ak135(zmin=0., zmax=660., ni=631, nj=601, hh=1000., hv=1000.)
rmodel.CylinderCosineSediment(x0=100000, y0=200000, R=200000, zmax=4000, vs=2000., qs=40.)

# rmodel.CylinderCosineAnomaly(x0=100000, y0=200000, R=50000,  dm=-0.2, mname='vs', zmin=0, zmax=25000, nb=2)
# hArr=rmodel.CylinderCosineSediment(x0=100000, y0=200000, R=50000, vs=1500., zmax=10000.)
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_trials/test_rmodel_001.rfile')
# rmodel.AddSingleBlock(ni=101, nj=101, nk=20, hh=20, hv=20, z0=None, #data=np.array([]),
#         vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
# rmodel.AddTopoBlock(ni=131, nj=101,hh=20)
# rmodel.AddSingleBlock(ni=131, nj=101, nk=20, hh=20, hv=20, z0=0, #data=np.array([]),
#         vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
# # rmodel.AddSingleBlock(ni=101, nj=101, nk=20, hh=20, hv=20, z0=0, #data=np.array([]),
# #         vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
# rmodel.AddSingleBlock(ni=131, nj=101, nk=40, hh=20, hv=20, z0=None, #data=np.array([]),
#         vs=3.5, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
# rmodel.AddSingleBlock(ni=131, nj=101, nk=11, hh=20, hv=20, z0=None, #data=np.array([]),
#         vs=4.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);

# rmodel.BlockAnomaly(xmin=600, xmax=1500, ymin=300, ymax=1200, zmin=300, zmax=800, dm=-0.3)
# rmodel.CylinderCosineAnomaly(x0=200000, y0=300000, R=100000, dm=-0.2, mname='vs', zmin=30000, zmax=120000, nb=None)
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')

# rmodel2=vmodel.rModel('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')


# for i in np.arange(rmodel.nb-1)+1:
#     rblock=rmodel.rblocks[i]
#     print rblock.z0/1000., (rblock.z0+(rblock.nk-1)*rblock.hv)/1000., (rblock.nk-1)*rblock.hv/1000.
#     print (rblock.data[:,:,:,1]).max(), (rblock.data[:,:,:,2]).max(), (rblock.data[:,:,:,0]).max()
#     print '\n'