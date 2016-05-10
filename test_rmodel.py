import vmodel

rmodel=vmodel.rModel();
rmodel.ak135(zmax=200., ni=631, nj=601, hh=1000., hv=1000.)
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
# rmodel.CylinderCosineAnomaly(x0=500, y0=500, R=400, dm=-0.2, mtype=2, zmin=200, zmax=800, nb=None)
rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')

# rmodel2=vmodel.rModel('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')
