import vmodel

rmodel=vmodel.rModel();
# rmodel.AddSingleBlock(ni=101, nj=101, nk=20, hh=20, hv=20, z0=None, #data=np.array([]),
#         vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
rmodel.AddTopoBlock(ni=101, nj=101,hh=20)
rmodel.AddSingleBlock(ni=101, nj=101, nk=20, hh=20, hv=20, z0=0, #data=np.array([]),
        vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
# rmodel.AddSingleBlock(ni=101, nj=101, nk=20, hh=20, hv=20, z0=0, #data=np.array([]),
#         vs=3.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
rmodel.AddSingleBlock(ni=101, nj=101, nk=40, hh=20, hv=20, z0=None, #data=np.array([]),
        vs=3.5, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
rmodel.AddSingleBlock(ni=101, nj=101, nk=11, hh=20, hv=20, z0=None, #data=np.array([]),
        vs=4.0, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None);
rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')

rmodel2=vmodel.rModel('/lustre/janus_scratch/life9360/sw4_working_dir/test_rmodel_001.rfile')
