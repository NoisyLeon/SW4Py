import vmodel, stations


SLst=stations.StaLst()
SLst.HomoStaLst(xmin=20000, Nx=149, dx=20000, ymin=20000, Ny=29, dy=20000);
# SLst.WriteStaList('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/station_ak135_4mem2d.lst')
# SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/mem2d_test_001.in')
# 
# rmodel=vmodel.rModel()
# rmodel.ak135(zmin=0., zmax=200., ni=3001, nj=601, hh=1000., hv=1000.)
# rmodel.checkInput('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/mem2d_test_001.in')
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.CynlinderRingBasin(x0=2300000, y0=300000, zmax=10000., Rmax=100000., vs=2000., outfname='./cpsin.txt')
# 
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.writeVprofile('./cpsinput_20km_0.1.txt')
# 
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/model_20km_0.1.rfile')
# 
