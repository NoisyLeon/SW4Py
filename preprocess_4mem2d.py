import vmodel, stations


# SLst=stations.StaLst()
# # SLst.HomoStaLst(xmin=20000, Nx=149, dx=20000, ymin=20000, Ny=29, dy=20000)
# SLst.LineStaLst(xmin=300000, Nx=133, dx=20000, y=300000)
# SLst.WriteStaList('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/station_4mem2d.lst')
# # SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/single_ring_basin.in')
# SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/single_homo_basin.in')
# # 
rmodel=vmodel.rModel()
rmodel.ak135(zmin=0., zmax=100., ni=3001, nj=601, hh=1000., hv=1000., CPS=True)
rmodel.checkInput('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/single_staircase_basin.in')
# rmodel.checkInput('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/single_homo_basin.in')
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderCosineAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='rho', zmin=0, zmax=20000, nb=2)

# rmodel.CynlinderRingBasin(x0=1100000, y0=300000, zmax=4000., Rmax=200000., vs=2000., outfname='./cpsin.txt')
# rmodel.CylinderHomoAnomaly(x0=1100000, y0=300000, zmax=4000., R=200000., dm=-0.1)
rmodel.CylinderLinearDepthAnomalyAll(x0=1100000, y0=300000, R=100000, vt=2000., vb=3000., zmax=5000, zmin=0, nb=None, outfname='cpsin_staircase.txt')
# rmodel.CylinderHomoSediment(x0=1100000, y0=300000, R=100000,  zmax=3000, vs=2000.)


# 
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=2300000, y0=300000, R=100000,  dm=0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vs', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='vp', zmin=0, zmax=20000, nb=2)
# # rmodel.CylinderHomoAnomaly(x0=700000, y0=300000, R=100000,  dm=-0.1, mname='rho', zmin=0, zmax=20000, nb=2)
# rmodel.writeVprofile('./cpsinput_4km_0.3.txt')
# 
rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/single_staircase_basin.rfile')
# 
