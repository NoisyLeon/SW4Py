import vmodel, stations


SLst=stations.StaLst()
# # SLst.HomoStaLst(xmin=300000, Nx=230, dx=20000, ymin=100000, Ny=29, dy=20000)
SLst.LineStaLst(xmin=300000, Nx=133, dx=20000, y=300000)
# # SLst.HomoStaLst(xmin=300000, Nx=230, dx=20000, ymin=100000, Ny=29, dy=20000)
SLst.WriteStaList('/lustre/janus_scratch/life9360/sw4_working_dir_Q/sta_Q.lst')
# SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir_Q/Q_run_001.in')
# # 
rmodel=vmodel.rModel(attenuation = 1 )
rmodel.ak135(zmin=0., zmax=200., ni=3001, nj=601, hh=1000., hv=1000.)
rmodel.checkInput('/lustre/janus_scratch/life9360/sw4_working_dir_Q/Q_run_001.in')

# rmodel.CylinderCosineSediment(x0=1100000, y0=300000, R=200000, zmax=4000, vs=2000., qs=40.)
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_Q/R_200km_zmax_4km_qs_40_vs_2000.rfile')
# rmodel.CylinderCosineSediment(x0=1100000, y0=300000, R=200000, zmax=4000, vs=2000.)
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_Q/R_200km_zmax_4km_qs_ak135_vs_2000.rfile')
# rmodel.CylinderCosineSediment(x0=1100000, y0=300000, R=200000, zmax=4000,  qs=40.)
# rmodel.write('/lustre/janus_scratch/life9360/sw4_working_dir_Q/R_200km_zmax_4km_qs_40_vs_ak135.rfile')
# # 
# # # rmodel.writeVprofile('./cpsinput_20km_0.1.txt')

