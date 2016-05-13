import vmodel, stations


SLst=stations.StaLst();
SLst.HomoStaLst(xmin=50000, Nx=38, dx=50000, ymin=50000, Ny=38, dy=50000);
# # SLst.LineStaLst(xmin=1000000., Nx=200, dx=5000, y=1000000.)
# # SLst.ReadStaList('station.lst');
SLst.WriteStaList('ak135_station.lst')
# # SLst.Write2Input(infname='sw4.input')
SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in')
# SLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in')


BLst=vmodel.BlockLst();
# BLst.AddSingle(vs=3.5)
BLst.ak135(zmax=200.);
BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in');
# BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in');
