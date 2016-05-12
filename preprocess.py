import vmodel, stations


# SLst=stations.StaLst();
# SLst.HomoStaLst(xmin=50000, Nx=190, dx=10000, ymin=50000, Ny=190, dy=10000);
# # SLst.LineStaLst(xmin=1000000., Nx=200, dx=5000, y=1000000.)
# # SLst.ReadStaList('station.lst');
# SLst.WriteStaList('ak135_station.lst')
# # SLst.Wrtie2Input(infname='sw4.input')
# SLst.Wrtie2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in')
# SLst.Wrtie2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in')


BLst=vmodel.BlockLst();
# BLst.AddSingle(vs=3.5)
BLst.ak135(zmax=400.);
BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in');
# BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in');
