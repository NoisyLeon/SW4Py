import vmodel, stations


SLst=stations.StaLst();
# SLst.HomoStaLst(xmin=200, Nx=100, dx=200, ymin=200, Ny=100, dy=200);
# SLst.LineStaLst(xmin=1000000., Nx=200, dx=5000, y=1000000.)
SLst.LineStaLst(xmin=1000000., Nx=200, dx=5000, y=1000000.)
# SLst.ReadStaList('station.lst_001');
# SLst.WriteStaList('station.lst_003')
# SLst.Wrtie2Input(infname='sw4.input')
SLst.Wrtie2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in')
# SLst.Wrtie2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in')


BLst=vmodel.BlockLst();
# BLst.AddSingle(vs=3500)
BLst.ak135(zmax=200.);
BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_test.in');
# BLst.Write2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/homo_test.in');
