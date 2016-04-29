import stations

SLst=stations.StaLst();
SLst.HomoStaLst(xmin=200, Nx=100, dx=200, ymin=200, Ny=100, dy=200);
# SLst.LineStaLst(xmin=0, Nx=100, dx=1000, y=33000)
# SLst.ReadStaList('station.lst_001');
# SLst.WriteStaList('station.lst_003')
# SLst.Wrtie2Input(infname='sw4.input')
SLst.Wrtie2Input(infname='/lustre/janus_scratch/life9360/sw4_working_dir/LOH.1-h50.in')