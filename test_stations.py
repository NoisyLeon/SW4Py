import stations

SLst=stations.StaLst();
SLst.ReadStaList('station_ak135_4mem2d.lst')
nSLst=SLst.GetLineStaLst(y=300000)
nSLst.WriteStaList('station_ak135_4mem2d_y_300.lst')