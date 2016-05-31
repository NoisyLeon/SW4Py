
import symdata

# dbase=symdata.sw4ASDF('/home/lili/sw4synthetics_ak135_EX_z_1km.h5')
dbase=symdata.sw4ASDF('../sw4synthetics_ak135_0.1_homo_y_300.h5')
# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
dbase.PlotStreamsDistance()
