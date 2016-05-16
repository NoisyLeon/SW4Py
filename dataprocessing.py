import symdata
# import noisepy as npy
import obspy
import obspy.geodetics as obsGeo
dbase=symdata.sw4ASDF('sw4synthetics.h5')
inftan=symdata.InputFtanParam()
inftan.ffact=5.
inftan.pmf=True
# del dbase.auxiliary_data.DISPbasic1
# del dbase.auxiliary_data.DISPpmf2
# dbase.aftan(inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)
# del dbase.events
# del dbase.auxiliary_data.DISPbasic1interp
# dbase.InterpDisp(data_type='DISPpmf2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# dbase.AddEvent(x=1000, y=1000, z=0)

# del dbase.auxiliary_data.FieldDISPbasic1interp
# dbase.GetField(outdir='./', fieldtype='amp', data_type='DISPpmf2')
# dbase.GetField(outdir='./', fieldtype='Vgr',  data_type='DISPpmf2')

flst=dbase.aftanParallel(inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)