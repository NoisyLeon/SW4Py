import symdata
# import noisepy as npy
import obspy
# dbase1=symdata.sw4ASDF('sw4synthetics_homo_001.h5')
dbase=symdata.sw4ASDF('sw4synthetics_homo_002.h5')
# inv=dbase.Readsac('ak135_station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir/homo_001', comptype='u', verbose=True)
# dbase.AddEvent(x=1000, y=1000, z=0)
inftan=symdata.InputFtanParam()
inftan.ffact=1.
inftan.pmf=True

# del dbase.auxiliary_data.DISPbasic1
# del dbase.auxiliary_data.DISPbasic2
# del dbase.auxiliary_data.DISPpmf1
# del dbase.auxiliary_data.DISPpmf2

# dbase.aftan(inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)
# # del dbase.events
# del dbase.auxiliary_data.DISPpmf2interp
# dbase.InterpDisp(data_type='DISPpmf2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# dbase.AddEvent(x=1000, y=1000, z=0)

# del dbase.auxiliary_data.FieldDISPbasic1interp
# dbase.GetField(outdir='./homo', fieldtype='amp', data_type='DISPpmf2')
# dbase.GetField(outdir='./homo', fieldtype='Vgr',  data_type='DISPpmf2')
# 
dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/sw4_working_dir/DISP', inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)