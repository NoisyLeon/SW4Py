import symdata
# import noisepy as npy
import obspy
# dbase1=symdata.sw4ASDF('sw4synthetics_homo_001.h5')
dbase=symdata.sw4ASDF('/lustre/janus_scratch/life9360/sw4_working_dir/sw4synthetics_ak135_EX_z_1km.h5')

# dbase.Readsac('ak135_station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_EX_z_1km', comptype='u', verbose=True)
# inv=dbase.Readsac('ak135_station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir/homo_001', comptype='u', verbose=True)
try:
    del dbase.events
except:
    pass
dbase.AddEvent(x=1000, y=1000, z=0)
inftan=symdata.InputFtanParam()
inftan.ffact=5.
inftan.pmf=True
try:
    del dbase.auxiliary_data.DISPbasic1
    del dbase.auxiliary_data.DISPbasic2
    del dbase.auxiliary_data.DISPpmf1
    del dbase.auxiliary_data.DISPpmf2
except:
    pass

# # dbase.aftan(inftan=inftan, basic2=True,
# #             pmf1=True, pmf2=True)
dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/sw4_working_dir/DISP', inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)

try:
    del dbase.auxiliary_data.DISPpmf2interp
except:
    pass
dbase.InterpDisp(data_type='DISPpmf2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# dbase.AddEvent(x=1000, y=1000, z=0)

# del dbase.auxiliary_data.FieldDISPbasic1interp
dbase.GetField(outdir='./ak135_EX', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='./ak135_EX', fieldtype='Vgr',  data_type='DISPpmf2')
# 
