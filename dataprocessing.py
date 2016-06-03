import symdata
# import noisepy as npy
import obspy

dbase=symdata.sw4ASDF('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/sw4synthetics_ak135_0.1_homo_y_300.h5')
# # # dbase=symdata.sw4ASDF('./sw4synthetics.h5')
dbase.Readsac('station_ak135_4mem2d_y_300.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/ak135_20km_dm_0.1',
            comptype='u', verbose=True)
inv=dbase.Readsac('ak135_station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir/homo_001', comptype='u', verbose=True)
try:
    del dbase.events
except:
    pass
dbase.AddEvent(x=1500, y=300, z=1)
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

# dbase.aftan(inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)
dbase.aftanMP(tb=-13.5, outdir='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/DISP', inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)

try:
    del dbase.auxiliary_data.DISPpmf2interp
except:
    pass
dbase.InterpDisp(data_type='DISPpmf2')
# ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')

# inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# dbase.AddEvent(x=1000, y=1000, z=0)

del dbase.auxiliary_data.FieldDISPbasic1interp
dbase.GetField(outdir='./ak135_4mem2d_20km_0.1', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='./ak135_4mem2d_20km_0.1', fieldtype='Vgr',  data_type='DISPpmf2')
# dbase.GetField(outdir='./ak135_VF', fieldtype='Vgr', data_type='DISPpmf2')
# dbase.GetField(outdir='./ak135_VF', fieldtype='Vph',  data_type='DISPpmf2')
# # 
