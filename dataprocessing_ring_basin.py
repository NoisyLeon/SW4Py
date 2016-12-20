import symdata
# import noisepy as npy
import obspy

dtype='single_staircase_basin'
dbase=symdata.sw4ASDF('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/sw4synthetics_'+dtype+'.h5')

dbase.Readsac('/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/station_4mem2d.lst',\
            datadir='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/'+dtype, comptype='u', verbose=True)
try:
    del dbase.events
except:
    pass
dbase.AddEvent(x=100, y=300, z=1)
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
# 
dbase.aftan(inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)
# dbase.aftanMP(tb=-13.5, outdir='/lustre/janus_scratch/life9360/sw4_working_dir_Q/DISP', inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)
# 
try:
    del dbase.auxiliary_data.DISPpmf2interp
except:
    pass
dbase.InterpDisp(data_type='DISPpmf2')
# # ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')
# 
# # inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# # dbase.AddEvent(x=1000, y=1000, z=0)
# 
# del dbase.auxiliary_data.FieldDISPbasic1interp
dbase.GetField(outdir='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/field_'+dtype,\
               fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/field_'+dtype,\
               fieldtype='Vgr',  data_type='DISPpmf2')

# # 
