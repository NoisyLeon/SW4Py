import symdata
# import noisepy as npy
import obspy

dbase=symdata.sw4ASDF('../sw4synthetics_ak135_VF_0km.h5')
# dbase.Readsac('ak135_station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir/ak135_VF_0km',
#             comptype='u', verbose=True)

# try:
#     del dbase.events
# except:
#     pass
# dbase.AddEvent(x=1000, y=1000, z=0)
inftan=symdata.InputFtanParam()
inftan.ffact=1.
inftan.pmf=True
try:
    del dbase.auxiliary_data.DISPbasic1
    del dbase.auxiliary_data.DISPbasic2
    del dbase.auxiliary_data.DISPpmf1
    del dbase.auxiliary_data.DISPpmf2
except:
    pass


dbase.aftanMP(tb=-8.48550289567, outdir='/lustre/janus_scratch/life9360/sw4_working_dir/DISP', inftan=inftan, basic2=True,
            pmf1=True, pmf2=True)

try:
    del dbase.auxiliary_data.DISPpmf2interp
except:
    pass
dbase.InterpDisp(data_type='DISPpmf2')

# dbase.GetField(outdir='./ak135_VF_0km', fieldtype='amp', data_type='DISPpmf2')
# dbase.GetField(outdir='./ak135_VF_0km', fieldtype='Vgr',  data_type='DISPpmf2')
dbase.GetField(outdir='./ak135_VF_0km', fieldtype='Vph',  data_type='DISPpmf2')
# dbase.GetField(outdir='./ak135_VF', fieldtype='Vgr', data_type='DISPpmf2')
# dbase.GetField(outdir='./ak135_VF', fieldtype='Vph',  data_type='DISPpmf2')
# # 
