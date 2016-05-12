import symdata
# import noisepy as npy
import obspy
import obspy.geodetics as obsGeo
dbase=symdata.sw4ASDF('sw4synthetics.h5')

# inv=dbase.Readsac('station.lst',datadir='/lustre/janus_scratch/life9360/sw4_working_dir_trials/ak135', comptype='u')
# dbase.AddEvent(x=300, y=300, z=0)
# dbase.aftan()
f='/lustre/janus_scratch/life9360/sw4_working_dir_trials/ak135/SW4.165S1000.u'
prefname='./ak135.disp'
st=obspy.read(f)
tr=st[0]
tr1=symdata.sw4trace(tr.data, tr.stats)
dist, az, baz=obsGeo.base.gps2dist_azimuth(tr1.stats.sac.evla, tr1.stats.sac.evlo, tr1.stats.sac.stla, tr1.stats.sac.stlo ) # distance is in m
tr1.stats.sac.dist=1677.98241945
tr1.aftan(piover4=-1., pmf=True, phvelname=prefname)