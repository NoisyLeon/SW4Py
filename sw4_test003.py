import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from obspy.geodetics import kilometer2degrees


infname = './ak135_VF_0km/Tgr_10.0.txt'
# infname = './syndata_dir_000_DD/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
fig1, ax1 = plt.subplots()
T=inArr[:,2]
DistArr=inArr[:,3]
T=T[DistArr<1340]
DistArr=DistArr[DistArr<1340]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
T=T[ind]
DeltaArr=kilometer2degrees(DistArr)
VgrArr=DistArr/T
ax1.plot(DistArr, VgrArr,'--o' , markersize=10);
plt.ylabel('Vgr (km/s) ', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);
plt.title('Group Velocity', fontsize=30);
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)


fig2, ax2 = plt.subplots()
infname = './ak135_VF_0km/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]*1000.
DistArr=inArr2[:,3]
AmpArr=AmpArr[DistArr<1340]
DistArr=DistArr[DistArr<1340]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
AmpArr=AmpArr[ind]
mindist=DistArr.min()
indexmin=DistArr.argmin()
maxamp=AmpArr[indexmin]
print mindist, maxamp
plt.ylabel('Amplitude', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);

ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)

CampArr=AmpArr*np.sqrt( DistArr/mindist )  / maxamp
# plt.plot(DistArr,AmpArr,'g--o', markersize=10 );
plt.plot(DistArr,CampArr-0.1,'r--o', markersize=10 );
# plt.title('Amplitude ', fontsize=30);
plt.title('Corrected Normalized Amplitude', fontsize=30);
plt.show()

#     
