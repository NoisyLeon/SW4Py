import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



Dx=600.
infname = '/lustre/janus_scratch/life9360/specfem2d_working_dir/2d_3d_staircase_basin/field_data/Amp_10.0.txt'
inArr=np.loadtxt(infname)

XArr=inArr[:,0]
distArr=inArr[:,3]
AmpArr=inArr[:,2]

ind=np.argsort(XArr)
XArr=XArr[ind]
AmpArr=AmpArr[ind]
distArr=distArr[ind]

XArr2= XArr
AmpArr2= np.ones(XArr.size) / np.sqrt(distArr)

index=np.where(XArr==1000)
amp=AmpArr[index]
index2=np.where(XArr2==1000.)
amp2=AmpArr2[index2]

Dx2=100.
infname = '/lustre/janus_scratch/life9360/sw4_working_dir_4mem2d/field_single_staircase_basin/Amp_10.0.txt'
inArr4=np.loadtxt(infname)
XArr3=inArr4[:,0]
AmpArr3=inArr4[:,2]
DistArr3=inArr4[:,3]
ind3=XArr3.argsort()
DistArr3=DistArr3[ind3]
AmpArr3=AmpArr3[ind3]
XArr3=XArr3[ind3]
index3=np.where(XArr3==500.)
amp3=AmpArr3[index3]


fig, ax=plt.subplots()
ax.plot(XArr-Dx, AmpArr/amp,'go', lw=3, markersize=10, label='2D' );
ax.plot(XArr2-Dx, AmpArr2/amp2,'r-', lw=5, markersize=10, label='1D ' );
ax.plot(XArr3-Dx2, AmpArr3/amp3,'bo--', lw=5, markersize=10, label='3D' );
# ax.plot(XArr2-Dx, AmpArr/amp/AmpArr2*amp2,'bo', lw=5, markersize=10);
# ax.plot(XArr, (TArr-TArr2-1.5),'ro', lw=3, markersize=10, label='difference')
# ax.plot(Dx, 1, 'y*', markersize=20)
plt.legend(loc='upper right', fontsize=25, numpoints = 1)
y1=900.
y2=1100.
ax.fill_betweenx(np.array([0, (AmpArr/amp).max()]), y1, y2, facecolor='red', alpha=0.5)
# y1=2200
# y2=2400
# ax.fill_betweenx(np.array([0.2, 1.1]), y1, y2, facecolor='blue', alpha=0.5)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.ylim([0, 1.05])
# plt.xlim([300., 3500.])
plt.xticks(np.arange(300., 3700., 200.))
# plt.ylim([(TArr-TArr2-1.5).min(), (TArr-TArr2).max()])
plt.ylabel('Normalized amplitude', fontsize=30);
# plt.xlabel('X position (km)', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);
# plt.title('Amplitude measurement with rectangle anomaly', fontsize=30);
# plt.title('Amplitude measurement with circular anomaly', fontsize=30);
plt.show()
