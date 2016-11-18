import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



infname = '/lustre/janus_scratch/life9360/sw4_working_dir_Q/field_R_200km_zmax_4km_qs_40_vs_2000/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
plt.figure();
T=inArr[:,2]
DistArr=inArr[:,3]
VgrArr=DistArr/T
mindist=DistArr.min()
indexmin=DistArr.argmin()
plt.plot(DistArr, T,'o' );
plt.ylabel('Vgr(km/s)');
plt.ylabel('Travel time(sec)', fontsize=20);
plt.xlabel('Distance(km)', fontsize=20);

infname = '/lustre/janus_scratch/life9360/sw4_working_dir_Q/field_R_200km_zmax_4km_qs_40_vs_2000/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]
DistArr=inArr2[:,3]
ind=DistArr.argsort()
DistArr=DistArr[ind]
AmpArr=AmpArr[ind]
T=T[ind]

infname = '/lustre/janus_scratch/life9360/sw4_working_dir_Q/field_R_200km_zmax_4km_qs_ak135_vs_2000/Amp_10.0.txt'
inArr3=np.loadtxt(infname)
AmpArr2=inArr3[:,2]
DistArr2=inArr3[:,3]
ind2=DistArr2.argsort()
DistArr2=DistArr2[ind2]
AmpArr2=AmpArr2[ind2]
Q=200

fig, ax=plt.subplots()
mindist=DistArr.min()
indexmin=DistArr.argmin()
maxamp=AmpArr[indexmin]
AmpArr0=np.ones(AmpArr.size)/np.sqrt(DistArr/mindist ) *np.exp(-np.pi/10.*(T-T[indexmin])/Q)
plt.ylabel('Amplitude', fontsize=20);
plt.xlabel('Distance(km)', fontsize=20);


AmpArr=AmpArr /AmpArr[indexmin]
AmpArr2=AmpArr2 /AmpArr2[indexmin]
plt.plot(DistArr, AmpArr/AmpArr0,'bo' , ms=10)
plt.plot(DistArr, AmpArr2/AmpArr0,'go' , ms=10)
# plt.plot(DistArr, AmpArr0,'r--', lw=3 )
y1=800
y2=1200
ax.fill_betweenx(np.array([0., 1.1]), y1, y2, facecolor='red', alpha=0.5)
plt.ylim([0, 1.1])
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.xticks(np.arange(0., 3000., 200.))
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err

