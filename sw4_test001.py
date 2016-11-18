import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



infname = '/lustre/janus_scratch/life9360/sw4_working_dir_Q/field_R_200km_zmax_4km_qs_ak135_vs_2000/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
plt.figure();
T=inArr[:,2]
DistArr=inArr[:,3]
VgrArr=DistArr/T
mindist=DistArr.min()
indexmin=DistArr.argmin()
plt.plot(DistArr, T,'o' );
# plt.plot(DistArr, (VgrArr-VgrArr[indexmin])/VgrArr[indexmin]*100.,'o' );
# plt.ylabel('Relative Difference in Vgr (%)');
plt.ylabel('Vgr(km/s)');
plt.ylabel('Travel time(sec)');
plt.xlabel('Distance(km)');

infname = '/lustre/janus_scratch/life9360/sw4_working_dir_Q/field_R_200km_zmax_4km_qs_ak135_vs_2000/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]
DistArr=inArr2[:,3]
fig, ax=plt.subplots()
mindist=DistArr.min()
indexmin=DistArr.argmin()
maxamp=AmpArr[indexmin]
# plt.plot(DistArr, AmpArr*1e9,'o' );
plt.ylabel('Amplitude');
plt.xlabel('Distance(km)')

CampArr=AmpArr*np.sqrt(DistArr/mindist )  /AmpArr[indexmin]
# CampArr=AmpArr*DistArr/ DistArr[0] 
# plt.plot(DistArr, (CampArr-CampArr[indexmin])/CampArr[indexmin]*100.,'o' );
plt.plot(DistArr, CampArr,'o' );
y1=900
y2=1300
ax.fill_betweenx(np.array([0.2, 1.1]), y1, y2, facecolor='red', alpha=0.5)
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err

