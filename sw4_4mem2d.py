import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



infname = './ak135_4mem2d_20km_0.1/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
plt.figure();
T=inArr[:,2]
XArr=inArr[:,0]
index100=np.where(XArr==1300.)
DistArr=inArr[:,3]
VgrArr=DistArr/T

# plt.plot(XArr, VgrArr,'o' );
plt.plot(XArr, T,'o' );
# plt.plot(XArr, (VgrArr-VgrArr[indexmin])/VgrArr[indexmin]*100.,'o' );
# plt.ylabel('Relative Difference in Vgr (%)');
# plt.ylabel('Vgr(km/s)');
# plt.xlabel('Distance(km)');



infname = './ak135_4mem2d_20km_0.1/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]
XArr=inArr[:,0]
plt.figure();
amp100=AmpArr[index100]
# mindist=XArr.min()
# indexmin=XArr.argmin()
# maxamp=AmpArr[indexmin]
# plt.plot(XArr, AmpArr*1e9,'o' );
plt.ylabel('Amplitude');
plt.xlabel('Distance(km)');
# 
# plt.figure();
# # plt.plot(XArr, VgrArr, 'x');
# # plt.plot(XArr, AmpArr*np.sqrt(np.sin(DeltaArr*np.pi/180.) ));
# # CampArr=AmpArr*np.sqrt(np.sin(DeltaArr*np.pi/180.) )/np.sqrt(np.sin(DeltaArr[0]*np.pi/180.) )
# CampArr=AmpArr*np.sqrt(XArr/mindist )  
# CampArr=AmpArr*XArr/ XArr[0] 
# plt.plot(XArr, (CampArr-CampArr[indexmin])/CampArr[indexmin]*100.,'o' );
plt.plot(XArr, AmpArr/amp100,'o' );
# # plt.
# plt.ylabel('Relative Difference in Corrected Amp (%)');
# plt.xlabel('Distance(km)');
# # plt.axis([ XArr.min(), XArr.max(), CampArr.min(), CampArr.max()])
plt.show()
DistArr=inArr[:,3]
# slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
# print slope, intercept, r_value, p_value, std_err
    
    
# InstaStream.plot(type='section', norm_method='stream', alpha=1)
#     
