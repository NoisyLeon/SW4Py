import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



infname = './ak135_VF_0km/Tph_10.0.txt'
inArr=np.loadtxt(infname)
plt.figure();
T=inArr[:,2]
DistArr=inArr[:,3]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
T=T[ind]
diffD=DistArr[1:] - DistArr[:-2]
diffT=T[1:] - T[:-2]
VphArr=diffD/diffT
DistArr=DistArrp[1:]

mindist=DistArr.min()
indexmin=DistArr.argmin()
plt.plot(DistArr, VphArr,'o' );
# plt.plot(DistArr, (VgrArr-VgrArr[indexmin])/VgrArr[indexmin]*100.,'o' );
# plt.ylabel('Relative Difference in Vgr (%)');
plt.ylabel('Vgr(km/s)');
plt.xlabel('Distance(km)');

plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err
    
    
# InstaStream.plot(type='section', norm_method='stream', alpha=1)
#     
