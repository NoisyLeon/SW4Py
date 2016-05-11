# import noisepy 
from pyasdf import ASDFDataSet
# import noisepy 
# import obspy
# f='/projects/life9360/code/CPSPy/sac_dir/B00511ZVF.sac'
# st=obspy.read(f)
# tr=st[0]
# tr1=noisepy.noisetrace(tr.data, tr.stats)
# tr1.aftan(piover4=-1., pmf=True, tmin=2.0, tmax=50.0)

import obspy.geodetics as obsGeo
import obspy.taup.taup
import obspy
import pyaftan as ftan  # Comment this line if you do not have pyaftan
import numpy as np
import glob, os
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.pylab as plb
import copy
import scipy.signal
import numexpr as npr
from functools import partial
import multiprocessing as mp
import math

# try:
#     import fftw3 # pyfftw3-0.2
#     useFFTW=True;
# except:
#     useFFTW=False;
useFFTW=False;
import time
import shutil
# import CURefPy as ref # Comment this line if you do not have CURefPy
from subprocess import call
from mpl_toolkits.basemap import Basemap
import warnings


#LDPATH = os.environ['LD_LIBRARY_PATH']
#sys.path.append(LDPATH)

class ftanParam(object):
    """
    Basic FTAN parameters:
    nfout1_1 - output number of frequencies for arr1, (integer*4)
    arr1_1   - preliminary results.
              Description: real*8 arr1(8,n), n >= nfin)
              arr1(1,:) -  central periods, s
              arr1(2,:) -  observed periods, s
              arr1(3,:) -  group velocities, km/s
              arr1(4,:) -  phase velocities, km/s or phase if nphpr=0, rad
              arr1(5,:) -  amplitudes, Db
              arr1(6,:) -  discrimination function
              arr1(7,:) -  signal/noise ratio, Db
              arr1(8,:) -  maximum half width, s
              arr1(9,:) -  amplitudes, nm/m
    arr2_1   - final results
    nfout2_1 - output number of frequencies for arr2, (integer*4)
              Description: real*8 arr2(7,n), n >= nfin)
              If nfout2 == 0, no final result.
              arr2(1,:) -  central periods, s
              arr2(2,:) -  observed periods, s
              arr2(3,:) -  group velocities, km/sor phase if nphpr=0, rad
              arr2(4,:) -  phase velocities, km/s
              arr2(5,:) -  amplitudes, Db
              arr2(6,:) -  signal/noise ratio, Db
              arr2(7,:) -  maximum half width, s
              arr2(8,:) -  amplitudes, nm/m
    tamp_1      -  time to the beginning of ampo table, s (real*8)
    nrow_1      -  number of rows in array ampo, (integer*4)
    ncol_1      -  number of columns in array ampo, (integer*4)
    amp_1       -  Ftan amplitude array, Db, (real*8)
    ierr_1   - completion status, =0 - O.K.,           (integer*4)
                                 =1 - some problems occures
                                 =2 - no final results
    ==========================================================
    Phase-Matched-Filtered FTAN parameters:
    nfout1_2 - output number of frequencies for arr1, (integer*4)
    arr1_2   - preliminary results.
             Description: real*8 arr1(8,n), n >= nfin)
             arr1(1,:) -  central periods, s (real*8)
             arr1(2,:) -  apparent periods, s (real*8)
             arr1(3,:) -  group velocities, km/s (real*8)
             arr1(4,:) -  phase velocities, km/s (real*8)
             arr1(5,:) -  amplitudes, Db (real*8)
             arr1(6,:) -  discrimination function, (real*8)
             arr1(7,:) -  signal/noise ratio, Db (real*8)
             arr1(8,:) -  maximum half width, s (real*8)
             arr1(9,:) -  amplitudes, nm/m
    arr2_2   - final results
    nfout2_2 - output number of frequencies for arr2, (integer*4)
             Description: real*8 arr2(7,n), n >= nfin)
             If nfout2 == 0, no final results.
             arr2(1,:) -  central periods, s (real*8)
             arr2(2,:) -  apparent periods, s (real*8)
             arr2(3,:) -  group velocities, km/s (real*8)
             arr1(4,:) -  phase velocities, km/s (real*8)
             arr2(5,:) -  amplitudes, Db (real*8)
             arr2(6,:) -  signal/noise ratio, Db (real*8)
             arr2(7,:) -  maximum half width, s (real*8)
             arr2(8,:) -  amplitudes, nm/m
    tamp_2      -  time to the beginning of ampo table, s (real*8)
    nrow_2      -  number of rows in array ampo, (integer*4)
    ncol_2      -  number of columns in array ampo, (integer*4)
    amp_2       -  Ftan amplitude array, Db, (real*8)
    ierr_2   - completion status, =0 - O.K.,           (integer*4)
                                =1 - some problems occures
                                =2 - no final results
    """
    def __init__(self):
        # Parameters for first iteration
        self.nfout1_1=0
        self.arr1_1=np.array([])
        self.nfout2_1=0
        self.arr2_1=np.array([])
        self.tamp_1=0.
        self.nrow_1=0
        self.ncol_1=0
        self.ampo_1=np.array([],dtype='float32')
        self.ierr_1=0
        # Parameters for second iteration
        self.nfout1_2=0
        self.arr1_2=np.array([])
        self.nfout2_2=0
        self.arr2_2=np.array([])
        self.tamp_2=0.
        self.nrow_2=0
        self.ncol_2=0
        self.ampo_2=np.array([])
        self.ierr_2=0
        # Flag for existence of predicted phase dispersion curve
        self.preflag=False

    def writeDISP(self, fnamePR):
        """
        Write FTAN parameters to DISP files given a prefix.
        fnamePR: file name prefix
        _1_DISP.0: arr1_1
        _1_DISP.1: arr2_1
        _2_DISP.0: arr1_2
        _2_DISP.1: arr2_2
        """
        if self.nfout1_1!=0:
            f10=fnamePR+'_1_DISP.0';
            
            # Lf10=self.nfout1_1;
            # outArrf10=np.arange(Lf10);
            # for i in np.arange(7):
            #     outArrf10=np.append(outArrf10, self.arr1_1[i,:Lf10]);
            # outArrf10=outArrf10.reshape((8,Lf10));
            # outArrf10=outArrf10.T;
            # np.savetxt(f10+'new', outArrf10, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
            
            f=open(f10,'w')
            for i in np.arange(self.nfout1_1):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr1_1[0,i] , self.arr1_1[1,i] , self.arr1_1[2,i] , self.arr1_1[3,i]  \
                    , self.arr1_1[4,i] , self.arr1_1[5,i] , self.arr1_1[6,i] )
                f.writelines(tempstr)
            f.close()
        if self.nfout2_1!=0:
            f11=fnamePR+'_1_DISP.1'
            
            # Lf11=self.nfout2_1;
            # outArrf11=np.arange(Lf11);
            # for i in np.arange(6):
            #     outArrf11=np.append(outArrf11, self.arr2_1[i,:Lf11]);
            # outArrf11=outArrf11.reshape((7,Lf11));
            # outArrf11=outArrf11.T;
            # np.savetxt(f11+'new', outArrf11, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
            
            f=open(f11,'w')
            for i in np.arange(self.nfout2_1):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr2_1[0,i], self.arr2_1[1,i] , self.arr2_1[2,i] , self.arr2_1[3,i]  \
                    , self.arr2_1[4,i] , self.arr2_1[5,i]  )
                f.writelines(tempstr)
            f.close()
        if self.nfout1_2!=0:
            f20=fnamePR+'_2_DISP.0';
            
            # Lf20=self.nfout1_2;
            # outArrf20=np.arange(Lf20);
            # for i in np.arange(7):
            #     outArrf20=np.append(outArrf20, self.arr1_2[i,:Lf20]);
            # outArrf20=outArrf20.reshape((8,Lf20));
            # outArrf20=outArrf20.T;
            # np.savetxt(f20+'new', outArrf20, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
            
            f=open(f20,'w')
            for i in np.arange(self.nfout1_2):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf \n' %( i, self.arr1_2[0,i], self.arr1_2[1,i] , self.arr1_2[2,i] , self.arr1_2[3,i]  \
                    , self.arr1_2[4,i] , self.arr1_2[5,i] , self.arr1_2[6,i] )
                f.writelines(tempstr)
            f.close()
        if self.nfout2_2!=0:
            f21=fnamePR+'_2_DISP.1';
            
            # Lf21=self.nfout2_2;
            # outArrf21=np.arange(Lf21);
            # for i in np.arange(6):
            #     outArrf21=np.append(outArrf21, self.arr2_2[i,:Lf21]);
            # outArrf21=outArrf21.reshape((7,Lf21));
            # outArrf21=outArrf21.T;
            # np.savetxt(f21+'new', outArrf21, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
            
            f=open(f21,'w')
            for i in np.arange(self.nfout2_2):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr2_2[0,i], self.arr2_2[1,i] , self.arr2_2[2,i] , self.arr2_2[3,i]  \
                    , self.arr2_2[4,i] , self.arr2_2[5,i]  )
                f.writelines(tempstr)
            f.close()
        return

    def FTANcomp(self, inftanparam, compflag=4):
        """
        Compare aftan results for two ftanParam objects.
        """
        fparam1=self
        fparam2=inftanparam
        if compflag==1:
            obper1=fparam1.arr1_1[1,:fparam1.nfout1_1]
            gvel1=fparam1.arr1_1[2,:fparam1.nfout1_1]
            phvel1=fparam1.arr1_1[3,:fparam1.nfout1_1]
            obper2=fparam2.arr1_1[1,:fparam2.nfout1_1]
            gvel2=fparam2.arr1_1[2,:fparam2.nfout1_1]
            phvel2=fparam2.arr1_1[3,:fparam2.nfout1_1]
        elif compflag==2:
            obper1=fparam1.arr2_1[1,:fparam1.nfout2_1]
            gvel1=fparam1.arr2_1[2,:fparam1.nfout2_1]
            phvel1=fparam1.arr2_1[3,:fparam1.nfout2_1]
            obper2=fparam2.arr2_1[1,:fparam2.nfout2_1]
            gvel2=fparam2.arr2_1[2,:fparam2.nfout2_1]
            phvel2=fparam2.arr2_1[3,:fparam2.nfout2_1]
        elif compflag==3:
            obper1=fparam1.arr1_2[1,:fparam1.nfout1_2]
            gvel1=fparam1.arr1_2[2,:fparam1.nfout1_2]
            phvel1=fparam1.arr1_2[3,:fparam1.nfout1_2]
            obper2=fparam2.arr1_2[1,:fparam2.nfout1_2]
            gvel2=fparam2.arr1_2[2,:fparam2.nfout1_2]
            phvel2=fparam2.arr1_2[3,:fparam2.nfout1_2]
        else:
            obper1=fparam1.arr2_2[1,:fparam1.nfout2_2]
            gvel1=fparam1.arr2_2[2,:fparam1.nfout2_2]
            phvel1=fparam1.arr2_2[3,:fparam1.nfout2_2]
            obper2=fparam2.arr2_2[1,:fparam2.nfout2_2]
            gvel2=fparam2.arr2_2[2,:fparam2.nfout2_2]
            phvel2=fparam2.arr2_2[3,:fparam2.nfout2_2]
        plb.figure()
        ax = plt.subplot()
        ax.plot(obper1, gvel1, '--k', lw=3) #
        ax.plot(obper2, gvel2, '-.b', lw=3)
        plt.xlabel('Period(s)')
        plt.ylabel('Velocity(km/s)')
        plt.title('Group Velocity Comparison')
        if (fparam1.preflag==True and fparam2.preflag==True):
            plb.figure()
            ax = plt.subplot()
            ax.plot(obper1, phvel1, '--k', lw=3) #
            ax.plot(obper2, phvel2, '-.b', lw=3)
            plt.xlabel('Period(s)')
            plt.ylabel('Velocity(km/s)')
            plt.title('Phase Velocity Comparison')
        return

class snrParam(object):
    """
    SNR parameters:
        suffix: p=positve lag;  n=negative lag; s=symmetric lag
        amp: largest amplitude measurement for each period
        snr: SNR measurement for each period
        nrms: noise rms measurement for each period
        oper: observed period

    """
    def __init__(self):
        self.amp_p=np.array([])
        self.snr_p=np.array([])
        self.nrms_p=np.array([])
        self.oper_p=np.array([])
        self.amp_n=np.array([])
        self.snr_n=np.array([])
        self.nrms_n=np.array([])
        self.oper_n=np.array([])
        self.amp_s=np.array([])
        self.snr_s=np.array([])
        self.nrms_s=np.array([])
        self.oper_s=np.array([])

    def writeAMPSNR(self, fnamePR):
        """
        writeAMPSNR:
        Write output SNR parameters to text files
        _pos_amp_snr - positive lag
        _neg_amp_snr - negative lag
        _amp_snr     - symmetric lag
        """
        len_p=len(self.amp_p)
        len_n=len(self.amp_n)
        len_s=len(self.amp_s)
        if len_p!=0:
            fpos=fnamePR+'_pos_amp_snr'
            f=open(fpos,'w')
            for i in np.arange(len_p):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(  self.oper_p[i] , self.amp_p[i],  self.snr_p[i] )
                f.writelines(tempstr)
            f.close()
        if len_n!=0:
            fneg=fnamePR+'_neg_amp_snr'
            f=open(fneg,'w')
            for i in np.arange(len_n):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(   self.oper_n[i] , self.amp_n[i],  self.snr_n[i] )
                f.writelines(tempstr)
            f.close()
        if len_s!=0:
            fsym=fnamePR+'_amp_snr'
            f=open(fsym,'w')
            for i in np.arange(len_s):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(   self.oper_s[i] , self.amp_s[i],  self.snr_s[i] )
                f.writelines(tempstr)
            f.close()
        return

class noisetrace(obspy.core.trace.Trace):
    """
    noisetrace:
    A derived class inherited from obspy.core.trace.Trace. This derived class have a variety of new member functions
    """
    def init_ftanParam(self):
        """
        Initialize ftan parameters
        """
        self.ftanparam=ftanParam()

    def init_snrParam(self):
        """
        Initialize SNR parameters
        """
        self.SNRParam=snrParam()

    def reverse(self):
        """
        Reverse the trace
        """
        self.data=self.data[::-1]
        return

    def makesym(self):
        """
        Turn the double lagged cross-correlation data to one single lag
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!');
        if self.stats.npts%2!=1:
            raise ValueError('Error: Incompatible begin and end time!');
        nhalf=(self.stats.npts-1)/2+1;
        neg=self.data[:nhalf];
        pos=self.data[nhalf-1:self.stats.npts];
        neg=neg[::-1];
        self.data=npr.evaluate( '(pos+neg)/2' );
        self.stats.npts=nhalf;
        self.stats.starttime=self.stats.starttime+self.stats.sac.e;
        self.stats.sac.b=0.;
        return

    def getneg(self):
        """
        Get the negative lag of a cross-correlation record
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!');
        negTr=self.copy();
        t=self.stats.starttime;
        L=(int)((self.stats.npts-1)/2)+1;
        negTr.data=negTr.data[:L];
        negTr.data=negTr.data[::-1];
        negTr.stats.npts=L;
        negTr.stats.sac.b=0.;
        negTr.stats.starttime=t-self.stats.sac.b;
        return negTr;

    def getpos(self):
        """
        Get the positive lag of a cross-correlation record
        """
        if abs(self.stats.sac.b+self.stats.sac.e)>self.stats.delta:
            raise ValueError('Error: Not symmetric trace!');
        posTr=self.copy();
        t=self.stats.starttime;
        L=(int)((self.stats.npts-1)/2)+1;
        posTr.data=posTr.data[L-1:];
        posTr.stats.npts=L;
        posTr.stats.sac.b=0.;
        posTr.stats.starttime=t-self.stats.sac.b;
        return posTr;

    def aftan(self, pmf=True, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, \
        tmax=30.0, tresh=20.0, ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0,phvelname=''):

        """ (Automatic Frequency-Time ANalysis) aftan analysis:
        This function read SAC file, make it symmtric (if it is not), and then do aftan analysis.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        pmf         - flag for Phase-Matched-Filtered output (default: True)
        piover4     - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin        - minimal group velocity, km/s
        vmax        - maximal group velocity, km/s
        tmin        - minimal period, s
        tmax        - maximal period, s
        tresh       - treshold for jump detection, usualy = 10, need modifications
        ffact       - factor to automatic filter parameter, usualy =1
        taperl      - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        snr         - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
        fmatch      - factor to length of phase matching window
        fname       - SAC file name
        phvelname   - predicted phase velocity file name
        
        Output:
        self.ftanparam, a object of ftanParam class, to store output aftan results
        -----------------------------------------------------------------------------------------------------
        References:
        Levshin, A. L., and M. H. Ritzwoller. Automated detection, extraction, and measurement of regional surface waves.
             Monitoring the Comprehensive Nuclear-Test-Ban Treaty: Surface Waves. Birkh?user Basel, 2001. 1531-1545.
        Bensen, G. D., et al. Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements.
             Geophysical Journal International 169.3 (2007): 1239-1260.
        """
        try:
            self.ftanparam
        except:
            self.init_ftanParam()
        try:
            dist=self.stats.sac.dist;
        except:
            dist=self.stats.distance;
        if (phvelname==''):
            phvelname='./ak135.disp';
        nprpv = 0
        phprper=np.zeros(300)
        phprvel=np.zeros(300)
        if os.path.isfile(phvelname):
            php=np.loadtxt(phvelname)
            phprper=php[:,0]
            phprvel=php[:,1]
            nprpv = php[:,0].size
            phprper=np.append(phprper,np.zeros(300-phprper.size))
            phprvel=np.append(phprvel,np.zeros(300-phprvel.size))
            self.ftanparam.preflag=True
        nfin = 64
        npoints = 5  #  only 3 points in jump
        perc    = 50.0 # 50 % for output segment
        tempsac=self.copy()
        if abs(tempsac.stats.sac.b+tempsac.stats.sac.e)<tempsac.stats.delta:
            tempsac.makesym()
        tb=tempsac.stats.sac.b
        length=len(tempsac.data)
        if length>32768:
            print "Warning: length of seismogram is larger than 32768!"
            nsam=32768
            tempsac.data=tempsac.data[:nsam]
            tempsac.stats.e=(nsam-1)*tempsac.stats.delta+tb
            sig=tempsac.data
        else:
            sig=np.append(tempsac.data,np.zeros( float(32768-tempsac.data.size) ) )
            nsam=int( float (tempsac.stats.npts) )### for unknown reasons, this has to be done, nsam=int(tempsac.stats.npts)  won't work as an input for aftan
        dt=tempsac.stats.delta
        try:
            dist=tempsac.stats.sac.dist;
        except:
            dist=tempsac.stats.distance;
        # Start to do aftan utilizing pyaftan
        self.ftanparam.nfout1_1,self.ftanparam.arr1_1,self.ftanparam.nfout2_1,self.ftanparam.arr2_1,self.ftanparam.tamp_1,\
        self.ftanparam.nrow_1,self.ftanparam.ncol_1,self.ftanparam.ampo_1, self.ftanparam.ierr_1= ftan.aftanpg(piover4, nsam, \
        sig, tb, dt, dist, vmin, vmax, tmin, tmax, tresh, ffact, perc, npoints, taperl, nfin, snr, nprpv, phprper, phprvel)
        if pmf==True:
            if self.ftanparam.nfout2_1<3:
                return
            npred = self.ftanparam.nfout2_1
            tmin2 = self.ftanparam.arr2_1[1,0]
            tmax2 = self.ftanparam.arr2_1[1,self.ftanparam.nfout2_1-1]
            pred=np.zeros((2,300))
            pred[:,0:100]=self.ftanparam.arr2_1[1:3,:]
            pred=pred.T
            self.ftanparam.nfout1_2,self.ftanparam.arr1_2,self.ftanparam.nfout2_2,self.ftanparam.arr2_2,self.ftanparam.tamp_2, \
            self.ftanparam.nrow_2,self.ftanparam.ncol_2,self.ftanparam.ampo_2, self.ftanparam.ierr_2=ftan.aftanipg(piover4,nsam, \
            sig,tb,dt,dist,vmin,vmax,tmin2,tmax2,tresh,ffact,perc,npoints,taperl,nfin,snr,fmatch,npred,pred,nprpv,phprper,phprvel)
        return

    def findcoda(self, tfactor=2, Lwin=1200, Tmin=5.0, Tmax=10.0, method='stehly' ):
        """findcoda:
        Find the coda for a given trace, return both negative and positive trace and begin/end time
        This function find the coda window for a given noisetrace, the arrival of surface wave package is determined by aftan analysis.
        Two methods are utilized to find the coda.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        tfactor - 
            method='stehly'
            begin time is the arrival time of surface wave multiplied by tfactor
            method='ma'
            begin time is the arrival time of surface wave plus tfactor( in sec)
        Lwin    - coda window length
        Tmin    - minimum period for FTAN analysis
        Tmax    - maximum period for FTAN analysis
        method  - 'stehly' or 'ma'
        
        Output:
        neg     - negative lag data (numpy array)
        pos     - positive lag data (numpy array)
        Tbeg    - begin time of coda window
        Tend    - end time of coda window
        -----------------------------------------------------------------------------------------------------
        References:
        Stehly, L., et al. Reconstructing Green's function by correlation of the coda of the correlation (C3) of ambient seismic noise.
             Journal of Geophysical Research: Solid Earth (1978?2012) 113.B11 (2008).
        Ma, Shuo, and Gregory C. Beroza. Ambient field Green's functions from asynchronous seismic observations. Geophysical Research Letters 39.6 (2012).
        """
        npts=self.stats.npts
        Nhalf=(npts-1)/2+1
        neg=self.data[:Nhalf]
        pos=self.data[Nhalf-1:]
        seis=self.copy()
        seis.aftan(tmin=Tmin, tmax=Tmax)
        self.ftanparam=seis.ftanparam
        vposition=np.where(np.logical_and(self.ftanparam.arr2_1[1,:]>Tmin, self.ftanparam.arr2_1[1,:]<Tmax))[0]
        if vposition.size==0:
            vposition=np.where(np.logical_and(self.ftanparam.arr1_1[1,:]>Tmin, self.ftanparam.arr1_1[1,:]<Tmax))[0]
            Vmean=self.ftanparam.arr1_1[2,vposition].mean()
            print "Warning: no jump corrected result for: "+self.stats.sac.kevnm.split('\x00')[0]+"_"+self.stats.station+" "+self.stats.channel+' V:' +str(Vmean)
        else:
            Vmean=self.ftanparam.arr2_1[2,vposition].mean()

        if method=='ma':
            Tbeg=max( (int)( (self.stats.sac.dist/Vmean)/self.stats.delta ) +(int)(tfactor/self.stats.delta) , 1)
            Tend=Tbeg+(int)(Lwin/self.stats.delta)
        else:
            Tbeg=(int)( (self.stats.sac.dist/Vmean)/self.stats.delta * tfactor ) +1
            Tend=Tbeg+(int)(Lwin/self.stats.delta)
        if Tend>Nhalf:
            Tend=Nhalf
            print "Warning: The default coda window end excess the record length!"
        neg=neg[::-1]
        neg=neg[Tbeg-1:Tend]
        pos=pos[Tbeg-1:Tend]
        return neg, pos, Tbeg, Tend

    def findcodaTime(self, tfactor=2, Lwin=1200, Tmin=5.0, Tmax=10.0, method='stehly' ):
        """findcodaTime:
        Simlar to findcoda but only return begin/end time.
        """
        npts=self.stats.npts
        Nhalf=(npts-1)/2+1
        seis=self.copy()
        seis.aftan(tmin=Tmin, tmax=Tmax)
        self.ftanparam=seis.ftanparam
        vposition=np.where(np.logical_and(self.ftanparam.arr2_1[1,:]>Tmin, self.ftanparam.arr2_1[1,:]<Tmax))[0]
        if vposition.size==0:
            vposition=np.where(np.logical_and(self.ftanparam.arr1_1[1,:]>Tmin, self.ftanparam.arr1_1[1,:]<Tmax))[0]
            Vmean=self.ftanparam.arr1_1[2,vposition].mean()
            print "Warning: no jump corrected result for: "+self.stats.sac.kevnm.split('\x00')[0]+"_"+self.stats.station+" "+self.stats.channel+' V:' +str(Vmean)
        else:
            Vmean=self.ftanparam.arr2_1[2,vposition].mean()
        if Vmean<0.5 or Vmean >6 or np.isnan(Vmean):
            return -1, -1
        if method=='ma':
            Tbeg=max( (int)( (self.stats.sac.dist/Vmean)/self.stats.delta ) +(int)(tfactor/self.stats.delta) , 1)
            Tend=Tbeg+(int)(Lwin/self.stats.delta)
        else:
            Tbeg=(int)( (self.stats.sac.dist/Vmean)/self.stats.delta * tfactor ) +1
            Tend=Tbeg+(int)(Lwin/self.stats.delta)
        if Tend>Nhalf:
            Tend=-1
            print "Warning: The default coda window end excess the record length!"
        return Tbeg, Tend;

    def findcommoncoda(self, Tbeg1, Tbeg2, Tend1, Tend2):
        """findcommoncoda:
        Return common coda trace for two noisetraces
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        Tbeg1, Tend1 - coda window time for the first trace
        Tbeg1, Tend1 - coda window time for the second trace
        
        Output:
        neg     - common negative lag data (numpy array)
        pos     - common positive lag data (numpy array)
        Tbeg    - begin time of common coda window
        Tend    - end time of common coda window
        -----------------------------------------------------------------------------------------------------
        """
        npts=self.stats.npts
        Nhalf=(npts-1)/2+1
        neg=self.data[:Nhalf]
        pos=self.data[Nhalf-1:]
        if Tbeg1>Tbeg2:
            Tbeg=Tbeg1
            Tend=Tend1
        else:
            Tbeg=Tbeg2
            Tend=Tend2
        neg=neg[::-1]
        neg=neg[Tbeg-1:Tend]
        pos=pos[Tbeg-1:Tend]
        return neg, pos, Tbeg, Tend;

    def plotftan(self, plotflag=3, sacname=''):
        """
        Plot ftan diagram:
        This function plot ftan diagram.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        plotflag -
            0: only Basic FTAN
            1: only Phase Matched Filtered FTAN
            2: both
            3: both in one figure
        sacname - sac file name than can be used as the title of the figure
        -----------------------------------------------------------------------------------------------------
        """
        try:
            fparam=self.ftanparam
            if fparam.nfout1_1==0:
                return "Error: No Basic FTAN parameters!"
            dt=self.stats.delta
            dist=self.stats.sac.dist
            if (plotflag!=1 and plotflag!=3):
                v1=dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
                ampo_1=fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
                obper1_1=fparam.arr1_1[1,:fparam.nfout1_1]
                gvel1_1=fparam.arr1_1[2,:fparam.nfout1_1]
                phvel1_1=fparam.arr1_1[3,:fparam.nfout1_1]
                plb.figure()
                ax = plt.subplot()
                p=plt.pcolormesh(obper1_1, v1, ampo_1, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_1, gvel1_1, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_1, phvel1_1, '--w', lw=3) #

                if (fparam.nfout2_1!=0):
                    obper2_1=fparam.arr2_1[1,:fparam.nfout2_1]
                    gvel2_1=fparam.arr2_1[2,:fparam.nfout2_1]
                    phvel2_1=fparam.arr2_1[3,:fparam.nfout2_1]
                    ax.plot(obper2_1, gvel2_1, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_1, phvel2_1, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin1=obper1_1[0]
                Tmax1=obper1_1[fparam.nfout1_1-1]
                vmin1= v1[fparam.ncol_1-1]
                vmax1=v1[0]
                plt.axis([Tmin1, Tmax1, vmin1, vmax1])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('Basic FTAN Diagram '+sacname,fontsize=15)

            if fparam.nfout1_2==0 and plotflag!=0:
                return "Error: No PMF FTAN parameters!"
            if (plotflag!=0 and plotflag!=3):
                v2=dist/(fparam.tamp_2+np.arange(fparam.ncol_2)*dt)
                ampo_2=fparam.ampo_2[:fparam.ncol_2,:fparam.nrow_2]
                obper1_2=fparam.arr1_2[1,:fparam.nfout1_2]
                gvel1_2=fparam.arr1_2[2,:fparam.nfout1_2]
                phvel1_2=fparam.arr1_2[3,:fparam.nfout1_2]
                plb.figure()
                ax = plt.subplot()
                p=plt.pcolormesh(obper1_2, v2, ampo_2, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_2, gvel1_2, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_2, phvel1_2, '--w', lw=3) #

                if (fparam.nfout2_2!=0):
                    obper2_2=fparam.arr2_2[1,:fparam.nfout2_2]
                    gvel2_2=fparam.arr2_2[2,:fparam.nfout2_2]
                    phvel2_2=fparam.arr2_2[3,:fparam.nfout2_2]
                    ax.plot(obper2_2, gvel2_2, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_2, phvel2_2, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin2=obper1_2[0]
                Tmax2=obper1_2[fparam.nfout1_2-1]
                vmin2= v2[fparam.ncol_2-1]
                vmax2=v2[0]
                plt.axis([Tmin2, Tmax2, vmin2, vmax2])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('PMF FTAN Diagram '+sacname,fontsize=15)

            if ( plotflag==3 ):
                v1=dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
                ampo_1=fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
                obper1_1=fparam.arr1_1[1,:fparam.nfout1_1]
                gvel1_1=fparam.arr1_1[2,:fparam.nfout1_1]
                phvel1_1=fparam.arr1_1[3,:fparam.nfout1_1]
                plb.figure(num=None, figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
                ax = plt.subplot(2,1,1)
                p=plt.pcolormesh(obper1_1, v1, ampo_1, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_1, gvel1_1, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_1, phvel1_1, '--w', lw=3) #
                if (fparam.nfout2_1!=0):
                    obper2_1=fparam.arr2_1[1,:fparam.nfout2_1]
                    gvel2_1=fparam.arr2_1[2,:fparam.nfout2_1]
                    phvel2_1=fparam.arr2_1[3,:fparam.nfout2_1]
                    ax.plot(obper2_1, gvel2_1, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_1, phvel2_1, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin1=obper1_1[0]
                Tmax1=obper1_1[fparam.nfout1_1-1]
                vmin1= v1[fparam.ncol_1-1]
                vmax1=v1[0]
                plt.axis([Tmin1, Tmax1, vmin1, vmax1])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('Basic FTAN Diagram '+sacname)

                v2=dist/(fparam.tamp_2+np.arange(fparam.ncol_2)*dt)
                ampo_2=fparam.ampo_2[:fparam.ncol_2,:fparam.nrow_2]
                obper1_2=fparam.arr1_2[1,:fparam.nfout1_2]
                gvel1_2=fparam.arr1_2[2,:fparam.nfout1_2]
                phvel1_2=fparam.arr1_2[3,:fparam.nfout1_2]

                ax = plt.subplot(2,1,2)
                p=plt.pcolormesh(obper1_2, v2, ampo_2, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_2, gvel1_2, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_2, phvel1_2, '--w', lw=3) #

                if (fparam.nfout2_2!=0):
                    obper2_2=fparam.arr2_2[1,:fparam.nfout2_2]
                    gvel2_2=fparam.arr2_2[2,:fparam.nfout2_2]
                    phvel2_2=fparam.arr2_2[3,:fparam.nfout2_2]
                    ax.plot(obper2_2, gvel2_2, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_2, phvel2_2, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin2=obper1_2[0]
                Tmax2=obper1_2[fparam.nfout1_2-1]
                vmin2= v2[fparam.ncol_2-1]
                vmax2=v2[0]
                plt.axis([Tmin2, Tmax2, vmin2, vmax2])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('PMF FTAN Diagram '+sacname)
        except AttributeError:
            print 'Error: FTAN Parameters are not available!'
        return

    def GaussianFilter(self, fcenter, fhlen=0.008):
        """
        Gaussian Filter designed for SNR analysis, utilize pyfftw to do fft
        exp( (-0.5/fhlen^2)*(f-fcenter)^2 )
        -----------------------------------------------------------------------------------------------------
        Input parameters:
        fcenter - central period
        fhlen   - half length of Gaussian width
        -----------------------------------------------------------------------------------------------------
        """
        npts=self.stats.npts
        Ns=1<<(npts-1).bit_length()
        df=1.0/self.stats.delta/Ns
        nhalf=Ns/2+1
        fmax=(nhalf-1)*df
        if fcenter>fmax:
            fcenter=fmax
        alpha = -0.5/(fhlen*fhlen)
        F=np.arange(Ns)*df
        gauamp = F - fcenter
        sf=npr.evaluate('exp(alpha*gauamp**2)')
        sp, Ns=FFTW(self.data, direction='forward')
        filtered_sp=npr.evaluate('sf*sp')
        filtered_seis, Ns=FFTW(filtered_sp, direction='backward')
        filtered_seis=filtered_seis[:npts].real
        return filtered_seis

    def getSNR(self, foutPR='', fhlen=0.008, pmf=True, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, \
        tmax=30.0, tresh=20.0, ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0,phvelname=''):
        """getSNR
        Get the SNR for signal window based on FTAN analysis.
        If input noisetrace is double-lagged, it will do SNR analysis for pos/neg lag; otherwise it will do SNR analysis for sym lag.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        foutPR      - Output file prefix for positive and negative lags, default is ''(NOT to save positive/negative files)
        fhlen       - half length of Gaussian width
        pmf         - flag for Phase-Matched-Filtered output (default: True)
        piover4     - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin        - minimal group velocity, km/s
        vmax        - maximal group velocity, km/s
        tmin        - minimal period, s
        tmax        - maximal period, s
        tresh       - treshold for jump detection, usualy = 10, need modifications
        ffact       - factor to automatic filter parameter, usualy =1
        taperl      - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        snr         - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
        fmatch      - factor to length of phase matching window
        fname       - SAC file name
        phvelname   - predicted phase velocity file name
        
        Output:
        self.SNRParam, a object of snrParam class to store output SNR results
        -----------------------------------------------------------------------------------------------------
        References:
        Tian, Ye, and Michael H. Ritzwoller. "Directionality of ambient noise on the Juan de Fuca plate:
            implications for source locations of the primary and secondary microseisms." Geophysical Journal International 201.1 (2015): 429-443.
        """
        try:
            self.SNRParam
        except:
            self.init_snrParam()
        try:
            dist=self.stats.sac.dist
            begT=self.stats.sac.b
            endT=self.stats.sac.e
            symFlag=False
            dt=self.stats.delta
            if ( abs(begT+endT) < dt ):
                symFlag=True
                begT=0.
                temp_pos=self.getpos()
                temp_neg=self.getneg()
                temp_pos.aftan(pmf=pmf, piover4=piover4, vmin=vmin, vmax=vmax, tmin=tmin, \
                    tmax=tmax, tresh=tresh, ffact=ffact, taperl=taperl, snr=snr, fmatch=fmatch,phvelname=phvelname)
                temp_neg.aftan(pmf=pmf, piover4=piover4, vmin=vmin, vmax=vmax, tmin=tmin, \
                    tmax=tmax, tresh=tresh, ffact=ffact, taperl=taperl, snr=snr, fmatch=fmatch,phvelname=phvelname)
                fparam_p=temp_pos.ftanparam
                fparam_n=temp_neg.ftanparam
                if foutPR!='':
                    temp_pos.ftanparam.writeDISP(foutPR+'_pos')
                    temp_neg.ftanparam.writeDISP(foutPR+'_neg')
            if symFlag==True:
                if fparam_p.nfout2_2!=0:
                    o_per_p=fparam_p.arr2_2[1,:]
                    g_vel_p=fparam_p.arr2_2[2,:]
                    self.SNRParam.oper_p=o_per_p
                ### positive lag
                for i in np.arange(fparam_p.nfout2_2):
                    filtered_tr=temp_pos.GaussianFilter(1./o_per_p[i], fhlen=fhlen)
                    minT = dist/g_vel_p[i]-o_per_p[i]/2.
                    maxT = dist/g_vel_p[i]+o_per_p[i]/2.
                    if g_vel_p[i]<0 or o_per_p[i]<0:
                         self.SNRParam.snr_p=np.append(self.SNRParam.snr_p, -1.);
                         self.SNRParam.amp_p=np.append(self.SNRParam.amp_p, -1.);
                         continue;
                    if(minT<begT):
                        minT=begT
                    if(maxT>endT):
                        maxT=endT
                    ib = (int)(minT/dt)
                    ie = (int)(maxT/dt)+2
                    tempTr_p=filtered_tr[ib:ie]
                    tempTr_p=npr.evaluate('abs(tempTr_p)')
                    tempmax_p=tempTr_p.max()
                    self.SNRParam.amp_p=np.append(self.SNRParam.amp_p, tempmax_p)
                    # Noise window
                    minT = maxT + o_per_p[i] * 5 + 500.
                    skipflag=False
                    if( (endT - minT) < 50. ):
                        self.SNRParam.snr_p=np.append(self.SNRParam.snr_p, -1.)
                        skipflag=True
                    elif( (endT - minT) < 1100. ):
                        maxT = endT - 10.
                    else:
                        minT = endT - 1100.
                        maxT = endT - 100.
                    if skipflag==False:
                        ib = (int)(minT/dt)
                        ie = (int)(maxT/dt)+2
                        tempnoise_p=filtered_tr[ib:ie]
                        tempnoise_p=npr.evaluate('tempnoise_p**2')
                        noiserms_p=math.sqrt(npr.evaluate('sum(tempnoise_p)')/(ie-ib-1.))
                        self.SNRParam.nrms_p=np.append(self.SNRParam.nrms_p, noiserms_p)
                        tempSNR_p=tempmax_p/noiserms_p
                        self.SNRParam.snr_p=np.append(self.SNRParam.snr_p, tempSNR_p)
                ### negative lag
                if fparam_n.nfout2_2!=0:
                    o_per_n=fparam_n.arr2_2[1,:]
                    g_vel_n=fparam_n.arr2_2[2,:]
                    self.SNRParam.oper_n=o_per_n
                for i in np.arange(fparam_n.nfout2_2):
                    minT = dist/g_vel_n[i]-o_per_n[i]/2.
                    maxT = dist/g_vel_n[i]+o_per_n[i]/2.
                    if g_vel_n[i]<0 or o_per_n[i]<0:
                         self.SNRParam.snr_n=np.append(self.SNRParam.snr_n, -1.);
                         self.SNRParam.amp_n=np.append(self.SNRParam.amp_n, -1.);
                         continue;
                    filtered_tr=temp_neg.GaussianFilter(1./o_per_n[i], fhlen=fhlen)
                    if(minT<begT):
                        minT=begT
                    if(maxT>endT):
                        maxT=endT
                    ib = (int)(minT/dt)
                    ie = (int)(maxT/dt)+2
                    # print ib,ie, minT, maxT, g_vel_n[i], o_per_n[i]
                    tempTr_n=filtered_tr[ib:ie]
                    tempTr_n=npr.evaluate('abs(tempTr_n)')
                    tempmax_n=tempTr_n.max()
                    self.SNRParam.amp_n=np.append(self.SNRParam.amp_n, tempmax_n)
                    # Noise window
                    minT = maxT + o_per_n[i] * 5 + 500.
                    skipflag=False
                    if( (endT - minT) < 50. ):
                        self.SNRParam.snr_n=np.append(self.SNRParam.snr_n, -1.)
                        skipflag=True
                    elif( (endT - minT) < 1100. ):
                        maxT = endT - 10.
                    else:
                        minT = endT - 1100.
                        maxT = endT - 100.
                    if skipflag==False:
                        ib = (int)(minT/dt)
                        ie = (int)(maxT/dt)+2
                        tempnoise_n=filtered_tr[ib:ie]
                        tempnoise_n=npr.evaluate('tempnoise_n**2')
                        noiserms_n=math.sqrt(npr.evaluate('sum(tempnoise_n)')/(ie-ib-1.))
                        self.SNRParam.nrms_n=np.append(self.SNRParam.nrms_n, noiserms_n)
                        tempSNR_n=tempmax_n/noiserms_n
                        self.SNRParam.snr_n=np.append(self.SNRParam.snr_n, tempSNR_n)
            else:
                fparam=self.ftanparam
                if fparam.nfout2_2!=0:
                    o_per=fparam.arr2_2[1,:]
                    g_vel=fparam.arr2_2[2,:]
                    self.SNRParam.oper_s=o_per
                for i in np.arange(fparam.nfout2_2):
                    filtered_tr=self.GaussianFilter(1./o_per[i], fhlen=fhlen)
                    minT = dist/g_vel[i]-o_per[i]/2.
                    maxT = dist/g_vel[i]+o_per[i]/2.
                    if g_vel[i]<0 or o_per[i]<0:
                         self.SNRParam.snr_s=np.append(self.SNRParam.snr_s, -1.);
                         self.SNRParam.amp_s=np.append(self.SNRParam.amp_s, -1.)
                         continue;
                    if(minT<begT):
                        minT=begT
                    if(maxT>endT):
                        maxT=endT
                    ib = (int)(minT/dt)
                    ie = (int)(maxT/dt)+2
                    tempTr_s=filtered_tr[ib:ie]
                    tempTr_s=npr.evaluate('abs(tempTr_s)')
                    tempmax_s=tempTr_s.max()
                    self.SNRParam.amp_s=np.append(self.SNRParam.amp_s, tempmax_s)
                    # Noise window
                    minT = maxT + o_per[i] * 5 + 500.
                    skipflag=False
                    if( (endT - minT) < 50. ):
                        self.SNRParam.snr_s=np.append(self.SNRParam.snr_s, -1.)
                        skipflag=True
                    elif( (endT - minT) < 1100. ):
                        maxT = endT - 10.
                    else:
                        minT = endT - 1100.
                        maxT = endT - 100.
                    if skipflag==False:
                        ib = (int)(minT/dt)
                        ie = (int)(maxT/dt)+2
                        tempnoise_s=filtered_tr[ib:ie]
                        tempnoise_s=npr.evaluate('tempnoise_s**2')
                        noiserms_s=math.sqrt(npr.evaluate('sum(tempnoise_s)')/(ie-ib-1.))
                        self.SNRParam.nrms_s=np.append(self.SNRParam.nrms_s, noiserms_s)
                        tempSNR_s=tempmax_s/noiserms_s
                        self.SNRParam.snr_s=np.append(self.SNRParam.snr_s, tempSNR_s)
        except AttributeError:
            print 'Error: FTAN Parameters are not available!'
        return

def FFTW(indata, direction, flags=['estimate']):
    """
    FFTW: a function utilizes fftw, a extremely fast library to do FFT computation (pyfftw3 need to be installed)
    -----------------------------------------------------------------------------------------------------
    Input Parameters:
    indata      - Input data
    direction   - direction of FFT
    flags       - list of fftw-flags to be used in planning
    -----------------------------------------------------------------------------------------------------
    Functions that using this function:
        noisetrace.GaussianFilter()
    """
    npts=indata.size
    Ns=1<<(npts-1).bit_length()
    INput = np.zeros((Ns), dtype=complex)
    OUTput = np.zeros((Ns), dtype=complex)
    fftw = fftw3.Plan(INput, OUTput, direction=direction, flags=flags)
    INput[:npts]=indata
    fftw()
    nhalf=Ns/2+1
    if direction == 'forward':
        OUTput[nhalf:]=0
        OUTput[0]/=2
        OUTput[nhalf-1]=OUTput[nhalf-1].real+0.j
    if direction =='backward':
        OUTput=2*OUTput/Ns
    return OUTput, Ns


f='/projects/life9360/code/CPSPy/sac_dir/B00511ZVF.sac'
st=obspy.read(f)
tr=st[0]
tr1=noisetrace(tr.data, tr.stats)
tr1.aftan(piover4=-1., pmf=True, tmin=2.0, tmax=50.0)

ds = ASDFDataSet("./observed.h5")
ds.add_auxiliary_data(tr1.ftanparam.arr1_1, data_type='DISP',path='AAA', parameters={"a": 1, "b": 2})
# ds.add_auxiliary_data_file('/lustre/janus_scratch/life9360/specfem2d_data/aftan_homo/MEM2D.45S37..SAC_1_DISP.1', path='/DISP')