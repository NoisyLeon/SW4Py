# -*- coding: utf-8 -*-

import numpy as np
import pyasdf 
import matplotlib.pyplot as plt
import stations
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
# import multiprocessing as mp
import math
import time
import shutil
from subprocess import call
import warnings


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
        self.station_id=None

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
            Lf10=self.nfout1_1;
            outArrf10=np.arange(Lf10);
            for i in np.arange(7):
                outArrf10=np.append(outArrf10, self.arr1_1[i,:Lf10]);
            outArrf10=outArrf10.reshape((8,Lf10));
            outArrf10=outArrf10.T;
            np.savetxt(f10, outArrf10, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
        if self.nfout2_1!=0:
            f11=fnamePR+'_1_DISP.1'
            Lf11=self.nfout2_1;
            outArrf11=np.arange(Lf11);
            for i in np.arange(6):
                outArrf11=np.append(outArrf11, self.arr2_1[i,:Lf11]);
            outArrf11=outArrf11.reshape((7,Lf11));
            outArrf11=outArrf11.T;
            np.savetxt(f11, outArrf11, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
        if self.nfout1_2!=0:
            f20=fnamePR+'_2_DISP.0';
            Lf20=self.nfout1_2;
            outArrf20=np.arange(Lf20);
            for i in np.arange(7):
                outArrf20=np.append(outArrf20, self.arr1_2[i,:Lf20]);
            outArrf20=outArrf20.reshape((8,Lf20));
            outArrf20=outArrf20.T;
            np.savetxt(f20, outArrf20, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
        if self.nfout2_2!=0:
            f21=fnamePR+'_2_DISP.1';
            Lf21=self.nfout2_2;
            outArrf21=np.arange(Lf21);
            for i in np.arange(6):
                outArrf21=np.append(outArrf21, self.arr2_2[i,:Lf21]);
            outArrf21=outArrf21.reshape((7,Lf21));
            outArrf21=outArrf21.T;
            np.savetxt(f21, outArrf21, fmt='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf');
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

class ftanLst(object):
    """
    An object contains a ftanparam list(a list of ftanParam object) information several methods for ftanparam list related analysis.
        ftanparams: list of ftanParam
    """
    def __init__(self,ftanparams=None):
        self.ftanparams=[]
        if isinstance(ftanparams, ftanParam):
            ftanparams = [ftanparams]
        if ftanparams:
            self.ftanparams.extend(ftanparams)

    def __add__(self, other):
        """
        Add two ftanLst with self += other.
        """
        if isinstance(other, ftanParam):
            other = ftanLst([other])
        if not isinstance(other, ftanLst):
            raise TypeError
        ftanparams = self.ftanparams + other.ftanparams
        return self.__class__(ftanparams=ftanparams)

    def __len__(self):
        """
        Return the number of Traces in the ftanLst object.
        """
        return len(self.ftanparams)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(ftanparams=self.ftanparams.__getitem__(index))
        else:
            return self.ftanparams.__getitem__(index)

    def append(self, ftanparam):
        """
        Append a single ftanParam object to the current ftanLst object.
        """
        if isinstance(ftanparam, ftanParam):
            self.ftanparams.append(ftanparam)
        else:
            msg = 'Append only supports a single ftanParam object as an argument.'
            raise TypeError(msg)
        return self


class sw4trace(obspy.core.trace.Trace):
    """
    sw4trace:
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
    def aftan(self, pmf=True, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, \
        tmax=30.0, tresh=20.0, ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0, phvelname='', predV=np.array([]) ):

        """ (Automatic Frequency-Time ANalysis) aftan analysis:
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
        predV          - predicted phase velocity curve
        
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
            dist=self.stats.sac.dist
        except:
            dist = np.sqrt( (self.stats.sac.evlo - self.stats.sac.stlo)**2 + (self.stats.sac.evla - self.stats.sac.stla)**2 )
            self.stats.sac.dist=dist
        if (phvelname==''):
            phvelname='./ak135.disp'
        nprpv = 0
        phprper=np.zeros(300)
        phprvel=np.zeros(300)
        if predV.size != 0:
            phprper=predV[:,0]
            phprvel=predV[:,1]
            nprpv = predV[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        elif os.path.isfile(phvelname):
            php=np.loadtxt(phvelname)
            phprper=php[:,0]
            phprvel=php[:,1]
            nprpv = php[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        nfin = 64
        npoints = 5  #  only 3 points in jump
        perc    = 50.0 # 50 % for output segment
        tempsac=self.copy()
        tb=self.stats.sac.b
        length=len(tempsac.data)
        if length>32768:
            print "Warning: length of seismogram is larger than 32768!"
            nsam=32768
            tempsac.data=tempsac.data[:nsam]
            tempsac.stats.e=(nsam-1)*tempsac.stats.delta+tb
            sig=tempsac.data
        else:
            sig=np.append(tempsac.data, np.zeros( float(32768-tempsac.data.size) ) )
            nsam=int( float (tempsac.stats.npts) )### for unknown reasons, this has to be done, nsam=int(tempsac.stats.npts)  won't work as an input for aftan
        dt=tempsac.stats.delta
        # Start to do aftan utilizing pyaftan
        self.ftanparam.nfout1_1,self.ftanparam.arr1_1,self.ftanparam.nfout2_1,self.ftanparam.arr2_1,self.ftanparam.tamp_1, \
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
                    self.ftanparam.nrow_2,self.ftanparam.ncol_2,self.ftanparam.ampo_2, self.ftanparam.ierr_2 = ftan.aftanipg(piover4,nsam, \
                        sig,tb,dt,dist,vmin,vmax,tmin2,tmax2,tresh,ffact,perc,npoints,taperl,nfin,snr,fmatch,npred,pred,nprpv,phprper,phprvel)
        self.ftanparam.station_id=self.stats.network+'.'+self.stats.station 
        return

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



class InputFtanParam(object): ###
    """
    A subclass to store input parameters for aftan analysis and SNR Analysis
    -----------------------------------------------------------------------------------------------------
    Parameters:
    pmf         - flag for Phase-Matched-Filtered output (default: Fasle)
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
    fhlen       - half length of Gaussian width
    dosnrflag   - whether to do SNR analysis or not
    predV          - predicted phase velocity curve
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self):
        self.pmf=False
        self.piover4=-1.0
        self.vmin=1.5
        self.vmax=5.0
        self.tmin=4.0
        self.tmax=30.0
        self.tresh=20.0
        self.ffact=1.0
        self.taperl=1.0
        self.snr=0.2
        self.fmatch=1.0
        self.fhlen=0.008
        self.dosnrflag=False
        self.predV=np.array([])

    def setInParam(self, pmf=False, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, \
        tmax=30.0, tresh=20.0, ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0, fhlen=0.008, dosnrflag=False, predV=np.array([]) ):
        """
        Set the parameters
        """
        self.pmf=pmf
        self.piover4=piover4
        self.vmin=vmin
        self.vmax=vmax
        self.tmin=tmin
        self.tmax=tmax
        self.tresh=tresh
        self.ffact=ffact
        self.taperl=taperl
        self.snr=snr
        self.fmatch=fmatch
        self.fhlen=fhlen
        self.dosnrflag=dosnrflag
        self.predV=predV
        return

class sw4ASDF(pyasdf.ASDFDataSet):
    
    def Readsac(self, stafile, datadir, comptype='u', datatype='displacement'):
        if comptype == 'all':
            comptype=['e', 'n', 'u']
        else:
            comptype=[comptype]
        SLst=stations.StaLst()
        SLst.ReadStaList(stafile=stafile)
        StaInv=SLst.GetInventory()
        self.add_stationxml(StaInv)
        for sta in SLst.stations:
            if sta.variables == 'displacement':
                sacsfx=''
            else:
                sacsfx='v'
            for comp in comptype:
                sacfname = datadir+'/'+sta.network+'.'+sta.stacode+'.'+comp+sacsfx
                tr=obspy.read(sacfname)[0]
                # print tr.stats.sac.o
                tr.stats.network=sta.network
                self.add_waveforms(tr, tag='sw4_raw')
        return
    
    def AddEvent(self, x, y, z):
        print 'Attention: Event Location unit is km!'
        origin=obspy.core.event.origin.Origin(longitude=x, latitude=y, depth=z)
        event=obspy.core.event.event.Event(origins=[origin])
        catalog=obspy.core.event.catalog.Catalog(events=[event])
        self.add_quakeml(catalog)
        return
    
    
    def aftan(self, compindex=0, tb=0., outdir=None, inftan=InputFtanParam(), phvelname ='./ak135.disp', basic1=True, basic2=False,
            pmf1=False, pmf2=False):
        """
        _1_DISP.0: arr1_1
        _1_DISP.1: arr2_1
        _2_DISP.0: arr1_2
        _2_DISP.1: arr2_2
        """
        # compdict={'e':0, 'n':1, 'u':2}
        evlo=self.events.events[0].origins[0].longitude
        evla=self.events.events[0].origins[0].latitude
        predV=np.loadtxt(phvelname)
        for stacode in self.waveforms.list():
            # Get data from ASDF dataset
            tr=self.waveforms[stacode].sw4_raw[compindex]
            tr.stats.sac={}
            tr.stats.sac.evlo=evlo
            tr.stats.sac.evla=evla
            tr.stats.sac.b=tb
            stlo=self.waveforms[stacode].coordinates['longitude']
            stla=self.waveforms[stacode].coordinates['latitude']
            tr.stats.sac.stlo=stlo*100.
            tr.stats.sac.stla=stla*100.
            # aftan analysis
            ntrace=sw4trace(tr.data, tr.stats)
            ntrace.aftan(pmf=inftan.pmf, piover4=inftan.piover4, vmin=inftan.vmin,
                vmax=inftan.vmax, tmin=inftan.tmin, tmax=inftan.tmax, tresh=inftan.tresh,
                ffact=inftan.ffact, taperl=inftan.taperl, snr=inftan.snr, fmatch=inftan.fmatch, predV=inftan.predV)
            print 'aftan analysis for', tr.stats.station, ntrace.stats.sac.dist
            # save aftan results to ASDF dataset
            if basic1==True:
                parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'dis': 5, 'snrdb': 6, 'mhw': 7, 'amp': 8, 'Np': ntrace.ftanparam.nfout1_1}
                self.add_auxiliary_data(data=ntrace.ftanparam.arr1_1, data_type='DISPbasic1', path=tr.stats.station, parameters=parameters)
            if basic2==True:
                parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'snrdb': 5, 'mhw': 6, 'amp': 7, 'Np': ntrace.ftanparam.nfout2_1}
                self.add_auxiliary_data(data=ntrace.ftanparam.arr2_1, data_type='DISPbasic2', path=tr.stats.station, parameters=parameters)
            if inftan.pmf==True:
                if pmf1==True:
                    parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'dis': 5, 'snrdb': 6, 'mhw': 7, 'amp': 8, 'Np': ntrace.ftanparam.nfout1_2}
                    self.add_auxiliary_data(data=ntrace.ftanparam.arr1_2, data_type='DISPpmf1', path=tr.stats.station, parameters=parameters)
                if pmf2==True:
                    parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'dis': 5, 'snrdb': 6, 'mhw': 7, 'amp': 8, 'Np': ntrace.ftanparam.nfout2_2}
                    self.add_auxiliary_data(data=ntrace.ftanparam.arr2_2, data_type='DISPpmf2', path=tr.stats.station, parameters=parameters)
        ### dbase.auxiliary_data.DISPbasic1['112S1000'].data.value[dbase.auxiliary_data.DISPbasic1['112S1000'].parameters['Vph']]
        return        
        
        
                
            
        
        
        
    
            
        
        
    
        
    





