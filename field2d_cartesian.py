#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sinter
from matplotlib.mlab import griddata
import numpy.ma as ma
import scipy.ndimage.filters 
from math import pi
# from skimage.filters import roberts, sobel, scharr, prewitt
from scipy.ndimage import convolve

X_diff_weight_2 = np.array([[1., 0., -1.]])/2.
Y_diff_weight_2 = X_diff_weight_2.T
X_diff_weight_4 = np.array([[-1., 8., 0., -8., 1.]])/12.
Y_diff_weight_4 = X_diff_weight_4.T
X_diff_weight_6 = np.array([[1./60., 	-3./20.,  3./4.,  0., -3./4., 3./20.,  -1./60.]])
Y_diff_weight_6 = X_diff_weight_6.T

X_diff2_weight_2 = np.array([[1., -2., 1.]])
Y_diff2_weight_2 = X_diff2_weight_2.T
X_diff2_weight_4 = np.array([[-1., 16., -30., 16., -1.]])/12.
Y_diff2_weight_4 = X_diff2_weight_4.T
X_diff2_weight_6 = np.array([[1./90., 	-3./20.,  3./2.,  -49./18., 3./2., -3./20.,  1./90.]])
Y_diff2_weight_6 = X_diff2_weight_6.T


class Field2d(object):
    def __init__(self, Nx, Ny, dx, dy):
        self.Nx=int(Nx)+1
        self.Ny=int(Ny)+1
        self.dx=dx
        self.dy=dy
        self.x=np.arange(Nx+1)*self.dx
        self.y=np.arange(Ny+1)*self.dy
        self.Xarr, self.Yarr = np.meshgrid(self.x, self.y)
        # self.Xarr=self.Xarr.reshape(Nx*Ny)
        # self.Yarr=self.Yarr.reshape(Nx*Ny)
        return
    
    def LoadFile(self, fname):
        try:
            Inarray=np.loadtxt(fname)
        except:
            Inarray=np.load(fname)
        self.XarrIn=Inarray[:,0]
        self.YarrIn=Inarray[:,1]
        self.ZarrIn=Inarray[:,2]
        try:
            self.DistarrIn=Inarray[:,3]
        except:
            print 'No distance information in the field file!'
        return
    
    def LoadField(self, inField):
        self.XarrIn=inField.Xarr
        self.YarrIn=inField.Yarr
        self.ZarrIn=inField.Zarr
        try:
            self.DistarrIn=inField.Distarr
        except:
            print 'No distance information in the input field!'
        return
    
    def SaveFile(self, fname, fmt='npy'):
        OutArr=np.append(self.Xarr, self.Yarr)
        OutArr=np.append(OutArr, self.Zarr)
        OutArr=OutArr.reshape(3, self.Nx*self.Ny)
        OutArr=OutArr.T
        if fmt=='npy':
            np.save(fname, OutArr)
        elif fmt=='txt':
            np.savetxt(fname, OutArr)
        else:
            raise TypeError('Wrong output format!')
        return
        
    def interpTest(self):
        points = np.array([self.XarrIn,self.YarrIn]).T
        # interpOP = sinter.CloughTocher2DInterpolator(points, self.ZarrIn, tol=1e-7, fill_value = 0.0)
        # interpOP = sinter.NearestNDInterpolator(points, self.ZarrIn)
        interpOP = sinter.LinearNDInterpolator(points, self.ZarrIn, fill_value=0.0)
        
        self.Zarr=interpOP(self.Xarr, self.Yarr)
        return 
        
    
    def interp(self, kind='cubic', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr-self.XarrIn)).sum() < 0.01 and (np.abs(self.Yarr-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn
                return
        Finterp=sinter.interp2d(self.XarrIn, self.YarrIn, self.ZarrIn,
            kind=kind, copy=copy, bounds_error=bounds_error, fill_value=fill_value)
        self.Zarr = Finterp (self.x, self.y)
        self.Zarr=self.Zarr.reshape(self.Nx*self.Ny)
        return
    
    def natgridInterp(self, interp='nn', copy=True, bounds_error=False, fill_value=np.nan):
        if self.Xarr.size == self.XarrIn.size:
            if (np.abs(self.Xarr.reshape(self.Nx*self.Ny)-self.XarrIn)).sum() < 0.01\
                and (np.abs(self.Yarr.reshape(self.Nx*self.Ny)-self.YarrIn)).sum() < 0.01:
                print 'No need to interpolate!'
                self.Zarr=self.ZarrIn
                return
        # self.Zarr = ma.getdata(griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp))
        self.Zarr = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, self.Xarr, self.Yarr, interp=interp)
        self.Zarr = ma.getdata(self.Zarr)
        # self.Zarr=self.Zarr.reshape(self.Nx*self.Ny)
        return
    
    def PlotInField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        xi, yi = np.meshgrid(self.x, self.y)
        zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(xi, yi, zi, cmap='gist_ncar_r', shading='gouraud')
        levels=np.linspace(zi.min(), zi.max(), 40)
        plt.contour(xi, yi, zi, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]])
        plt.xlabel('km')
        plt.ylabel('km')
        plt.show()
        return
        
    def PlotField(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr, self.Yarr, self.Zarr, cmap='gist_ncar_r', shading='gouraud')
        levels=np.linspace(self.Zarr.min(), self.Zarr.max(), 100)
        plt.contour(self.Xarr, self.Yarr, self.Zarr, colors='k', levels=levels)
        plt.axis([0, self.x[-1], 0, self.y[-1]])
        plt.xlabel('km')
        plt.ylabel('km')
        plt.show()
        return
    
    def CuttingEdges(self, nx, ny, fieldtype='TravelT'):
        self.Nx=self.Nx-2*nx
        self.Ny=self.Ny-2*ny
        self.x=np.arange(self.Nx)*self.dx
        self.y=np.arange(self.Ny)*self.dy
        self.Xarr, self.Yarr = np.meshgrid(self.x, self.y)
        self.Zarr=self.Zarr[nx:-nx, ny:-ny]
        return
        
    def fftDiff(self, m, n):
        try:
            h = np.fft.fft2(ma.getdata(self.Zarr))
        except:
            h = np.fft.fft2(self.Zarr)
        hshift = np.fft.fftshift(h)
        Nx=self.Nx
        Ny=self.Ny
        if Nx % 2 ==0:
            u=np.arange(Nx) - Nx/2.
        else:
            u=np.arange(Nx) - (Nx-1)/2.
        if Ny % 2 ==0:
            v=np.arange(Ny) - Ny/2.
        else:
            v=np.arange(Ny) - (Ny-1)/2.
        U,V=np.meshgrid(u,v)
        hdiff =  ((1j*2*np.pi*U/Nx)**m)*((1j*2*np.pi*V/Ny)**n) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )/(self.dx**m)/(self.dy**n)
        out_diff = out_diff[:self.Ny, :self.Nx]
        return out_diff
    
    def fftDiff2(self, m, n):
        Nx=1<<(self.Nx-1).bit_length()
        Ny=1<<(self.Ny-1).bit_length()
        # h = np.fft.fft2(self.Zarr, s=[Nx, Ny] )
        h = np.fft.fft2(ma.getdata(self.Zarr), s=[Nx, Ny] )
        hshift = np.fft.fftshift(h)
        u = np.arange(Nx) - Nx/2.
        v = np.arange(Ny) - Ny/2.
        U,V = np.meshgrid(u,v)
        hdiff = ( (1j*2*np.pi*U/Nx)**m )*( (1j*2*np.pi*V/Ny)**n ) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )/(self.dx**m)/(self.dy**n)
        out_diff = out_diff[:self.Ny, :self.Nx]
        return out_diff
    
    
    def Gradient(self, edge_order=1, method='default', order=2):
        if method=='default':
            self.grad=np.gradient( ma.getdata(self.Zarr), self.dx, self.dy, edge_order=edge_order)
        elif method=='freq':
            diff_x=self.fftDiff(m=1, n=0)
            diff_y=self.fftDiff(m=0, n=1)
            self.grad=[]
            self.grad.append(diff_y)
            self.grad.append(diff_x)
        elif method == 'convolve':
            if order==2:
                diff_x=convolve( ma.getdata(self.Zarr), X_diff_weight_2)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_2)/self.dy
            elif order==4:
                diff_x=convolve(ma.getdata(self.Zarr), X_diff_weight_4)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_4)/self.dy
            elif order==6:
                diff_x=convolve(ma.getdata(self.Zarr), X_diff_weight_6)/self.dx
                diff_y=convolve(ma.getdata(self.Zarr), Y_diff_weight_6)/self.dy
            self.grad=[]
            self.grad.append(diff_y)
            self.grad.append(diff_x)
        return
    
    def GetApparentV(self):
        self.AppV = np.sqrt ( self.grad[0] ** 2 + self.grad[1] ** 2)
        self.AppV[ np.where(self.AppV==0) ] = -1.
        self.AppV=1./self.AppV
        return
    
    def PlotAppV(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.pcolormesh(self.Xarr, self.Yarr, self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetDoT(self, enx, eny):
        Darr=np.sqrt( (self.Xarr-enx)**2 + (self.Yarr-eny)**2)
        self.DoT = Darr/ma.getdata(self.Zarr)
        return
    
    def PlotDoT(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        # plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.AppV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.pcolormesh(self.Xarr, self.Yarr, self.DoT, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def LaplacianEqualXY(self):
        if self.dx!=self.dy:
            raise ValueError('grid spacing not equal!')
        self.lplc=scipy.ndimage.filters.laplace(ma.getdata(self.Zarr) ) / (self.dx*self.dy)
        self.lplc=self.lplc[1:-1, 1:-1]
        return
    
    def Laplacian(self, method='default', order=2):
        Zarr=ma.getdata(self.Zarr)
        if method == 'default':
            Zarr_yp=Zarr[2:, 1:-1]
            Zarr_yn=Zarr[:-2, 1:-1]
            Zarr_xp=Zarr[1:-1, 2:]
            Zarr_xn=Zarr[1:-1, :-2]
            Zarr=Zarr[1:-1, 1:-1]
            self.lplc=(Zarr_yp+Zarr_yn-2*Zarr) / (self.dy**2) + (Zarr_xp+Zarr_xn-2*Zarr) / (self.dx**2)
        elif method == 'convolve':
            if order==2:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_2)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_2)/self.dy/self.dy
            elif order==4:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_4)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_4)/self.dy/self.dy
            elif order==6:
                diff2_x=convolve( ma.getdata(self.Zarr), X_diff2_weight_6)/self.dx/self.dx
                diff2_y=convolve(ma.getdata(self.Zarr), Y_diff2_weight_6)/self.dy/self.dy
            self.lplc=diff2_x+diff2_y
            self.lplc=self.lplc[1:-1, 1:-1]
        return
    
    
    
    def PlotLaplacian(self):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplc, cmap='seismic', shading='gouraud' )
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetLplcCorrection(self, per):
        omega=2.*np.pi/per
        Zarr=ma.getdata(self.Zarr)
        self.lplcCo=self.lplc/Zarr[1:-1, 1:-1]/(omega**2)
        return
    
    def PlotLplcCo(self):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.lplcCo, cmap='seismic', shading='gouraud' , vmin=-0.01, vmax=0.01)
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    def GetCorV(self, inAmpField):
        lplcCo=inAmpField.lplcCo
        try:
            self.CorV=1./ np.sqrt(1./(self.AppV**2) - lplcCo)
        except:
            self.CorV=1./ np.sqrt(1./(self.AppV[1:-1, 1:-1]**2) - lplcCo)
        return
    
    def PlotCorV(self, vmin=2.9, vmax=3.1):
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        plt.subplots()
        # xi, yi = np.meshgrid(self.x, self.y)
        # zi = griddata(self.XarrIn, self.YarrIn, self.ZarrIn, xi, yi )
        plt.pcolormesh(self.Xarr[1:-1, 1:-1], self.Yarr[1:-1, 1:-1], self.CorV, cmap='seismic_r', shading='gouraud', vmin=vmin, vmax= vmax)
        plt.axis('equal')
        plt.colorbar()
        plt.axis([self.x[1], self.x[-2], self.y[1], self.y[-2]])
        plt.xlabel('km')
        plt.ylabel('km')
        # plt.show()
        return
    
    

