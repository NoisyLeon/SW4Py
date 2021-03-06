import numpy as np
import os, sys
import warnings
flush = sys.stdout.flush()

hdr_dtype = np.dtype([
    ('magic'       , 'int32'   ),
    ('precision'   , 'int32'   ),
    ('attenuation' , 'int32'   ),
    ('az'          , 'float64' ),
    ('lon0'        , 'float64' ),
    ('lat0'        , 'float64' ),
    ('mlen'        , 'int32'   )
    ])

block_hdr_dtype = np.dtype([
    ('hh'  ,  'float64' ),
    ('hv'  ,  'float64' ),
    ('z0'  ,  'float64' ),
    ('nc'  ,  'int32'   ),
    ('ni'  ,  'int32'   ),
    ('nj'  ,  'int32'   ),
    ('nk'  ,  'int32'   )
    ])

ak135model = np.array([
    0.00, 2.7200, 5.8000, 3.4600, 1478.30, 599.99, 0.000,
    20.00, 2.7200, 5.8000, 3.4600, 1368.02, 599.99, 0.000,
    20.00, 2.9200, 6.5000, 3.8500, 1368.02, 599.99, 0.000,
    35.00, 2.9200, 6.5000, 3.8500, 950.50, 394.62, 0.000,
    35.00, 3.3200, 8.0400, 4.4800, 972.77, 403.93, 0.000,
    77.50, 3.3450, 8.0450, 4.4900, 1008.71, 417.59, 0.000,
    120.00, 3.3710, 8.0500, 4.5000, 182.57, 76.06, 0.000,
    120.00, 3.4268, 8.0505, 4.5000, 182.57, 76.06, 0.000,
    165.00, 3.3711, 8.1750, 4.5090, 188.72, 76.55, 0.000,
    210.00, 3.3243, 8.3007, 4.5184, 200.97, 79.40, 0.000,
    210.00, 3.3243, 8.3007, 4.5184, 338.47, 133.72, 0.000,
    260.00, 3.3663, 8.4822, 4.6094, 346.37, 136.38, 0.000,
    310.00, 3.4110, 8.6650, 4.6964, 355.85, 139.38, 0.000,
    360.00, 3.4577, 8.8476, 4.7832, 366.34, 142.76, 0.000,
    410.00, 3.5068, 9.0302, 4.8702, 377.93, 146.57, 0.000,
    410.00, 3.9317, 9.3601, 5.0806, 413.66, 162.50, 0.000,
    460.00, 3.9273, 9.5280, 5.1864, 417.32, 164.87, 0.000,
    510.00, 3.9233, 9.6962, 5.2922, 419.94, 166.80, 0.000,
    560.00, 3.9218, 9.8640, 5.3989, 422.55, 168.78, 0.000,
    610.00, 3.9206, 10.0320, 5.5047, 425.51, 170.82, 0.000,
    660.00, 3.9201, 10.2000, 5.6104, 428.69, 172.93, 0.000])

ak135modelCPS = np.loadtxt('ak135_cps.mod')
ak135modelCPS.reshape(ak135modelCPS.size)

### block model
class Blockmodel(object):
    """
    An object to handle single block model for SW4
    =========================================================================================
    Parameters:
    vp, vs          - P/S wave velocity real > 0 m/s none
    rho             - density real > 0 kg/m^3 none
    vpgrad, vsgrad  - vertical gradient for vp/vs real m/s/m none
    rhograd         - vertical gradient for rho real kg/m^4 none
    Qp, Qs          - P/S wave quality factor real > 0 none none
    x1, x2          - minimum/maximum x-dim for the box shaped sub-region 
    y1, y2          - minimum/maximum y-dim for the box shaped sub-region 
    z1, z2          - minimum/maximum z-dim for the box shaped sub-region 
    absdepth        - z1 and z2 relative to topography (0), or absolute z-coordinate (1)
    =========================================================================================
    """
    def __init__(self, vp=None, vs=None, rho=None, Qp=None, Qs=None, vpgrad=None, vsgrad=None, rhograd=None,
            absdepth=None, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
        self.vp=vp
        self.vs=vs
        self.rho=rho
        self.Qp=Qp
        self.Qs=Qs
        self.vpgrad=vpgrad
        self.vsgrad=vsgrad
        self.rhograd=rhograd
        self.absdepth=absdepth
        self.x1=x1
        self.x2=x2
        self.y1=y1
        self.y2=y2
        self.z1=z1
        self.z2=z2
        return
        
class BlockLst(object):
    """
    An object to handle a list of block object for SW4
    """
    def __init__(self,blocks=None):
        self.blocks=[]
        if isinstance(blocks, Blockmodel):
            blocks = [blocks]
        if blocks:
            self.blocks.extend(blocks)

    def __add__(self, other):
        """
        Add two BlockLst with self += other.
        """
        if isinstance(other, Blockmodel):
            other = BlockLst([other])
        if not isinstance(other, BlockLst):
            raise TypeError
        blocks = self.blocks + other.blocks
        return self.__class__(blocks=blocks)

    def __len__(self):
        """
        Return the number of Traces in the BlockLst object.
        """
        return len(self.blocks)

    def __getitem__(self, index):
        """
        __getitem__ method of BlockLst objects.
        :return: Blockmodel objects
        """
        if isinstance(index, slice):
            return self.__class__(blocks=self.blocks.__getitem__(index))
        else:
            return self.blocks.__getitem__(index)

    def append(self, block):
        """
        Append a single Blockmodel object to the current BlockLst object.
        """
        if isinstance(block, Blockmodel):
            self.blocks.append(block)
        else:
            msg = 'Append only supports a single Blockmodel object as an argument.'
            raise TypeError(msg)
        return self
    
    def AddSingle(self, vs, vp=None,  rho=None, Qp=None, Qs=None, vpgrad=None, vsgrad=None, rhograd=None,
            absdepth=None, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
        """
        Add single block to the model.
        ===================================================================================================================
        Input Parameters:
        vs, vp, rho             - S/P wave velocity, density
                                    If vp/rho not specified, vp and rho will be calculated by assuming Brocher crust.
        Qp, Qs                  - P/S wave quality factor
        vpgrad, vsgrad, rhograd - S/P wave velocity, density gradient( C(z) = C + z*Cgrad, z=0 at surface, NOT at z1)
        absdepth                - z1/z2 are relative to topography(0) or absolute z-coordinate
        x1, x2, y1, y2, z1, z2  - define the block region
        ===================================================================================================================
        """
        if vp ==None:
            vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
        if rho==None:
            rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5
        vp=vp*1000.
        vs=vs*1000.
        rho=rho*1000.
        self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad, rhograd=rhograd,
            absdepth=absdepth,x1=x1, x2=x2, y1=y1, y2=y2, z1=z1, z2=z2))
        return

    def ak135(self, zmax=410., CPS=True):
        """ ak135 model
        =================================================================
        Input Parameters:
        zmax    - maximum depth
        CPS     - use layered ak135 model used in CPS(homogeneous in each layer) or not
        =================================================================
        """
        if zmax>660.:
            raise ValueError('Depth is too large for Cartesian simulation.')
        zmax=zmax*1000.
        if CPS==True:
            modelArr=ak135modelCPS.reshape((210,7))
            Np=210
        else:
            modelArr=ak135model.reshape(21, 7)
            Np=21
            Qkappa=modelArr[:,4]
            Qmu=modelArr[:,5]
            alpha=modelArr[:,2]
            beta=modelArr[:,3]
            Qp_inv=4*(beta/alpha)**2/3/Qmu+(1-4*(beta/alpha)**2/3)/Qkappa
            modelArr[:,4]=1/Qp_inv
        modelArr[:,0]=modelArr[:,0]*1000.
        modelArr[:,1]=modelArr[:,1]*1000.
        modelArr[:,2]=modelArr[:,2]*1000.
        modelArr[:,3]=modelArr[:,3]*1000.
        for i in np.arange(Np):
            if modelArr[i, 0]==modelArr[i+1, 0]:
                continue
            if zmax <= modelArr[i+1,0]:
                if i==0:
                    rho=modelArr[0,1]
                    vp=modelArr[0,2]
                    vs=modelArr[0,3]
                    Qp=modelArr[0,4]
                    Qs=modelArr[0,5]
                    self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, z1=0, z2=zmax))
                else:
                    z1=modelArr[i,0]
                    rho1=modelArr[i,1]
                    vp1=modelArr[i,2]
                    vs1=modelArr[i,3]
                    Qp=modelArr[i,4]
                    Qs=modelArr[i,5]
                    z2=modelArr[i+1,0]
                    rho2=modelArr[i+1,1]
                    vp2=modelArr[i+1,2]
                    vs2=modelArr[i+1,3]
                    vpgrad=(vp2-vp1)/(z2-z1)
                    vsgrad=(vs2-vs1)/(z2-z1)
                    rhograd=(rho2-rho1)/(z2-z1)
                    vp=vp1-z1*vpgrad
                    vs=vs1-z1*vsgrad
                    rho=rho1-z1*rhograd
                    self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad,
                            rhograd=rhograd, z1=z1, z2=zmax))
                break
            z1=modelArr[i,0]
            rho1=modelArr[i,1]
            vp1=modelArr[i,2]
            vs1=modelArr[i,3]
            Qp=modelArr[i,4]
            Qs=modelArr[i,5]
            z2=modelArr[i+1,0]
            rho2=modelArr[i+1,1]
            vp2=modelArr[i+1,2]
            vs2=modelArr[i+1,3]
            vpgrad=(vp2-vp1)/(z2-z1)
            vsgrad=(vs2-vs1)/(z2-z1)
            rhograd=(rho2-rho1)/(z2-z1)
            vp=vp1-z1*vpgrad
            vs=vs1-z1*vsgrad
            rho=rho1-z1*rhograd
            self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad,
                            rhograd=rhograd, z1=z1, z2=z2))
        return
        
    def Write2Input(self, infname):
        """
        Write block model to input parameter file for sw4
        """
        InputStr = "\n\n#------------------- block models: -------------------\n"
        blockstring = 'block vp=%f vs=%f rho=%f '
        for block in self.blocks:
            tailstr=''
            if block.Qp!=None:
                if block.Qp<0:
                    raise ValueError('Negative Qp!')
                tailstr+=' qp=%f ' %(block.Qp)
            if block.Qs!=None:
                if block.Qs<0:
                    raise ValueError('Negative Qs!')
                tailstr+=' qs=%f ' %(block.Qs)
            if block.vpgrad!=None:
                tailstr+=' vpgrad=%f ' %(block.vpgrad)
            if block.vsgrad!=None:
                tailstr+=' vsgrad=%f ' %(block.vsgrad)
            if block.rhograd!=None:
                tailstr+=' rhograd=%f ' %(block.rhograd)
            if block.absdepth!=None:
                tailstr+=' absdepth=%d ' %(block.absdepth)
            if block.x1!=None:
                tailstr+=' x1=%f ' %(block.x1)
            if block.x2!=None:
                tailstr+=' x2=%f ' %(block.x2)
            if block.y1!=None:
                tailstr+=' y1=%f ' %(block.y1)
            if block.y2!=None:
                tailstr+=' y2=%f ' %(block.y2)
            if block.z1!=None:
                tailstr+=' z1=%f ' %(block.z1)
            if block.z2!=None:
                tailstr+=' z2=%f ' %(block.z2)
            
            headstr=( blockstring %(block.vp, block.vs, block.rho) )
            InputStr += headstr
            InputStr += tailstr
            InputStr += '\n'
        with open(infname, "a") as f:
            f.write(InputStr)
        f.close()
        return
    
### rfile model
class rBlock(object):
    """
    An object to handle single rfile block model for SW4
    ===================================================================================================
    Parameters:
    number                          - number(index) of rfile block
    xyextent, xzextent, yzextent    - extent in xy, xz, yz
    data                            - data array
    
    ---------------------------------------- header variables ------------------------------------------
    hh, hv                          - horizontal/vertical grid spacing
    z0                              - upper depth 
    nc                              - number of component
    ni, nj, nk                      - number of grid point in x, y, z direction
    ===================================================================================================
    Modified from pysw4 by Shahar Shani-Kadmiel
    """
    def __init__(self, number=0, hh=None, hv=None, z0=0.0, nc=3, ni=None, nj=None, nk=None,
            xyextent=(), xzextent=(), yzextent = (), data=np.array([])):
        self.number = number
        self.hh = hh
        if hh !=None and hv==None:
            self.hv = hh
        else:
            self.hv = hv
        self.z0  = z0
        self.nc  = nc
        self.ni  = ni
        self.nj  = nj
        self.nk = nk
        self.xyextent = xyextent
        self.xzextent = xzextent
        self.yzextent = yzextent
        self.data = data
        return

    def __str__(self):
        string =   'Block information:\n'
        string +=  '-----------------\n'
        string +=  '          Number : %s\n' %self.number
        string +=  ' Horizontal h, m : %s\n' %self.hh
        string +=  '   Vertical h, m : %s\n' %self.hv
        string +=  '          z 0, m : %s\n' %self.z0
        string +=  '              ni : %s\n' %self.ni
        string +=  '              nj : %s\n' %self.nj
        string +=  '              nk : %s\n' %self.nk
        string +=  '              nc : %s\n' %self.nc
        return string

    def vp(self):
        """ Return values of vp for block"""
        return self.data[:,:,:,1]

    def vs(self):
        """ Return values of vs for block"""
        return self.data[:,:,:,2]

    def rho(self):
        """ Return values of density for block"""
        return self.data[:,:,:,0]

    def qp(self):
        """ Return values of Qp for block"""
        return self.data[:,:,:,3]

    def qs(self):
        """ Return values of Qs for block"""
        return self.data[:,:,:,4]
    
    def read_block_hdr(self, f):
        """Read rfile block header"""
        hh,hv,z0,nc,ni,nj,nk = np.fromfile(f, block_hdr_dtype, 1)[0]
        self.hh = hh
        self.hv = hv
        self.z0  = z0
        self.nc  = nc
        self.ni  = ni
        self.nj  = nj
        self.nk = nk
        return hh,hv,z0,nc,ni,nj,nk
    
    def write_block_hdr(self, f):
        """Write rfile block header"""
        hh = np.float64(self.hh)
        hv = np.float64(self.hv)
        z0 = np.float64(self.z0)
        nc = np.int32(self.nc)
        ni = np.int32(self.ni)
        nj = np.int32(self.nj)
        nk = np.int32(self.nk)
    
        block_hdr = [hh,hv,z0,nc,ni,nj,nk]
        for val in block_hdr:
            f.write(val)
        return
        
        
class vprofile(object):
    def __init__(self, vsArr=np.array([]), vpArr=np.array([]), rhoArr=np.array([]), RmaxArr=np.array([]), RminArr=np.array([]),
                 z0Arr=np.array([]), HArr=np.array([]), xArr=np.array([]), yArr=np.array([]), dtypeArr=np.array([]) ):
        """
        An object to handle vertical profile, will be used by CPSPy
        =============================================================================
        output txt format:
        vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y   dtype
        
        dtype:
        0. absolute value for vs/vp/rho
        1. percentage value for dvs/dvp/drho 
        =============================================================================
        """
        self.vsArr=vsArr
        self.vpArr=vpArr
        self.rhoArr=rhoArr
        self.RmaxArr=RmaxArr
        self.RminArr=RminArr
        self.z0Arr=z0Arr
        self.HArr=HArr
        self.xArr=xArr
        self.yArr=yArr
        self.dtypeArr=dtypeArr
        return
    
    def read(self, infname):
        inArr=np.loadtxt(infname)
        self.vsArr=inArr[:,0]
        self.vpArr=inArr[:,1]
        self.rhoArr=inArr[:,2]
        self.RmaxArr=inArr[:,3]
        self.RminArr=inArr[:,4]
        self.z0Arr=inArr[:,5]
        self.HArr=inArr[:,6]
        self.xArr=inArr[:,7]
        self.yArr=inArr[:,8]
        self.dtypeArr=inArr[:,9]
        return
    
    def write(self, outfname):
        N=self.vsArr.size
        outArr=np.append(self.vsArr, self.vpArr)
        outArr=np.append(outArr, self.rhoArr)
        outArr=np.append(outArr, self.RmaxArr)
        outArr=np.append(outArr, self.RminArr)
        outArr=np.append(outArr, self.z0Arr)
        outArr=np.append(outArr, self.HArr)
        outArr=np.append(outArr, self.xArr)
        outArr=np.append(outArr, self.yArr)
        outArr=np.append(outArr, self.dtypeArr)
        outArr=outArr.reshape(10, N)
        outArr=outArr.T
        self.check()
        np.savetxt(outfname, outArr, fmt='%g', header='vs/dvs vp/dvp rho/drho Rmax Rmin z0 H x y dtype(0: m, 1: dm)')
        return
    
    def check(self):
        N = self.vsArr.size
        for i in xrange(N):
            x=self.xArr[i]
            y=self.yArr[i]
            Rmax=self.RmaxArr[i]
            H=self.HArr[i]
            index = (self.xArr==x)* (self.yArr==y) * (self.RmaxArr==Rmax) * (self.HArr==H)
            if (self.vsArr[index]).size >1:
                warnings.warn('Profile of same geometry and location exists.', UserWarning, stacklevel=1)
                return
        return
    
    def add(self, Rmax, Rmin, z0, H, x, y, dtype, vs=None, vp=None, rho=None):
        """
        Two kinds of input:
        1. Rmax or x is an array, the function will assume N vertical profile input(N=Rmax.size or x.size)
        2. If not, the function will first check if there is any previous profile with same x, y, Rmax. If yes, change
            vs, vp or rho in the existing profile; if not, add a new profile
        """
        if isinstance(x, np.ndarray) or isinstance(Rmax, np.ndarray):
            if vs ==None or vp == None or rho == None:
                raise ValueError('For array input, all model parameters must be assigned !')
            try:
                N=x.size
            except:
                N=Rmax.size
            ##################################################
            # Making all the input variables to be numpy array of size N
            ##################################################
            if not isinstance(Rmax, np.ndarray):
                Rmax=np.ones(N)*Rmax
            if not isinstance(Rmin, np.ndarray):
                Rmin=np.ones(N)*Rmin
            if not isinstance(z0, np.ndarray):
                z0=np.ones(N)*z0
            if not isinstance(H, np.ndarray):
                H=np.ones(N)*H
            if not isinstance(x, np.ndarray):
                x=np.ones(N)*x
            if not isinstance(y, np.ndarray):
                y=np.ones(N)*y
            if not isinstance(dtype, np.ndarray):
                dtype=np.ones(N)*dtype
            if not isinstance(vs, np.ndarray):
                vs=np.ones(N)*vs
            if not isinstance(vp, np.ndarray):
                vp=np.ones(N)*vp
            if not isinstance(rho, np.ndarray):
                rho=np.ones(N)*rho
        else:
            index = (self.xArr==x)* (self.yArr==y) * (self.RmaxArr==Rmax) * (self.HArr==H)
            # change model parameters if profile with same geometry/location exist
            if np.any(index):
                if vs !=None:
                    self.vsArr[index]=vs
                if vp !=None:
                    self.vpArr[index]=vp
                if rho !=None:
                    self.rhoArr[index]=rho
                return
            if dtype==0 and (vs ==None or vp == None or rho == None):
                raise ValueError('For absolute value profile, all model parameters must be assigned !')
            if vs==None:
                vs=0.
            if vp == None:
                vp=0.
            if rho == None:
                rho=0.
        self.vsArr=np.append( self.vsArr, vs)
        self.vpArr=np.append( self.vpArr, vp)
        self.rhoArr=np.append( self.rhoArr, rho)
        self.xArr=np.append( self.xArr, x)
        self.yArr=np.append( self.yArr, y)
        self.RmaxArr=np.append( self.RmaxArr, Rmax)
        self.RminArr=np.append( self.RminArr, Rmin)
        self.HArr=np.append( self.HArr, H)
        self.z0Arr=np.append( self.z0Arr, z0)
        self.dtypeArr=np.append( self.dtypeArr, dtype)
        return
    
class rModel(object):
    """
    An object to rfile header and blocks for SW4
    ===================================================================================================
    Parameters:
    rblocks                 - list of rfile block (rBlock objects)

    ---------------------------------------- header variables ------------------------------------------
    magic                   - horizontal/vertical grid spacing
    precision               - upper depth 
    attenuation             - number of component
    az, lon0, lat0          - azimuth / origin of model coordinate
    proj_str                - projection string
    nb                      - number of blocks
    ===================================================================================================
    """
    def __init__(self, filename=None, magic=1, precision = 4, attenuation = 0, az = 0., lon0= -118, lat0=37,
                 proj_str  = '+proj=utm +ellps=WGS84 +no_defs', nb = 0, rblocks=None):
        self.filename  = filename
        self.magic = magic
        self.precision  = precision
        self.attenuation = attenuation
        self.az  = az
        self.lon0  = lon0
        self.lat0 = lat0
        self.proj_str = proj_str
        self.nb = nb
        self.rblocks = []
        self.Vprofile=vprofile()
        if isinstance(rblocks, rBlock):
            rblocks = [rblocks]
        if rblocks:
            self.rblocks.extend(rblocks)
        if filename!=None:
            self.read(filename)
        return

    def __add__(self, other):
        """
        Add two rModel with self += other.
        """
        if isinstance(other, rBlock):
            other = rModel([other])
        if not isinstance(other, rModel):
            raise TypeError
        rblocks = self.rblocks + other.rblocks
        return self.__class__(rblocks=rblocks)

    def __len__(self):
        """
        Return the number of Traces in the rModel object.
        """
        return len(self.rblocks)

    def __getitem__(self, index):
        """
        __getitem__ method of rModel objects.
        :return: rBlock objects
        """
        if isinstance(index, slice):
            return self.__class__(rblocks=self.rblocks.__getitem__(index))
        else:
            return self.rblocks.__getitem__(index)

    def append(self, rblock):
        """
        Append a single rBlock object to the current rModel object.
        """
        if isinstance(rblock, rBlock):
            self.rblocks.append(rblock)
        else:
            msg = 'Append only supports a single rBlock object as an argument.'
            raise TypeError(msg)
        return self

    def __str__(self):
        string = '\nModel information:\n'
        string +=  '-----------------\n'
        string +=  '        Filename : %s\n' %self.filename
        string +=  '           lon 0 : %s\n' %self.lon0
        string +=  '           lat 0 : %s\n' %self.lat0
        string +=  '         Azimuth : %s\n' %self.az
        string +=  '    Proj4 string : %s\n' %self.proj_str
        string +=  'Number of blocks : %s\n' %self.nb

        return string
    
    def AddTopoBlock(self, ni, nj, hh=None, hv=None, data=np.array([])):
        """
        Add topography block as the first block
        ====================================================
        Input Parameters:
        ni, nj  - number of x/y points
        hh, hv  - horizontal/vertical grid spacing
        data    - input topogrphy data array
        ====================================================
        """
        nk=1
        nc=1
        z0=0.0
        if self.nb != 0:
            raise ValueError('Error: Topography block has to be the first!')
        if hh==None and hv == None :
            raise ValueError('Error: Gird spacing not specified!')
        if hh !=None and hv == None:
            hv = hh
        if hh == None and hv !=None:
            hh = hv
        xyextent = (0, (nj-1)*hh*1e-3, 0, (ni-1)*hh*1e-3)
        xzextent = (0, (ni-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)
        yzextent = (0, (nj-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)
        Ndp=ni*nj*nk*nc
        if data.size !=0 and data.size !=Ndp:
            data=np.array([])
        if data.size == 0:
            data=np.zeros(Ndp)
        if self.precision == 4:
            data=np.float32( data.reshape(ni,nj) )
        elif self.precision == 8:
            data=np.float64( data.reshape(ni,nj) )
        self.nb+=1
        self.rblocks.append( rBlock(number=self.nb, hh=hh, hv=hv, z0=z0, nc=nc, ni=ni, nj=nj, nk=nk,
                xyextent=xyextent, xzextent=xzextent, yzextent = yzextent, data=data) )
        return
        
    def AddSingleBlock(self, ni, nj, nk, hh=None, hv=None, z0=None, data=np.array([]),
        vs=None, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None):
        """
        Add single rfile block to the model
        ===========================================================================================================
        Input Parameters:
        ni, nj, nk      - number of x/y/z points
        hh, hv          - horizontal/vertical grid spacing
        z0              - upper z coordinate
        data            - input velocity data array
        vs, vp, rho     - S/P wave velocity, density
                             If vp/rho not specified, vp and rho will be calculated by assuming Brocher crust.
        Qp, Qs          - P/S wave quality factor
        *grad           - gradient ( C(z) = C + (z-z0)*Cgrad )
        ===========================================================================================================
        """
        if self.attenuation ==1:
            nc=5
        else:
            nc=3
        if data.size==0 and vs==None :
            raise ValueError('Error: Velocity model not specified!')
        if self.nb == 0 and hh==None and hv == None :
            raise ValueError('Error: Gird spacing not specified!')
        if self.nb != 0 and hh==None and hv == None :
            hh=self.rblocks[-1].hh
            hv=self.rblocks[-1].hv
        if hh !=None and hv == None:
            hv = hh
        if hh == None and hv !=None:
            hh = hv
        if self.nb == 0:
            self.AddTopoBlock(ni=ni, nj=nj, hh=hh)
            print 'No topography block, added automatically!'
            # raise ValueError('Error: First block must be topography!')    
        if self.nb==1 and z0==None:
            z0=0
        if self.nb >1 and z0==None:
            z0=self.rblocks[-1].z0 + self.rblocks[-1].nk*self.rblocks[-1].hv
        xyextent = (0, (nj-1)*hh*1e-3, 0, (ni-1)*hh*1e-3)
        xzextent = (0, (ni-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)
        yzextent = (0, (nj-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)

        if nc==5 and (Qs==None or Qp==None):
            raise ValueError('Error: Q factors not specified!')
        Ndp=ni*nj*nk*nc
        if data.size !=0 and data.size !=Ndp:
            if vs==None:
                raise ValueError('Error: Incompatible input data array and vs not specified!')
            warnings.warn('Incompatible input data array, using input vs/vp/rho.', UserWarning, stacklevel=1)
            data=np.array([])
        
        if data.size == 0:
            if vs < 100.: warnings.warn('Double check vs, the unit is m/s!', UserWarning, stacklevel=1)
            else: vs=vs/1000.
            if vp ==None: vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4 #Brocher Crust
            if vp < 100.: warnings.warn('Double check vp, the unit is m/s!', UserWarning, stacklevel=1)
            else: vp=vp/1000.
            if rho==None: rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5 #Brocher Crust
            if vs < 100.: vs=vs*1000.
            if vp < 100.: vp=vp*1000.
            if rho < 100.: rho=rho*1000.
            if vsgrad==None: vsArr = np.ones(nk)*vs
            else: vsArr = np.ones(nk)*vs + np.arange(nk)*hv*vsgrad
            if vpgrad==None: vpArr = np.ones(nk)*vp
            else: vpArr = np.ones(nk)*vp + np.arange(nk)*hv*vpgrad
            if rhograd==None: rhoArr = np.ones(nk)*rho
            else: rhoArr = np.ones(nk)*rho + np.arange(nk)*hv*rhograd
            
            Vprofile = np.append(rhoArr, vpArr)
            Vprofile = np.append(Vprofile, vsArr)
            if nc == 5:
                if Qsgrad==None:
                    QsArr=np.ones(nk)*Qs
                else:
                    QsArr = np.ones(nk)*Qs + np.arange(nk)*hv*Qsgrad
                if Qpgrad==None:
                    QpArr=np.ones(nk)*Qp
                else:
                    QpArr = np.ones(nk)*Qp + np.arange(nk)*hv*Qpgrad
                Vprofile = np.append(Vprofile, QpArr)
                Vprofile = np.append(Vprofile, QsArr)
            Vprofile = Vprofile.reshape(nc, nk)
            Vprofile = Vprofile.T
            Vprofile = Vprofile.reshape(nc*nk)
            data = np.tile(Vprofile, ni*nj)
            if self.precision == 4:
                data=np.float32( data.reshape(ni,nj,nk,nc) ) 
            elif self.precision == 8:
                data=np.float64( data.reshape(ni,nj,nk,nc) ) 
        self.nb+=1
        self.rblocks.append( rBlock(number=self.nb, hh=hh, hv=hv, z0=z0, nc=nc, ni=ni, nj=nj, nk=nk,
                xyextent=xyextent, xzextent=xzextent, yzextent = yzextent, data=data) )
        return
    
    def ak135(self, ni, nj, zmin=0., zmax=410.,  hh=None, hv=None, CPS=True):
        """
        Implement ak135 model
        ============================================================
        Input Parameters:
        ni, nj      - number of x/y points
        zmin, zmax  - manimum/maximum depth
        hh, hv      - horizontal/vertical grid spacing
        CPS         - whether to use layered ak135 from CPS or not
        ============================================================
        """
        if zmax>660.:
            raise ValueError('Depth is too large for Cartesian simulation.')
        zmax=zmax*1000.
        if self.nb == 0 and hh==None and hv == None :
            raise ValueError('Error: Gird spacing not specified!')
        if self.nb != 0 and hh==None and hv == None :
            hh=self.rblocks[-1].hh
            hv=self.rblocks[-1].hv
        if hh !=None and hv == None:
            hv = hh
        if hh == None and hv !=None:
            hh = hv
        if self.nb == 0:
            self.AddTopoBlock(ni=ni, nj=nj, hh=hh)
            print 'No topography block, added automatically!'
        if CPS==True:
            modelArr=ak135modelCPS.reshape((210,7))
        else:
            modelArr=ak135model.reshape(21, 7)
            Qkappa=modelArr[:,4]
            Qmu=modelArr[:,5]
            alpha=modelArr[:,2]
            beta=modelArr[:,3]
            Qp_inv=4*(beta/alpha)**2/3/Qmu+(1-4*(beta/alpha)**2/3)/Qkappa
            modelArr[:,4]=1/Qp_inv
        modelArr[:,0]=modelArr[:,0]*1000.
        modelArr[:,1]=modelArr[:,1]*1000.
        modelArr[:,2]=modelArr[:,2]*1000.
        modelArr[:,3]=modelArr[:,3]*1000.
        for i in np.arange(21):
            if modelArr[i,0]==modelArr[i+1,0]: continue
            if modelArr[i+1,0]<zmin: continue
            if zmax <= modelArr[i+1,0]:
                if i==0:
                    rho=modelArr[0,1]
                    vp=modelArr[0,2]
                    vs=modelArr[0,3]
                    Qp=modelArr[0,4]+(modelArr[1,4]-modelArr[0,4])/modelArr[1,0] * zmin
                    Qs=modelArr[0,5]
                    nk=int( (zmax-zmin) /hv) + 1
                    Qpgrad = (modelArr[1,4]-modelArr[0,4])/modelArr[1,0]
                    self.AddSingleBlock(ni=ni, nj=nj, nk=nk, hh=hh, hv=hv, z0=zmin, vs=vs, vp=vp, rho=rho, Qs=Qs, Qp=Qp, Qpgrad=Qpgrad)
                else:
                    z1=modelArr[i,0]
                    rho1=modelArr[i,1]
                    vp1=modelArr[i,2]
                    vs1=modelArr[i,3]
                    Qp1=modelArr[i,4]
                    Qs1=modelArr[i,5]
                    z2=modelArr[i+1,0]
                    rho2=modelArr[i+1,1]
                    vp2=modelArr[i+1,2]
                    vs2=modelArr[i+1,3]
                    Qp2=modelArr[i+1,4]
                    Qs2=modelArr[i+1,5]
                    vpgrad=(vp2-vp1)/(z2-z1)
                    vsgrad=(vs2-vs1)/(z2-z1)
                    rhograd=(rho2-rho1)/(z2-z1)
                    Qpgrad=(Qp2-Qp1)/(z2-z1)
                    Qsgrad=(Qs2-Qs1)/(z2-z1)
                    if zmin <= z1:
                        nk=int( (zmax-z1) /hv) + 1
                        self.AddSingleBlock(ni=ni, nj=nj, nk=nk, hh=hh, hv=hv, z0=z1, vs=vs1, vp=vp1, rho=rho1, Qs=Qs1, Qp=Qp1,
                                vsgrad=vsgrad, vpgrad=vpgrad, rhograd=rhograd, Qsgrad=Qsgrad, Qpgrad=Qpgrad)
                    else:
                        vp = vp1 + vpgrad * (zmin-z1)
                        vs = vs1 + vsgrad * (zmin-z1)
                        rho = rho1 + rhograd * (zmin-z1)
                        Qp = Qp1 + Qpgrad * (zmin-z1)
                        Qs = Qs1 + Qsgrad * (zmin-z1)
                        nk=int( (zmax-zmin) /hv) + 1
                        self.AddSingleBlock(ni=ni, nj=nj, nk=nk, hh=hh, hv=hv, z0=zmin, vs=vs, vp=vp, rho=rho, Qs=Qs, Qp=Qp,
                                vsgrad=vsgrad, vpgrad=vpgrad, rhograd=rhograd, Qsgrad=Qsgrad, Qpgrad=Qpgrad)
                break
            z1=modelArr[i,0]
            rho1=modelArr[i,1]
            vp1=modelArr[i,2]
            vs1=modelArr[i,3]
            Qp1=modelArr[i,4]
            Qs1=modelArr[i,5]
            z2=modelArr[i+1,0]
            rho2=modelArr[i+1,1]
            vp2=modelArr[i+1,2]
            vs2=modelArr[i+1,3]
            Qp2=modelArr[i+1,4]
            Qs2=modelArr[i+1,5]
            vpgrad=(vp2-vp1)/(z2-z1)
            vsgrad=(vs2-vs1)/(z2-z1)
            rhograd=(rho2-rho1)/(z2-z1)
            Qpgrad=(Qp2-Qp1)/(z2-z1)
            Qsgrad=(Qs2-Qs1)/(z2-z1)
            
            if z1 < zmin:
                vp = vp1 + vpgrad * (zmin-z1)
                vs = vs1 + vsgrad * (zmin-z1)
                rho = rho1 + rhograd * (zmin-z1)
                Qp = Qp1 + Qpgrad * (zmin-z1)
                Qs = Qs1 + Qsgrad * (zmin-z1)
                nk=int( (z2-zmin) /hv) + 1
                self.AddSingleBlock(ni=ni, nj=nj, nk=nk, hh=hh, hv=hv, z0=zmin, vs=vs, vp=vp, rho=rho, Qs=Qs, Qp=Qp,
                                    vsgrad=vsgrad, vpgrad=vpgrad, rhograd=rhograd, Qsgrad=Qsgrad, Qpgrad=Qpgrad)
            else:
                nk=int( (z2-z1) /hv) + 1
                self.AddSingleBlock(ni=ni, nj=nj, nk=nk, hh=hh, hv=hv, z0=z1, vs=vs1, vp=vp1, rho=rho1, Qs=Qs1, Qp=Qp1,
                                    vsgrad=vsgrad, vpgrad=vpgrad, rhograd=rhograd, Qsgrad=Qsgrad, Qpgrad=Qpgrad)
        return

    def BlockAnomaly(self, xmin, xmax, ymin, ymax, dm, mname='vs', zmin=0, zmax=None, nb=None):
        """
        Implement block anomaly
        ============================================================================
        Input Parameters:
        xmin, xmax, ymin, ymax, zmin, zmax  - defines the bound (in meter)
        dm                                  - model parameter anomaly in percentage
        mname                               - model variable name
        nb                                  - block number(index)
        ============================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        mtype=dictparam[mname]
        if nb!=1:
            print 'Adding block anomaly to',mname,'!'
        else:
            print 'Adding block anomaly to topography !'
        if nb==None:
            for b in np.arange(self.nb-1)+1:
                xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
                yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
                zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv + self.rblocks[b].z0
                # meshgrid notes: 
                # In the 3-D case with inputs of length M, N and P,
                # outputs are of shape (N, M, P) for 'xy' indexing and (M, N, P) for 'ij' indexing.
                xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
                tempdata=self.rblocks[b].data[:, :, :, mtype]
                xindex=(xgrid>=xmin)*(xgrid<=xmax)
                yindex=(ygrid>=ymin)*(ygrid<=ymax)
                if zmax==None:
                    zindex=zgrid>=zmin
                else:
                    zindex=(zgrid>=zmin)*(zgrid<=zmax)
                self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata*dm*xindex*yindex*zindex
        else:
            b=nb-1
            xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
            yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
            zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv
            xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
            tempdata=self.rblocks[b].data[:, :, :, mtype]
            xindex=(xgrid>=xmin)*(xgrid<=xmax)
            yindex=(ygrid>=ymin)*(ygrid<=ymax)
            if b!=0:
                if zmax==None:
                   zindex=zgrid>=zmin
                else:
                   zindex=(zgrid>=zmin)*(zgrid<=zmax)
            else:
                zindex=1.;
            self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata*dm*xindex*yindex*zindex
        return
    
    def CylinderHomoAnomaly(self, x0, y0, R, dm, va=None, mname='vs', zmin=0, zmax=None, nb=None):
        """
        Inplement homogeneous cylinder anomaly
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        dm          - model parameter anomaly in percentage
        va          - absolute value of velocity of the anomaly
        mname       - model variable name
        zmin, zmax  - defines the upper/lower bound (in meter)
        nb          - block number(index)
        ------------------------------------------------------------------------
        Version History:
                    Nov 28th, 2016  - add va to the function
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        mtype=dictparam[mname]
        if nb!=1: print 'Adding homogeneous cynlinder anomaly to',mname,'!'
        else: print 'Adding homogeneous cynlinder to topography !'
        if nb==None:
            for b in np.arange(self.nb-1)+1:
                xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
                yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
                zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv + self.rblocks[b].z0
                # meshgrid notes: 
                # In the 3-D case with inputs of length M, N and P,
                # outputs are of shape (N, M, P) for 'xy' indexing and (M, N, P) for 'ij' indexing.
                xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
                tempdata=self.rblocks[b].data[:, :, :, mtype]
                dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
                Rindex = dArr < R
                if zmax==None: zindex=zgrid>=zmin
                else: zindex=(zgrid>=zmin)*(zgrid<zmax)
                # # # self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata*dm*Rindex*zindex OLD
                if va ==None: tempdata[Rindex*zindex]=tempdata[Rindex*zindex]*(1+dm)
                else: tempdata[Rindex*zindex]=va
                self.rblocks[b].data[:, :, :, mtype] = tempdata
        else:
            b=nb-1
            xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
            yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
            zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv
            xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
            tempdata=self.rblocks[b].data[:, :, :, mtype]
            dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
            Rindex = dArr <= R;
            if b!=0:
                if zmax==None: zindex=zgrid>=zmin
                else: zindex=(zgrid>=zmin)*(zgrid<zmax)
            else: zindex=1.;
            # # # self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata*dm*Rindex*zindex OLD
            if va ==None: tempdata[Rindex*zindex]=tempdata[Rindex*zindex]*(1+dm)
            else: tempdata[Rindex*zindex]=va
            self.rblocks[b].data[:, :, :, mtype] = tempdata
        ### Add anomaly to vprofile object
        if zmax==None: zmax=9999.
        if mname=='vs':
            self.Vprofile.add(Rmax=R, Rmin=0., z0=zmin, H=zmax-zmin, x=x0, y=y0, dtype=1, vs=dm)
        elif mname=='vp':
            self.Vprofile.add(Rmax=R, Rmin=0., z0=zmin, H=zmax-zmin, x=x0, y=y0, dtype=1, vp=dm)
        elif mname=='rho':
            self.Vprofile.add(Rmax=R, Rmin=0., z0=zmin, H=zmax-zmin, x=x0, y=y0, dtype=1, rho=dm)
        return
        
    def CylinderLinearAnomaly(self, x0, y0, R, dm, mname='vs', zmin=0, zmax=None, nb=None):
        """
        Inplement linear varying cylinder anomaly
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        dm          - model parameter anomaly in percentage
        mname       - model variable name
        zmin, zmax  - defines the upper/lower bound (in meter)
        nb          - block number(index)
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        mtype=dictparam[mname]
        if nb!=1:
            print 'Adding linear cynlinder anomaly to',mname,'!'
        else:
            print 'Adding linear cynlinder to topography !'
        if nb==None:
            for b in np.arange(self.nb-1)+1:
                xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
                yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
                zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv + self.rblocks[b].z0
                # meshgrid notes: 
                # In the 3-D case with inputs of length M, N and P,
                # outputs are of shape (N, M, P) for 'xy' indexing and (M, N, P) for 'ij' indexing.
                xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
                tempdata=self.rblocks[b].data[:, :, :, mtype]
                delD = R - np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
                Rindex = (delD>=0)
                if zmax==None:
                    zindex=zgrid>=zmin
                else:
                    zindex=(zgrid>=zmin)*(zgrid<=zmax)
                self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata * dm * zindex * Rindex * delD/R 
        else:
            b=nb-1
            xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
            yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
            zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv
            xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
            tempdata=self.rblocks[b].data[:, :, :, mtype]
            delD = R - np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
            Rindex = (delD>=0)
            if b!=0:
                if zmax==None:
                   zindex=zgrid>=zmin
                else:
                   zindex=(zgrid>=zmin)*(zgrid<=zmax)
            else:
                zindex=1.;
            self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata * dm * zindex * Rindex * delD/R 
        return
    
    def CylinderCosineAnomaly(self, x0, y0, R, dm, mname='vs', zmin=0, zmax=None, nb=None):
        """
        Inplement cosine varying cylinder anomaly.
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        dm          - model parameter anomaly in percentage
        mname       - model variable name
        zmin, zmax  - defines the upper/lower bound (in meter)
        nb          - block number(index)
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        mtype=dictparam[mname]
        if nb>1:
            print 'Adding cosine cynlinder anomaly to',mname,'!'
        else:
            print 'Adding cosine cynlinder to topography !'
        if nb==None:
            for b in np.arange(self.nb-1)+1:
                xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
                yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
                zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv + self.rblocks[b].z0
                # meshgrid notes: 
                # In the 3-D case with inputs of length M, N and P,
                # outputs are of shape (N, M, P) for 'xy' indexing and (M, N, P) for 'ij' indexing.
                xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
                tempdata=self.rblocks[b].data[:, :, :, mtype]
                dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
                Rindex = (dArr <= R)
                if zmax==None:
                    zindex=zgrid>=zmin
                else:
                    zindex=(zgrid>=zmin)*(zgrid<=zmax)
                IndexIn=zindex*Rindex
                self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata * dm * IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. 
        else:
            b=nb-1
            xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
            yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
            zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv
            xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
            tempdata=self.rblocks[b].data[:, :, :, mtype]
            dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
            Rindex = (dArr <= R)
            if b!=0:
                if zmax==None:
                   zindex=zgrid>=zmin
                else:
                   zindex=(zgrid>=zmin)*(zgrid<=zmax)
            else:
                zindex=1.;
            IndexIn=zindex*Rindex
            self.rblocks[b].data[:, :, :, mtype] = tempdata + tempdata * dm * IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. 
        return
    
    def CylinderHomoSediment(self, x0, y0, R, zmax, vs, vp=None, rho=None ):
        """
        Implement cylindrical sedimentary basin with constant depth to the first rblock.
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        zmax        - maximum depth of the basin (in meter)
        vs          - vs for the basin (in m/s)
        vp, rho     - vp/rho for the basin (default is Brocher Crust)
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        xArr=np.arange(self.rblocks[1].ni)*self.rblocks[1].hh
        yArr=np.arange(self.rblocks[1].nj)*self.rblocks[1].hh
        zArr=np.arange(self.rblocks[1].nk)*self.rblocks[1].hv
        if zmax > zArr.max():
            raise ValueError('Maximum depth of sedimentary basin is too large!')
        xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij')
        vs=vs/1000.
        if vp ==None:
            vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
        if rho==None:
            rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5
        vs=vs*1000.
        vp=vp*1000.
        rho=rho*1000.
        tempdataVs=self.rblocks[1].data[:, :, :, 2]
        tempdataVp=self.rblocks[1].data[:, :, :, 1]
        tempdataRho=self.rblocks[1].data[:, :, :, 0]
        dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
        Rindex = (dArr <= R)
        zmaxArr= Rindex * zmax
        zIndexT = (zgrid <= zmaxArr)*Rindex
        zIndexB = np.logical_not(zIndexT)
        self.rblocks[1].data[:, :, :, 0] = zIndexT * rho+ tempdataRho*zIndexB
        self.rblocks[1].data[:, :, :, 1] = zIndexT * vp+ tempdataVp*zIndexB
        self.rblocks[1].data[:, :, :, 2] = zIndexT * vs+ tempdataVs*zIndexB
        self.Vprofile.add(Rmax=R, Rmin=0, z0=0., H=zmax, x=x0, y=y0, dtype=0, vs=vs, vp=vp, rho=rho)
        return 
    
    def CylinderCosineSediment(self, x0, y0, R, zmax, vs=None, vp=None, rho=None, qs=None ):
        """
        Implement cosine varying cylindrical sedimentary basin to the first rblock.
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        zmax        - maximum depth of the basin (in meter)
        vs          - vs for the basin (in m/s)
        vp, rho     - vp/rho for the basin (default is Brocher Crust)
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        xArr=np.arange(self.rblocks[1].ni)*self.rblocks[1].hh
        yArr=np.arange(self.rblocks[1].nj)*self.rblocks[1].hh
        zArr=np.arange(self.rblocks[1].nk)*self.rblocks[1].hv
        if zmax > zArr.max():
            raise ValueError('Maximum depth of sedimentary basin is too large!')
        xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij')
        dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
        Rindex = (dArr <= R)
        zmaxArr= Rindex * ( 1+np.cos( np.pi* dArr / R ) ) /2.* zmax
        zIndexT = (zgrid <= zmaxArr)*Rindex
        zIndexB = np.logical_not(zIndexT)
        if vs!=None:
            vs=vs/1000.
            if vp ==None:
                vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
            if rho==None:
                rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5
            vs=vs*1000.
            vp=vp*1000.
            rho=rho*1000.
            tempdataVs=self.rblocks[1].data[:, :, :, 2]
            tempdataVp=self.rblocks[1].data[:, :, :, 1]
            tempdataRho=self.rblocks[1].data[:, :, :, 0]
            self.rblocks[1].data[:, :, :, 0] = zIndexT * rho+ tempdataRho*zIndexB
            self.rblocks[1].data[:, :, :, 1] = zIndexT * vp+ tempdataVp*zIndexB
            self.rblocks[1].data[:, :, :, 2] = zIndexT * vs+ tempdataVs*zIndexB
        else:
            vs=self.rblocks[1].data[0, 0, 0, 2]
            vp=self.rblocks[1].data[0, 0, 0, 1]
        if self.attenuation==1 and qs!=None:
            tempdataQs=self.rblocks[1].data[:, :, :, 4]
            tempdataQp=self.rblocks[1].data[:, :, :, 3]
            qp=3/4.*qs*(vp/vs)**2
            self.rblocks[1].data[:, :, :, 4] = zIndexT * qs+ tempdataQs*zIndexB
            self.rblocks[1].data[:, :, :, 3] = zIndexT * qp+ tempdataQp*zIndexB
        return
    
    def CynlinderRingBasin(self, x0, y0, zmax, Rmax, vs, vp=None, rho=None, nr=None, dR=None, Rmin=0, outfname=None):
        """
        Implement cosine varying cylindrical sedimentary basin to the first rblock.
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle ( in meter)
        zmax        - maximum depth of the basin (in meter)
        Rmax        - maximum radius (in meter)
        Rmin        - minimum radius (in meter)
        nr          - number of rings
        dR          - ring interval
        vs          - vs for the basin (in m/s)
        vp, rho     - vp/rho for the basin (default is Brocher Crust)
        outfname    - output txt file for CPSPy 
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        xArr=np.arange(self.rblocks[1].ni)*self.rblocks[1].hh
        yArr=np.arange(self.rblocks[1].nj)*self.rblocks[1].hh
        zArr=np.arange(self.rblocks[1].nk)*self.rblocks[1].hv
        if zmax > zArr.max():
            raise ValueError('Maximum depth of sedimentary basin is too large!')
        xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij')
        vs=vs/1000.
        if vp ==None: vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
        if rho==None: rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5
        vs=vs*1000.
        vp=vp*1000.
        rho=rho*1000.
        tempdataVs=self.rblocks[1].data[:, :, :, 2]
        tempdataVp=self.rblocks[1].data[:, :, :, 1]
        tempdataRho=self.rblocks[1].data[:, :, :, 0]
        if nr == None or dR == None:
            nr = int( (Rmax-Rmin) /self.rblocks[1].hh/10 ) -1
            dR = 10.*self.rblocks[1].hh
        dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
        zIndexT=np.zeros(tempdataVs.shape, dtype='bool')
        HArr=np.array([])
        RmaxArr=np.array([])
        RminArr=np.array([])
        ### 
        for ir in np.arange(nr):
            cRmax = Rmax - ir * dR
            cRmin = Rmax - (ir+1) * dR
            cRmean=(cRmax+cRmin)/2.
            Rindex = (dArr <= cRmax) * (dArr > cRmin)
            H =( np.int( ( 1+np.cos( np.pi* cRmean / Rmax ) ) /2.* zmax / self.rblocks[1].hv) +1.)*  self.rblocks[1].hv
            # print cRmax, cRmin, cRmean, H
            if H > zmax:
                H=zmax
            zmaxArr= Rindex * H
            zIndexT = zIndexT + (zgrid <= zmaxArr)*Rindex
            HArr=np.append(HArr, H)
            RmaxArr=np.append(RmaxArr, cRmax)
            RminArr=np.append(RminArr, cRmin)
        Rindex = (dArr <= cRmin)
        zIndexT = zIndexT + (zgrid <= zmax)*Rindex
        zIndexB = np.logical_not(zIndexT)
        HArr=np.append(HArr, zmax)
        RmaxArr=np.append(RmaxArr, cRmin)
        RminArr=np.append(RminArr, 0)
        self.rblocks[1].data[:, :, :, 0] = zIndexT * rho+ tempdataRho*zIndexB
        self.rblocks[1].data[:, :, :, 1] = zIndexT * vp+ tempdataVp*zIndexB
        self.rblocks[1].data[:, :, :, 2] = zIndexT * vs+ tempdataVs*zIndexB
        self.Vprofile.add(Rmax=RmaxArr, Rmin=RminArr, z0=0., H=HArr, x=x0, y=y0, dtype=0, vs=vs, vp=vp, rho=rho)
        if outfname != None: self.Vprofile.write(outfname=outfname)
        return
    
    def CylinderLinearDepthAnomaly(self, x0, y0, R, vt, vb, zmax, mname='vs', zmin=0, nb=None, outfname=None):
        """
        Inplement linear varying cylinder anomaly
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        dm          - model parameter anomaly in percentage
        mname       - model variable name
        zmin, zmax  - defines the upper/lower bound (in meter)
        nb          - block number(index)
        ========================================================================
        """
        dictparam={'rho':0, 'vp': 1 , 'vs': 2 , 'qp': 3 , 'qs': 4 }
        mtype=dictparam[mname]
        if nb!=1: print 'Adding homogeneous cynlinder anomaly to',mname,'!'
        else: print 'Adding homogeneous cynlinder to topography !'
        tzArr=np.array([])
        if nb==None:
            for b in np.arange(self.nb-1)+1:
                xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
                yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
                zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv + self.rblocks[b].z0
                tzArr = np.append(tzArr, zArr[(zArr>=zmin)*(zArr<=zmax)])
                # meshgrid notes: 
                # In the 3-D case with inputs of length M, N and P,
                # outputs are of shape (N, M, P) for 'xy' indexing and (M, N, P) for 'ij' indexing.
                xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
                tempdata=self.rblocks[b].data[:, :, :, mtype]
                dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
                Rindex = dArr < R
                if zmax==None: zindex=zgrid>=zmin
                else: zindex=(zgrid>=zmin)*(zgrid<=zmax)
                vArr = zgrid / (zmax-zmin) * (vb-vt) + vt
                tempdata[Rindex*zindex]=vArr[Rindex*zindex]
                self.rblocks[b].data[:, :, :, mtype] = tempdata
        else:
            b=nb-1
            xArr=np.arange(self.rblocks[b].ni)*self.rblocks[b].hh
            yArr=np.arange(self.rblocks[b].nj)*self.rblocks[b].hh
            zArr=np.arange(self.rblocks[b].nk)*self.rblocks[b].hv
            tzArr = np.append(tzArr, zArr[(zArr>=zmin)*(zArr<=zmax)])
            xgrid, ygrid, zgrid = np.meshgrid(xArr, yArr, zArr, indexing='ij') 
            tempdata=self.rblocks[b].data[:, :, :, mtype]
            dArr = np.sqrt( (xgrid-x0)**2 + (ygrid-y0)**2)
            Rindex = dArr <= R;
            if b!=0:
                if zmax==None: zindex=zgrid>=zmin
                else: zindex=(zgrid>=zmin)*(zgrid<=zmax)
            else: zindex=1.;
            vArr = zgrid / (zmax-zmin) * (vb-vt) + vt
            tempdata[Rindex*zindex]=vArr[Rindex*zindex]
            self.rblocks[b].data[:, :, :, mtype] = tempdata
        if outfname !=None:
            zArrT=tzArr[:-1]; zArrB=tzArr[1:]
            zArrm = (zArrB + zArrT)/2.
            HArr = zArrB - zArrT
            vHArr = zArrm / (zmax-zmin) * (vb-vt) + vt
            outArr = np.append(zArrT/1000., HArr/1000.)
            outArr = np.append(outArr, vHArr/1000.)
            outArr = outArr.reshape(3, HArr.size)
            outArr = outArr.T
            np.savetxt(outfname, outArr, fmt='%g', header='z0 H vs')
        else:
            self.tzArr=tzArr
        return
    
    def CylinderLinearDepthAnomalyAll(self, x0, y0, R, vt, vb, zmax, zmin=0, nb=None, outfname=None):
        """
        Inplement linear varying cylinder anomaly
        ========================================================================
        Input Parameters:
        x0, y0      - the center of the circle( in meter)
        R           - radius (in meter)
        dm          - model parameter anomaly in percentage
        mname       - model variable name
        zmin, zmax  - defines the upper/lower bound (in meter)
        nb          - block number(index)
        ========================================================================
        """
        # top model parameters
        vst=vt/1000.
        vpt=0.9409+2.0947*vst-0.8206*vst**2+0.2683*vst**3-0.0251*vst**4
        rhot=1.6612*vpt-0.4721*vpt**2+0.0671*vpt**3-0.0043*vpt**4+0.000106*vpt**5
        vst=vst*1000.
        vpt=vpt*1000.
        rhot=rhot*1000.
        # bottom model parameters
        vsb=vb/1000.
        vpb=0.9409+2.0947*vsb-0.8206*vsb**2+0.2683*vsb**3-0.0251*vsb**4
        rhob=1.6612*vpb-0.4721*vpb**2+0.0671*vpb**3-0.0043*vpb**4+0.000106*vpb**5
        vsb=vsb*1000.
        vpb=vpb*1000.
        rhob=rhob*1000.
        self.CylinderLinearDepthAnomaly(x0=x0, y0=y0, R=R, vt=vst, vb=vsb, zmax=zmax, mname='vs', zmin=zmin, nb=nb)
        self.CylinderLinearDepthAnomaly(x0=x0, y0=y0, R=R, vt=vpt, vb=vpb, zmax=zmax, mname='vp', zmin=zmin, nb=nb)
        self.CylinderLinearDepthAnomaly(x0=x0, y0=y0, R=R, vt=rhot, vb=rhob, zmax=zmax, mname='rho', zmin=zmin, nb=nb)
        if outfname !=None:
            zArrT  = self.tzArr[:-1]; zArrB=self.tzArr[1:]
            zArrm  = (zArrB + zArrT)/2.
            HArr   = zArrB - zArrT
            vsHArr = zArrm / (zmax-zmin) * (vsb-vst) + vst
            vpHArr = zArrm / (zmax-zmin) * (vpb-vpt) + vpt
            rhoHArr= zArrm / (zmax-zmin) * (rhob-rhot) + rhot
            outArr = np.append(zArrT/1000., HArr/1000.)
            outArr = np.append(outArr, vsHArr/1000.)
            outArr = np.append(outArr, vpHArr/1000.)
            outArr = np.append(outArr, rhoHArr/1000.)
            outArr = outArr.reshape(5, HArr.size)
            outArr = outArr.T
            np.savetxt(outfname, outArr, fmt='%g', header='z0 H vs vp rho')
        return
    
    def read_hdr(self, f):
        """
        Read rfile header
        """
        (magic,precision,
            attenuation,az,lon0,lat0,mlen) = np.fromfile(f, hdr_dtype, 1)[0]
        proj_str_dtype = 'S' + str(mlen)
        proj_str = np.fromfile(f, proj_str_dtype, 1)[0]
        nb = np.fromfile(f, 'int32', 1)[0]
        self.magic=magic
        self.precision  = precision
        self.attenuation    = attenuation
        self.az  = az
        self.lon0  = lon0
        self.lat0 = lat0
        self.mlen=mlen
        self.proj_str = proj_str
        self.nb = nb
        return magic,precision,attenuation,az,lon0,lat0,mlen,proj_str,nb
    
    def write_hdr(self, f):
       """
       Write rfile header
       """
       magic        = np.int32(self.magic)
       precision    = np.int32(self.precision)
       attenuation  = np.int32(self.attenuation)
       az           = np.float64(self.az)
       lon0         = np.float64(self.lon0)
       lat0         = np.float64(self.lat0)
       mlen         = np.int32(len(self.proj_str))
       nb           = np.int32(self.nb)
       proj_str = self.proj_str
   
       hdr = [magic,precision,attenuation,az,lon0,lat0,mlen,proj_str,nb]
       for val in hdr:
           f.write(val)
       return
    
    def read(self, filename, verbose=False ):
        """
        Read rfile into Model object
        """
        self.filename = filename
        with open(filename, 'rb') as f:
            self.read_hdr(f)
            if verbose:
                print self
                flush
            for b in range(self.nb):
                rblock = rBlock()
                rblock.number = b+1
                rblock.read_block_hdr(f)
                self.append(rblock)
                if verbose:
                    print block
                    flush
            for b in np.arange(self.nb):
                hh = self.rblocks[b].hh
                hv = self.rblocks[b].hv
                z0 = self.rblocks[b].z0
                nc = self.rblocks[b].nc
                ni = self.rblocks[b].ni
                nj = self.rblocks[b].nj
                nk = self.rblocks[b].nk
                z = np.fromfile(f, np.float32, ni*nj*nk*nc)
                if nc == 1: # topo is independant of k
                    self.rblocks[b].data = z.reshape(ni,nj)
                else:
                    # C-order reshape
                    self.rblocks[b].data = z.reshape(ni,nj,nk,nc)
                self.rblocks[b].xyextent = (0, (nj-1)*hh*1e-3, 0, (ni-1)*hh*1e-3)
                self.rblocks[b].xzextent = (0, (ni-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)
                self.rblocks[b].yzextent = (0, (nj-1)*hh*1e-3, (z0+(nk-1)*hv)*1e-3, z0*1e-3)
        return

    def write(self, filename, verbose=False ):
        """
        Write Model object into rfile
        """
        self.filename = filename
        with open(filename, 'wb') as f:
            self.write_hdr(f)
            if verbose:
                print self
                flush
            for b in np.arange(self.nb):
                rblock = self.rblocks[b]
                if b==0 and ( rblock.nc!=1 or rblock.nk!=1):
                    raise ValueError('Error: First block must be topography!')
                rblock.write_block_hdr(f)
                if verbose:
                    print block
                    flush
            for b in np.arange(self.nb):                
                z=self.rblocks[b].data.reshape(self.rblocks[b].data.size)
                z.tofile(f)
        return
    
    def writeVprofile(self, filename):
        self.Vprofile.write(outfname=filename)
        return
        
    def checkInput(self, infname):
        print '=========== Checking input with rfile ==========='
        with open(infname, 'rb') as f:
            for line in f.readlines():
                cline=line.split()
                if cline[0]=='grid':
                    cline = line
                    break
        h=float( (cline.split('h=')[1].split())[0] )
        nx=float( (cline.split('nx=')[1].split())[0] )
        ny=float( (cline.split('ny=')[1].split())[0] )
        nz=float( (cline.split('nz=')[1].split())[0] )
        for rblock in self.rblocks:
            if rblock.hh*rblock.ni != h *nx or rblock.hv*rblock.nj != h *ny:
                raise ValueError('Incompatible x or y size!')
            if h != rblock.hv or h != rblock.hh:
                warnings.warn('Incompatible grid spacing: hv='+str(rblock.hv)+' hh='+str(rblock.hh)+' h='+str(h), UserWarning, stacklevel=1)
            if rblock.ni != nx or rblock.nj !=ny:
                warnings.warn('Incompatible ni='+str(rblock.ni) +' nx='+str(nx)+' nj='+str(rblock.nj)+' ny='+str(ny), UserWarning, stacklevel=1)
        if rblock.hv * (rblock.nk-1) + rblock.z0 < (nz-1) * h:
            print rblock.hv * (rblock.nk-1) + rblock.z0 , (nz-1) * h
            raise ValueError('Depth size too small !')
        print '=========== Checked input with rfile  ==========='
        return 
        

    
    
    