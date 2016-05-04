import numpy as np
import os, sys

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
    660.00, 3.9201, 10.2000, 5.6104, 428.69, 172.93, 0.000]);

### block model
class Blockmodel(object):
    def __init__(self, vp=None, vs=None, rho=None, Qp=None, Qs=None, vpgrad=None, vsgrad=None, rhograd=None,
            absdepth=None, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
        self.vp=vp;
        self.vs=vs;
        self.rho=rho;
        self.Qp=Qp;
        self.Qs=Qs;
        self.vpgrad=vpgrad;
        self.vsgrad=vsgrad;
        self.rhograd=rhograd;
        self.absdepth=absdepth;
        self.x1=x1;
        self.x2=x2;
        self.y1=y1;
        self.y2=y2;
        self.z1=z1;
        self.z2=z2;
        return;
        
class BlockLst(object):
    
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
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
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
        return self;
    
    def AddSingle(self, vs, vp=None,  rho=None, Qp=None, Qs=None, vpgrad=None, vsgrad=None, rhograd=None,
            absdepth=None, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
        """
        Add single block to the model.
        If not specified, vp and rho will be calculated by assuming Brocher crust.
        """
        if vp ==None:
            vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4;
        if rho==None:
            rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5;
        vp=vp*1000.;
        vs=vs*1000.;
        rho=rho*1000.;
        self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad, rhograd=rhograd,
            absdepth=absdepth,x1=x1, x2=x2, y1=y1, y2=y2, z1=z1, z2=z2));
        return;

    def ak135(self, zmax=410.):
        if zmax>660.:
            raise ValueError('Depth is too large for Cartesian simulation.');
        zmax=zmax*1000.
        modelArr=ak135model.reshape((21,7));
        modelArr[:,0]=modelArr[:,0]*1000.;
        modelArr[:,1]=modelArr[:,1]*1000.;
        modelArr[:,2]=modelArr[:,2]*1000.;
        modelArr[:,3]=modelArr[:,3]*1000.;
        for i in np.arange(21):
            if modelArr[i,0]==modelArr[i+1,0]:
                continue;
            if zmax <= modelArr[i+1,0]:
                if i==0:
                    rho=modelArr[0,1];
                    vp=modelArr[0,2];
                    vs=modelArr[0,3];
                    Qp=modelArr[0,4];
                    Qs=modelArr[0,5];
                    self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, z1=0, z2=zmax));
                else:
                    z1=modelArr[i,0];
                    rho=modelArr[i,1];
                    vp=modelArr[i,2];
                    vs=modelArr[i,3];
                    Qp=modelArr[i,4];
                    Qs=modelArr[i,5];
                    rho2=modelArr[i+1,1];
                    vp2=modelArr[i+1,2];
                    vs2=modelArr[i+1,3];
                    vpgrad=(vp2-vp)/(zmax-z1);
                    vsgrad=(vs2-vs)/(zmax-z1);
                    rhograd=(rho2-rho)/(zmax-z1);
                    self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad,
                            rhograd=rhograd, z1=z1, z2=zmax));
                break;
            z1=modelArr[i,0];
            rho=modelArr[i,1];
            vp=modelArr[i,2];
            vs=modelArr[i,3];
            Qp=modelArr[i,4];
            Qs=modelArr[i,5];
            z2=modelArr[i+1,0];
            rho2=modelArr[i+1,1];
            vp2=modelArr[i+1,2];
            vs2=modelArr[i+1,3];
            vpgrad=(vp2-vp)/(zmax-z1);
            vsgrad=(vs2-vs)/(zmax-z1);
            rhograd=(rho2-rho)/(zmax-z1);
            self.append(Blockmodel(vp=vp, vs=vs, rho=rho, Qp=Qp, Qs=Qs, vpgrad=vpgrad, vsgrad=vsgrad,
                            rhograd=rhograd, z1=z1, z2=z2));
        return;
        
    def Write2Input(self, infname):
        InputStr = "\n\n#------------------- block models: -------------------\n";
        blockstring = 'block vp=%f vs=%f rho=%f ';
        for block in self.blocks:
            tailstr='';
            if block.Qp!=None:
                if block.Qp<0:
                    raise ValueError('Negative Qp!');
                tailstr+=' qp=%f ' %(block.Qp);
            if block.Qs!=None:
                if block.Qs<0:
                    raise ValueError('Negative Qs!');
                tailstr+=' qs=%f ' %(block.Qs);
            if block.vpgrad!=None:
                tailstr+=' vpgrad=%f ' %(block.vpgrad);
            if block.vsgrad!=None:
                tailstr+=' vsgrad=%f ' %(block.vsgrad);
            if block.rhograd!=None:
                tailstr+=' rhograd=%f ' %(block.rhograd);
            if block.absdepth!=None:
                tailstr+=' absdepth=%d ' %(block.absdepth);
            if block.x1!=None:
                tailstr+=' x1=%f ' %(block.x1);
            if block.x2!=None:
                tailstr+=' x2=%f ' %(block.x2);
            if block.y1!=None:
                tailstr+=' y1=%f ' %(block.y1);
            if block.y2!=None:
                tailstr+=' y2=%f ' %(block.y2);
            if block.z1!=None:
                tailstr+=' z1=%f ' %(block.z1);
            if block.z2!=None:
                tailstr+=' z2=%f ' %(block.z2);
            
            headstr=( blockstring %(block.vp, block.vs, block.rho) );
            InputStr += headstr;
            InputStr += tailstr;
            InputStr += '\n';
        with open(infname, "a") as f:
            f.write(InputStr)
        f.close()
        return
    
    
### rfile model
class rBlock(object):
    """
    A class to hold rfile block data and header
    Modified from pysw4 by Shahar Shani-Kadmiel
    """
    def __init__(self, number=0, hh=None, hv=None, z0=0.0, nc=3, ni=None, nj=None, nk=None,
            xyextent=(), xzextent=(), yzextent = (), data=np.array([])):
        self.number = number;
        self.hh = hh;
        if hh !=None and hv==None:
            self.hv = hh;
        else:
            self.hv = hv;
        self.z0  = z0;
        self.nc  = nc;
        self.ni  = ni;
        self.nj  = nj;
        self.nk = nk;
        self.xyextent = xyextent;
        self.xzextent = xzextent;
        self.yzextent = yzextent;
        self.data = data;
        return;

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
        hh,hv,z0,nc,ni,nj,nk = np.fromfile(f, block_hdr_dtype, 1)[0];
        self.hh = hh;
        self.hv = hv;
        self.z0  = z0;
        self.nc  = nc;
        self.ni  = ni;
        self.nj  = nj;
        self.nk = nk;
        return hh,hv,z0,nc,ni,nj,nk;
    
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
            f.write(val);
        return;
        
class rModel(object):
    """
    A class to hold rfile header and blocks
    """
    def __init__(self, filename=None, magic=1, precision = 4, attenuation = 0, az = 135, lon0= -118, lat0=37,
                 proj_str  = '+proj=longlat +datum=WGS84 +no_defs', nb = 0, rblocks=None):
        self.filename  = filename;
        self.magic = magic;
        self.precision  = precision;
        self.attenuation    = attenuation;
        self.az  = az;
        self.lon0  = lon0;
        self.lat0 = lat0;
        self.proj_str = proj_str;
        self.nb = nb;
        self.rblocks = [];
        if isinstance(rblocks, rBlock):
            rblocks = [rblocks];
        if rblocks:
            self.rblocks.extend(rblocks);
        if filename!=None:
            self.read(filename);
        return;

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
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
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
    
    def AddSingleBlock(self, ni=0, nj=0, nk=0, hh=None, hv=None, z0=0.0, nc=3, xyextent=(), xzextent=(), yzextent = (), data=np.array([]),
        vs=None, vp=None, rho=None, Qs=None, Qp=None, vsgrad=None, vpgrad=None, rhograd=None, Qsgrad=None, Qpgrad=None):
        ya
        if self.nb == 0 and hh==None and hv == None :
            raise ValueError('Error: Please set gird spacing and initial depth!')
        if self.nb != 0 and hh==None and hv == None :
            hh=self.rblocks[-1].hh;
            hv=self.rblocks[-1].hv;
        
        
        rblock=rBlock(number=self.nb+1, hh=hh, hv=hv, z0=z0, nc=nc, ni=ni, nj=nj, nk=nk,
            xyextent=(), xzextent=(), yzextent = (), data=np.array([]));
        
        
        self.nb+=1;
    
    
    def read_hdr(self, f):
        (magic,precision,
            attenuation,az,lon0,lat0,mlen) = np.fromfile(f, hdr_dtype, 1)[0];
        proj_str_dtype = 'S' + str(mlen)
        proj_str = np.fromfile(f, proj_str_dtype, 1)[0]
        nb = np.fromfile(f, 'int32', 1)[0]
        self.magic=magic;
        self.precision  = precision;
        self.attenuation    = attenuation;
        self.az  = az;
        self.lon0  = lon0;
        self.lat0 = lat0;
        self.mlen=mlen;
        self.proj_str = proj_str;
        self.nb = nb;
        return magic,precision,attenuation,az,lon0,lat0,mlen,proj_str,nb;
    
    def write_hdr(self, f):
       """Write rfile header"""
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
            self.read_hdr(f);
            if verbose:
                print self;
                flush;
            for b in range(self.nb):
                rblock = rBlock()
                rblock.number = b+1
                rblock.read_block_hdr(f);
                self.append(rblock);
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
                self.rblocks[b].xyextent = (0                  , (nj-1)*hh*1e-3,
                                            0                  , (ni-1)*hh*1e-3)
                self.rblocks[b].xzextent = (0                  , (ni-1)*hh*1e-3,
                                            (z0+(nk-1)*hv)*1e-3, z0*1e-3       )
                self.rblocks[b].yzextent = (0                  , (nj-1)*hh*1e-3,
                                            (z0+(nk-1)*hv)*1e-3, z0*1e-3       )
        return

    def write(self, filename, verbose=False ):
        """
        Write rfile into Model object
        """
        self.filename = filename
        with open(filename, 'wb') as f:
            self.write_hdr(f);
            if verbose:
                print self;
                flush;
            for b in np.arange(self.nb):
                rblock = self.rblocks[b]; 
                rblock.write_block_hdr(f);
                if verbose:
                    print block
                    flush;
            for b in np.arange(self.nb):                
                z=self.rblocks[b].data.reshape(self.rblocks[b].data.size);
                z.tofile(f)
        return 
    
    
    
    
    
