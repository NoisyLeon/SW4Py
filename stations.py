import os, shutil
import numpy as np


class StaInfo(object):
    """
    An object contains a station information several methods for station related analysis.
    -----------------------------------------------------------------------------------------------------
    General Parameters:
    stacode     - station name
    network     - network
    chan        - channels for analysis
    x,y,z     - position for station
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self, stacode=None, network='SW4', x=None, y=None, z=0,
            usgsformat=0, sacformat=1, writeEvery=1000, nsew=1, variables='displacement'):

        self.stacode=stacode;
        self.network=network;
        self.x=x;
        self.y=y;
        self.z=z;
        self.usgsformat=usgsformat;
        self.sacformat=sacformat;
        self.writeEvery=writeEvery;
        self.nsew=nsew;
        self.variables=variables;
        return;
    
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of Traces in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self

    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """        
        f = open(stafile, 'r')
        Sta=[]
        for lines in f.readlines():
            lines=lines.split()
            stacode=lines[0];
            network=lines[1];
            x=float(lines[2]);
            y=float(lines[3]);
            z=float(lines[4]);
            usgsformat=float(lines[5]);
            sacformat=float(lines[6]);
            writeEvery=float(lines[7]);
            nsew=float(lines[8]);
            variables=lines[9];
            netsta=network+'.'+stacode
            if Sta.__contains__(netsta):
                index=Sta.index(netsta)
                if abs(self[index].x-x) >0.01 and abs(self[index].y-y) and abs(self[index].z-z)>0.01:
                    raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                else:
                    print 'Warning: Repeated Station:' +netsta+' in Station List!'
                    continue
            Sta.append(netsta);
            self.append(StaInfo (stacode=stacode, network=network, x=x, y=y, z=z, usgsformat=usgsformat,
                sacformat=sacformat, writeEvery=writeEvery, nsew=nsew, variables=variables ));
            f.close()
        return
    
    def WriteStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x y z
        """
        f = open(stafile, 'w')
        for sta in self.stations:
            tempstr='%s %s %1.5e %1.5e %1.5e %d %d %d %d %s  \n' %( sta.stacode, sta.network, sta.x, sta.y, sta.z,
                sta.usgsformat, sta.sacformat, sta.writeEvery, sta.nsew, sta.variables)
            f.writelines(tempstr)
        f.close()
        return
    
    def AddSingle(self, x, y, stacode=None, network='SW4',  z=0,
            usgsformat=0, sacformat=1, writeEvery=1000, nsew=1, variables='displacement'):
        if stacode==None:
            stacode=str(int(x/1000))+'S'+str(int(y/1000));
        self.append(StaInfo (stacode=stacode, network=network, x=x, y=y, z=z, usgsformat=usgsformat,
                sacformat=sacformat, writeEvery=writeEvery, nsew=nsew, variables=variables ));
        return;
    
    def HomoStaLst(self, xmin, Nx, dx, ymin, Ny, dy, network='SW4', z=0,
            usgsformat=0, sacformat=1, writeEvery=1000, nsew=1, variables='displacement'):
        for i in np.arange(Nx):
            for j in np.arange(Ny):
                x=xmin+dx*i;
                y=ymin+dy*j;
                stacode=str(int(i))+'S'+str(int(j));
                self.append(StaInfo (stacode=stacode, network=network, x=x, y=y, z=z, usgsformat=usgsformat,
                    sacformat=sacformat, writeEvery=writeEvery, nsew=nsew, variables=variables ));
        return;
    
    def LineStaLst(self, xmin=None, Nx=None, dx=None, y=None, ymin=None, Ny=None, dy=None, x=None, network='SW4', z=0,
            usgsformat=0, sacformat=1, writeEvery=1000, nsew=1, variables='displacement'):
        if xmin!=None and Nx!=None and dx!=None and y!=None:
            for i in np.arange(Nx):
                x=xmin+dx*i;
                j=int(y/1000);
                stacode=str(int(i))+'S'+str(int(j));
                self.append(StaInfo (stacode=stacode, network=network, x=x, y=y, z=z, usgsformat=usgsformat,
                    sacformat=sacformat, writeEvery=writeEvery, nsew=nsew, variables=variables ));
        elif ymin!=None and Ny!=None and dy!=None and x!=None:
            for j in np.arange(Ny):
                y=ymin+dy*j;
                i=int(x/1000);
                stacode=str(int(i))+'S'+str(int(j));
                self.append(StaInfo (stacode=stacode, network=network, x=x, y=y, z=z, usgsformat=usgsformat,
                    sacformat=sacformat, writeEvery=writeEvery, nsew=nsew, variables=variables ));
        return;

    def Wrtie2Input(self, infname, cartesian=True):
        if cartesian==True:
            sac_string = 'rec x=%.3f y=%.3f depth=%.3f file=%s sta=%s writeEvery=%d usgsformat=%d sacformat=%d nsew=%d variables=%s\n';
        else:
            sac_string = 'rec lon=%.10f lat=%.10f depth=%.3f file=%s sta=%s writeEvery=%d usgsformat=%d sacformat=%d nsew=%d variables=%s\n';
        InputStr = "\n\n#------------------- stations at locations: -------------------\n";
        for station in self.stations:
            netsta=station.network+'.'+station.stacode;
            InputStr += (sac_string %( station.x, station.y, station.z, netsta, station.stacode, station.writeEvery,
                station.usgsformat, station.sacformat, station.nsew, station.variables));
        InputStr += "\n\n#------------------- end of stations -------------------\n";
        with open(infname, "a") as f:
            f.write(InputStr)
        f.close()
        return 
