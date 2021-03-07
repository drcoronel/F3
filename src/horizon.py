import pathlib
import numpy as np
import pandas as pd


class Horizons:

    """
    Horizons class:

    methods :
        - load_horizons: load bulk horizons based on the parameters defined by
        the dictionary params. The result is 1 or 2 dictionaries with surface
        names as keys, and ixt and/or xyt as values
        -get_horizon: dictionary with surface name as key and ixt and/or xyt as values
        -plot_horizon: grid the xyt or ixt values of a horizon, and returns a plot.
    
    attributes:
        -horizons_ixt: [dict] horizon data in (3,n) arrays with IL/XL format
        -horizons_xyt: [dict] horizon data in (3,n) arrays with X/Y format
    
    """
    def __init__(self):
        pass

    def load_horizons(self, file,params):
        """
        params: a dictionary that contains:
            origin : software of origin of the horizons. Right now just "OpendTect"
            spatial : True if the file contains X/Y coordinates
            il_xl : True if the file contains IL/XL coordinates 
    
        """
        self.path = pathlib.Path(file)

        origin = params['origin']
        spatial = params['spatial']
        il_xl = params['il_xl']

        if (origin == "OpendTect") and (spatial == True) and (il_xl == False):
            
            xyt, surfaces = parse_horizons_3D(self.path,params)
            self.horizons_xyt = sort_horizons(xyt, surfaces)
            self.keys = np.unique(surfaces)
        
        if (origin == "OpendTect") and (spatial == True) and (il_xl == True):

            xyt, ixt,surfaces = parse_horizons_3D(self.path,params)
            self.horizons_ixt = sort_horizons(ixt, surfaces)
            self.horizons_xyt = sort_horizons(xyt, surfaces)
            self.keys = np.unique(surfaces)

        if (origin == "OpendTect") and (spatial == False) and (il_xl == True):
            
            ixt,surfaces = parse_horizons_3D(self.path,params)
            self.horizons_ixt = sort_horizons(ixt, surfaces)
            self.keys = np.unique(surfaces)


    def get_horizon(self, horizon_name,kind = 'ixt'):
        horizon = Horizon()
        horizon.name = horizon_name

        if kind == 'ixt':
            horizon.data_ixt = self.horizons_ixt[horizon_name]
        if kind == 'xyt':
            horizon.data_xyt = self.horizons_xyt[horizon_name]
        if kind == 'both':
            horizon.data_ixt = self.horizons_ixt[horizon_name]
            horizon.data_xyt = self.horizons_xyt[horizon_name]

        return horizon

    def plot_horizon(self, horizon_name,kind = 'ixt'):
        import matplotlib.pyplot as plt

        horizon = self.get_horizon(horizon_name,kind = kind)
        horizon.plot_horizon(kind = kind)


class Horizon(Horizons):

    """
    Horizon class:

    Inherits some methods from Horizons class
    methods :
        - load_horizons: load bulk horizons based on the parameters defined by
        the dictionary params. The result is 1 or 2 dictionaries with surface
        names as keys, and ixt and/or xyt as values
        -get_horizon: dictionary with surface name as key and ixt and/or xyt as values
        -plot_horizon: grid the xyt or ixt values of a horizon, and returns a plot.
    
    attributes:
        -horizons_ixt: [dict] horizon data in (3,n) arrays with IL/XL format
        -horizons_xyt: [dict] horizon data in (3,n) arrays with X/Y format
    
    """


    def __init__(self):
        self.seismic_attribute = {}
        pass

    def grid_horizon(self, kind ='ixt', nx=200, ny=200):
        from scipy.interpolate import griddata

        if kind == 'ixt':
            xs, ys, zs = self.data_ixt
        if kind == 'xyt':
            xs, ys, zs = self.data_xyt

        xi = np.linspace(xs.min(), xs.max(), nx)
        yi = np.linspace(ys.min(), ys.max(), ny)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((xs, ys), zs, (X, Y))

        return X, Y, Z

    def plot_horizon(self,kind = 'ixt',nx=200,ny=200):
        import matplotlib.pyplot as plt

        X, Y, Z = self.grid_horizon(kind = kind,nx=nx,ny=ny)

        fig, ax = plt.subplots()
        c = ax.pcolormesh(X, Y, Z, cmap="terrain_r")
        fig.colorbar(c, orientation="vertical", label="Z")
        ax.legend()
        ax.set_title(f"{self.name}")

        plt.show()
    
    def extract_seismic(self,seismic,twt_range=None,method='raw',attribute='seismic'):
        """

        args:
        seismic: Seismic object
        twt_range: a tuple with (twt.min,twt.max). If non, it is calculated from the horizon
        method : method to extract the amplitude. Raw (no interpolation) or Spl (Spline interp)
        attribute : if 'seismic': Extract seismic amplitude, if 'envelope', extract envelope attribute.

        returns: 
        update self.seismic_attribute dictionary.
        ndarray (4,N). Ilines,Xlines,Z,Seismic.

        """

        il = seismic.seismic.ilines.tolist()
        xl = seismic.seismic.xlines.tolist()
        tl = seismic.seismic.samples.tolist()
        
        horizon = self.data_ixt
        
        if twt_range == None:
            twt_range = (self.data_ixt[2].min()-50,self.data_ixt[2].max()+50) 
            idx_min = np.abs(tl - twt_range[0]).argmin()
            idx_max = np.abs(tl - twt_range[1]).argmin()

        else:
            idx_min = np.abs(tl - twt_range[0]).argmin()
            idx_max = np.abs(tl - twt_range[1]).argmin()
        
        twt = tl[idx_min:idx_max+1]
        if attribute == 'seismic':
            data = seismic.seismic.data[...,idx_min:idx_max+1]
        elif attribute == 'envelope':
            from .seismic import seis_envelope
            data = seis_envelope(seismic.seismic.data[...,idx_min:idx_max+1])
        
        hor_extr = np.zeros((horizon.shape[1],4))
        for i in range(horizon.shape[1]):
            ii_idx = il.index(int(horizon[0,i]))
            xx_idx = xl.index(int(horizon[1,i]))
            zz_idx = np.abs(twt - horizon[2,i]).argmin()
            if method == 'raw':
                amp = data[ii_idx,xx_idx,zz_idx].flatten()
            elif method =='spl':
                from scipy.interpolate import splev,splrep
                trace = data[ii_idx,xx_idx,:].flatten()
                temp = splrep(twt,trace)
                amp = splev(horizon[2,i],temp)
            
            hor_extr[i,0] = horizon[0,i]
            hor_extr[i,1] = horizon[1,i]
            hor_extr[i,2] = horizon[2,i]
            hor_extr[i,3] = amp
        temp = {attribute:hor_extr.T}
        self.seismic_attribute.update(temp)
        return hor_extr.T

 



def parse_horizons_3D(path,params):

    """
    This function parse a .dat file that containes different 3D horizons exported from 
    OpendTect

    params: a dictionary that contains:
        origin : software of origin of the horizons. Right now just "OpendTect"
        spatial : True if the file contains X/Y coordinates
        il_xl : True if the file contains IL/XL coordinates 

    returns: xyt (3,n),surface array   
    """
    origin = params['origin']
    spatial = params['spatial']
    il_xl = params['il_xl']

    if (origin == "OpendTect") and (spatial == True) and (il_xl == False):
        with open(path) as f:
            lines = f.readlines()

        name, X, Y, Z = [], [], [], []
        for l in lines[4:-1]:
            chars = l.split()
            # if len(chars) != 6:
            #     pass
            if (l == lines[0]) or (l == lines[1]) or (l == lines[2]) or (l == lines[3]):
                continue

            name.append(chars[2].rstrip('"'))
            X.append(chars[3].rstrip("'"))
            Y.append(chars[4].rstrip("'"))
            Z.append(chars[5].rstrip("'"))

        np.set_printoptions(suppress=True)
        surfaces = np.array(name)
        xyt = np.array([X, Y, Z]).astype(np.float)
        xyt = xyt.T

        return xyt, surfaces.squeeze()
    
    if (origin == "OpendTect") and (spatial == True) and (il_xl == True):
        
        with open(path) as f:
            lines = f.readlines()

        name, X, Y, il, xl , Z = [], [], [], [], [], []
        for l in lines:
            chars = l.split() 
            name.append(chars[2].rstrip('"'))
            X.append(chars[3].rstrip("'"))
            Y.append(chars[4].rstrip("'"))
            il.append(chars[5].rstrip("'"))
            xl.append(chars[6].rstrip("'"))
            Z.append(chars[7].rstrip("'"))


        np.set_printoptions(suppress=True)
        surfaces = np.array(name)
        xyt = np.array([X, Y, Z]).astype(np.float)
        xyt = xyt.T
        ixt = np.array([il,xl,Z]).astype(np.float)
        ixt = ixt.T

        return xyt, ixt, surfaces.squeeze()

    if (origin == "OpendTect") and (spatial == False) and (il_xl == True):
        
        with open(path) as f:
            lines = f.readlines()

        name, il, xl , Z = [], [], [], []
        for l in lines:
            chars = l.split() 

            name.append(chars[2].rstrip('"'))
            il.append(chars[3].rstrip("'"))
            xl.append(chars[4].rstrip("'"))
            Z.append(chars[5].rstrip("'"))

        np.set_printoptions(suppress=True)
        surfaces = np.array(name)
        ixt = np.array([il,xl,Z]).astype(np.float)
        ixt = ixt.T

        return ixt, surfaces.squeeze()
        


def sort_horizons(xyt, surfaces):

    """
    Function that takes the returns from the horizon parser
    and sort the horizons by surfaces in a dictionary.
    inputs: 
        xyt: (3,n) array representing the data of a horizon, it can be ixt or xyt
        surfaces: np.array with the surfaces names (n,)
    
    returns:
        horizon_dict: {surface1: arr1,surface2: arr2,....}
    """

    key_horizons = list(np.unique(surfaces))

    horizon_dict = {}

    for k in key_horizons:
        idx = np.where(surfaces == k)[0]
        arr = xyt[idx, :].T
        horizon_dict[k] = arr

    return horizon_dict
