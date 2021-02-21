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
        
        if (origin == "OpendTect") and (spatial == True) and (il_xl == True):

            xyt, ixt,surfaces = parse_horizons_3D(self.path,params)
            self.horizons_ixt = sort_horizons(ixt, surfaces)
            self.horizons_xyt = sort_horizons(xyt, surfaces)

        if (origin == "OpendTect") and (spatial == False) and (il_xl == True):
            
            ixt,surfaces = parse_horizons_3D(self.path,params)
            self.horizons_ixt = sort_horizons(ixt, surfaces)


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

    def plot_horizon(self,kind = 'ixt'):
        import matplotlib.pyplot as plt

        X, Y, Z = self.grid_horizon(kind = kind)

        fig, ax = plt.subplots()
        c = ax.pcolormesh(X, Y, Z, cmap="terrain_r")
        fig.colorbar(c, orientation="vertical", label="Z")
        ax.legend()
        ax.set_title(f"{self.name}")

        plt.show()


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
