import pathlib
import numpy as np
import pandas as pd


class Horizons:
    def __init__(self):
        pass

    def load_horizons(self, file):

        self.path = pathlib.Path(file)
        ixt, surfaces = parse_horizons_3D(self.path)
        self.horizons = sort_horizons(ixt, surfaces)

    def get_horizon(self, horizon_name):
        horizon = Horizon()
        horizon.name = horizon_name
        horizon.data = self.horizons[horizon_name]
        return horizon

    def plot_horizon(self, horizon_name):
        import matplotlib.pyplot as plt

        horizon = self.get_horizon(horizon_name)
        horizon.plot_horizon()


class Horizon(Horizons):
    def __init__(self):
        pass

    def grid_horizon(self, nx=200, ny=200):
        from scipy.interpolate import griddata

        xs, ys, zs = self.data
        xi = np.linspace(xs.min(), xs.max(), nx)
        yi = np.linspace(ys.min(), ys.max(), ny)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((xs, ys), zs, (X, Y))

        return X, Y, Z

    def plot_horizon(self):
        import matplotlib.pyplot as plt

        X, Y, Z = self.grid_horizon()

        fig, ax = plt.subplots()
        c = ax.pcolormesh(X, Y, Z, cmap="terrain_r")
        fig.colorbar(c, orientation="vertical", label="Z")
        ax.legend()
        ax.set_title(f"{self.name}")

        plt.show()


def parse_horizons_3D(path):
    """
    This function parse a .dat file that containes different 3D horizons exported from 
    OpendTect

    returns: xyt (3,n),surface array   
    """
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


def sort_horizons(xyt, surfaces):

    key_horizons = list(np.unique(surfaces))

    horizon_dict = {}

    for k in key_horizons:
        idx = np.where(surfaces == k)[0]
        arr = xyt[idx, :].T
        horizon_dict[k] = arr

    return horizon_dict
