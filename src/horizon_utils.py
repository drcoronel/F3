import pathlib
import numpy as np
import pandas as pd


# def read_horizons_3D(path):

#     with open(path) as f:
#         lines = f.readlines()
    
#     dfs =[]
#     for l in lines[4:-1]:
#         chars = l.split()
#         if len(chars) != 6:
#             pass
#         df = pd.DataFrame()
#         df["Surface"] = chars[2]
#         df["X"] = chars[3]
#         df["Y"] = chars[4]
#         df["Z"] = chars[5]
#         dfs.append(df)
    
#     return  pd.concat(dfs)



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


    surfaces = np.array(name)    
    xyt = np.array([X,Y,Z]).astype(np.float)
    xyt = xyt.T
    
    return  xyt , surfaces.squeeze()


def sort_horizons(xyt,surfaces):

    key_horizons = list(np.unique(surfaces))

    horizon_dict = {}

    for k in key_horizons:
        idx = np.where(surfaces == k)[0]
        x,y,z = xyt[idx,:].T
        horizon_dict[k] = [x,y,z]

    return horizon_dict

def get_horizon(horizon_dict,horizon_key):

    return horizon_dict[horizon_dict]
