import pathlib
import numpy as np
import pandas as pandas
import matplotlib.pyplot as plt
import segyio


class Seismic:
    def __init__(self):
        pass

    def load_seismic(self, file, iline=189, xline=193):

        self.path = pathlib.Path(file)
        self.seismic = segyio.open(self.path, iline=iline, xline=xline)
        inlines = self.seismic.ilines
        xlines = self.seismic.xlines
        samples = self.seismic.samples
        self.seismic.data = self.seismic.trace.raw[:].reshape((len(inlines),len(xlines),len(samples)))
        self.seismic.shape = self.seismic.data.shape

        self.bounds = self.calculate_bounds()

        #Calculate indices:
        il_idx = np.arange(0,len(inlines))
        xl_idx = np.arange(0,len(xlines))
        tl_idx = np.arange(0,len(samples))
        self.seismic.indices = np.array([il_idx,xl_idx,tl_idx])

        # # Add CDP_X, CDP_Y sorting

        # cdpx_ = self.bounds['CDP_X']
        # cdpy_ = self.bounds['CDP_Y']
        
        # cdpx = np.arange(cdpx_['min'],cdpx_['max'],step=cdpx_['step'])
        # cdpy = np.arange(cdpy_['min'],cdpy_['max'],step=cdpy_['step'])

        # self.seismic.CDP_X = cdpx
        # self.seismic.CDP_Y = cdpy
        # self.seismic.rotation = segyio.tools.rotation(self.seismic)



    def calculate_bounds(self):
        """
        Calculate seismic survey bounds. Assumes 4 corner points.
        corner points are given by the minimum and maximum of the ilines and xlines.

        returns: a dictionary that contains the trace indices, (iline,xline), (CDP_X,CDP_y)
        coordinates for each corner point.
        """
        inlines = self.seismic.ilines
        xlines = self.seismic.xlines

        # Find min and max ilines, xlines.
        IL_extremes = (inlines.min(), inlines.max())
        XL_extremes = (xlines.min(), xlines.max())

        ## Assume boundary has 4 points. Found the il,xl values of the boundary.
        lines_extremes = [
            (IL_extremes[i], XL_extremes[j]) for i in range(2) for j in range(2)
        ]

        trace_header = self.seismic.header

        # Find traces and CDP_X, CDP_Y locations

        bounds_dict = {"indices": [], "lines": lines_extremes, "CDP_coords": [],
                'ilines':{},
                'xlines':{},
                'CDP_X':{},
                'CDP_Y':{}           
        }

        for line in lines_extremes:
            for i, h in enumerate(trace_header):
                il = h.get(segyio.TraceField.INLINE_3D)
                xl = h.get(segyio.TraceField.CROSSLINE_3D)
                bound = (il, xl)
                if bound == line:
                    cdpx = h.get(segyio.TraceField.CDP_X)
                    cdpy = h.get(segyio.TraceField.CDP_Y)
                    bounds_dict["indices"].append(i)
                    bounds_dict["CDP_coords"].append((cdpx, cdpy))

        il_dict = {'min':inlines.min(), 'max': inlines.max(),
                   'step': np.diff(inlines).max() }
        xl_dict = {'min':xlines.min(), 'max': xlines.max(),
                   'step': np.diff(xlines).max() }  

        bounds_dict['ilines'].update(il_dict)
        bounds_dict['xlines'].update(xl_dict)


        cdpxs,cdpys = [],[]

        for h in trace_header:
            cdpxs.append(h.get(segyio.TraceField.CDP_X))
            cdpys.append(h.get(segyio.TraceField.CDP_Y))
        
        cdpx = np.unique(cdpxs)
        cdpy = np.unique(cdpys)

        I,X,_ = self.seismic.shape
        self.seismic.cdpx = np.array(cdpxs).reshape(I,X)
        self.seismic.cdpy = np.array(cdpys).reshape(I,X)

        cdpx_dict = {'min':cdpx.min(), 'max': cdpx.max(),
                   'step': np.diff(cdpx).max() }
        cdpy_dict = {'min':cdpy.min(), 'max': cdpy.max(),
                   'step': np.diff(cdpy).max() }

        bounds_dict['CDP_X'].update(cdpx_dict)
        bounds_dict['CDP_Y'].update(cdpy_dict)

        return bounds_dict

    # def extract_values(self,horizon,attribute = None):
    #     """
    #     This function is based on the work of Alessandro Amato del Monte
    #     https://github.com/aadm/geophysical_notes/blob/master/seismic_amplitude_extraction.ipynb
    #     """
    #     if attribute == None:
    #         data = self.seismic.data
        
    #     hor = horizon.data
    #     horizon_extraction = np.zeros((4,hor.size))
    #     seis_x,seis_y,seis_t = data
    #     hor_x, hor_y, hor_z = horizon.data

    #     for i in range(hor.shape[1]):

    #         xx_idx = self.seismic.CDP_X.tolist()





    def calculate_seismic_attribute(self,attribute_name = 'envelope',
        subset = None):

        if subset ==  None:
            data = self.seismic.seismic.data
        
        else: 
            data = self.seismic.seismic.data[subset]

        if 'envelope': 
            return seis_envelope(data)
        if 'phase':
            return seis_phase(data)
        if 'cosine_of_phase':
            return cosine_of_phase(data)





## Seismic attributes
from scipy.signal import hilbert as hilbert_transform


def seismic_hilbert_transform(seismic):
    
    hilbert = np.zeros_like(seismic)
    i,x,_=seismic.shape
    for ils in range(i):
        for xls in range(x): 
            hilbert[ils,xls,:] = np.imag(hilbert_transform(seismic[ils,xls,:])[:])

    return hilbert 

def seis_envelope(seismic):

    hilbert=seismic_hilbert_transform(seismic)
    z = seismic + 1j*hilbert
    return np.abs(z)

def seis_phase(seismic):

    hilbert=seismic_hilbert_transform(seismic)
    z = seismic + 1j*hilbert
    return np.angle(z)

def cosine_of_phase(seismic):

    phase=seis_phase(seismic)
    return np.cos(phase)
