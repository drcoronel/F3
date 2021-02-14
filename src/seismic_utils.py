import pathlib
import numpy as n
import pandas as pandas
import matplotlib.pyplot as plt
import segyio


class Seismic:
    def __init__(self):
        pass

    def load_seismic(self, file, iline=189, xline=193):

        self.path = pathlib.Path(file)
        self.seismic = segyio.open(self.path, iline=iline, xline=xline)

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

        bounds_dict = {"indices": [], "lines": lines_extremes, "CDP_coords": []}

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

        self.bounds = bounds_dict
