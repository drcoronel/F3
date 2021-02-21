import pathlib
import numpy as np
import pandas as pd


class Facies:
    """
    Facies are discribe as point clouds of xyt
    """

    def __init__(self):
        pass

    def from_horizons(self,horizons,name_base,name_top):
        """
        Uses a Horizons object, to extract top (T) and base (B) horizons and from
        them, create a facies object (F). 
        
        X  X  T  T           X  X  F  F  
        T  T  X  X           F  F  F  F  
        X  X  X  X           F  F  F  F   
        B  B  B  B           F  F  F  F 
        """
        
        base_horizon = horizons.get_horizon(name_base)
        top_horizon = horizons.get_horizon(name_top)
    

    ### TODOOO
        pass 
        