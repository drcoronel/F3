import pathlib
import numpy as np
import xarray as xr
import string

class Head:
    """
    A head class that holds the different fields from a Head
    file.
    """
    
    def __init__(self):
        pass
    
    def get_headsg_from_file(self,file):
        self.file = file
        for k,v in read_well_head(self.file).items():
            key = k.translate({ord(c): None for c in string.whitespace})
            setattr(self,key,v)

class Well(Head):
    """
    Well class inherits attributes from Head. 
    
    Methods:

    - get_head: Users a Head object, and the Well.Name to extract the well head
    from a Head object

    """

    def __init__(self,name,file):
        self.file = file
        self.path = pathlib.Path(file)
        self.Name = name
    
    def get_head(self,Head):
        
        idx = Head.Name.index(self.Name)
        head_dict = {k:v[idx] for k,v in Head.__dict__.items()}

        for k,v in head_dict.items():
            setattr(self,k,v)        

    def read_logs(self,logs_file):
        self.logs_file = logs_file
        self.logs = well_reader(self.logs_file)



def well_reader(path,nan_value = 1e30):
    """
    Function that reads a .data well logs file exported that from OpendTect

    returns : xr.DataArray that contains the corresponding well logs, and units as attribute
    """
    raw_names = np.loadtxt(path, delimiter='\t',max_rows=1,dtype=str)
    names = [n.split('(')[0] for n in raw_names]
    units = [n.split('(')[-1].split(')')[0] for n in raw_names]

    logs = np.loadtxt(path,delimiter = '\t', skiprows=1) 
    attributes = {names[i]:units[i] for i in range(len(names)) }
    depth = logs[:,0]
    
    temp = xr.DataArray(data = logs, coords=[depth,names],attrs=attributes)
    arr = temp.where(temp != nan_value)

    return arr

def read_well_head(file):
    

    """Reads a head file and returns a dictionary with the Head fields as keys
    and the fields values for each well as values
    """

    path = pathlib.Path(file)
    with open(path) as f:
        lines = f.readlines()

    # Ignore non relevant lines
    start = lines.index('BEGIN HEADER\n') 
    end = lines.index('END HEADER\n')
    
    head = lines[start + 1: end]

    # Make a list with the head fields
    raw_fields = [l.split('\n')[0] for l in head]
    fields = []  # Make a cleaner list of the head fields.
    for l in raw_fields:
        chars = l.split(',')
        if chars[0] == 'INT' or chars[0] =='FLOAT':
            fields.append(chars[1])
        else:
            fields.append(l)
    
    # Take the values part of the file, after the head is finished
    val_lines = lines[end + 1:]
    vals = [ l.split(' ')[:-1] for l in val_lines] # Clean the values list, it is space separated

    # Make the dictionary {field: [well1,well2,well3.....]}
    head_dict = {}
    for i,f in enumerate(fields):
        val = []
        for v in vals:
            val.append(v[i].replace('"',''))
        head_dict[f]=val
    


    return head_dict
