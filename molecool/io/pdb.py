"""

"""
import numpy as np

def open_pdb(file_location):
    """
    Open and read coordinates and symbols from a pdb file.

    The pdb file must specify the atom elements in the last column, 
    and follow the conventions outlline in the PDB format specification.

    Parameters
    ---------

    """
    # This function reads in a pdb file and returns the atom names and coordinates.
    with open(file_location) as f:
        data = f.readlines()
    coordinates = []
    symbols = []
    for line in data:
        if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
            sym.append(line[76:79].strip())
            atom_coordinates = [float(x) for x in line[30:55].split()]
            coordinates.append(atom_coordinates)

    # convert list to numpy array
    coordinates = np.array(coordinates)
    return symbols, coordinates