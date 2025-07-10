'''
XYZ MODULE --- DEFINES OPERATIONS FOR XYZ FILES
'''

import numpy as np
import typing
from constants import *


def read_xyz_geometry(file_path: str) -> tuple:
    """
    Reads an XYZ file and returns the atomic symbols, coordinates, and lattice vectors.
    
    Parameters:
    file_path (str): Path to the XYZ file.
    
    Returns:
    tuple: A tuple containing a list of atomic symbols, a numpy array of coordinates, and a numpy array of lattice vectors.
    """
    atoms = []
    coordinates = []
    lattice_vectors = np.zeros((3, 3))  # No lattice vectors in XYZ

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines[2:]:  # Skip first two lines
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        symbol = parts[0]
        coords = list(map(float, parts[1:4]))
        try:
            atoms.append(ATOM_NUMBER_DICT[symbol])
        except:
            atoms.append(int(symbol))
        coordinates.append(coords)
    
    if file_path.find("/"):
        name = file_path[-file_path[::-1].find("/")-1:-4]
    else:
        name = file_path[:-4]

    return atoms, coordinates, lattice_vectors


def read_ani(file_name, name="trajectory"):
    """
    Reads an ANI file and returns the atomic symbols, coordinates, and lattice vectors.
    Parameters:
    file_name (str): Path to the ANI file.
    name (str): Optional name for the structure, defaults to the file name without extension.
    
    Returns:
    tuple: A tuple containing a list of atomic symbols, a numpy array of coordinates, and a numpy array of lattice vectors.
    """
    atoms = []
    dynamics = []
    lattice_vectors = np.zeros((3, 3))
    name = name or file_name[:-4]  # Default name if not provided

    with open(file_name, 'r') as f:
        lines = f.readlines()

    num_atoms = int(lines[0].strip())
    
    for step in range(int(len(lines)/(num_atoms+2))):
        i = 0
        while i < len(lines):
            i += 1  # Skip comment line
            i += 1
    
            step_coords = []
            for _ in range(num_atoms):
                parts = lines[i].strip().split()
                symbol = parts[0]
                coords = list(map(float, parts[1:4]))
                step_coords.append(coords)
                if len(dynamics) == 0:
                    try:
                        atoms.append(ATOM_NUMBER_DICT[symbol])
                    except:
                        atoms.append(int(symbol))
                i += 1
            step_coords = np.array(step_coords)
        dynamics.append(step_coords)
    dynamics = np.array(dynamics)
    return atoms, dynamics, lattice_vectors, name
    

def write_xyz(atoms: typing.List[typing.Union[str, int]], coordinates: np.ndarray, file_name: str, symbols: bool): 
    """
    Writes a xyz file.

    Parameters:
    atoms (list): of atomic symbols/atomic numbers 
    coortinates (ndarray): with xyz coordinates of contained atoms
    file_name (str): name of the output xyz file
    symbols (bool): true if the first column of xyz file should be written as symbol (Au, C, H), false if atomic numbers (79, 6, 1)

    """
    if file_name.__contains__(".xyz"):
        fn = open(file_name,'w')  
    else:
        fn = open(file_name + '.xyz','w')
        
    fn.write(str(len(atoms)) + "\n")
    fn.write("\n")  #write obligatory blank line
    #fn.seek(0)    #rewind
      
    for i, site in enumerate(coordinates):
        if symbols:
            try:
                atom_symbol = ATOM_SYMBOL_DICT[atoms[i]]
            except KeyError:
                atom_symbol = str(atoms[i]) 
        else:
            try:
                atom_symbol = int(atoms[i])
            except KeyError:
                atom_symbol = ATOM_NUMBER_DICT[atoms[i]]
        try:
            if atom_symbol <= 250:
                fn.write(str(atom_symbol)+ "   ")
            else: # in case of surface metallic atoms
                fn.write(str(atom_symbol//10)+ "   ")
        except:
            fn.write(str(atom_symbol)+ "   ")
        fn.write("{:15.11f}".format(site[0]) + "   ")
        fn.write("{:15.11f}".format(site[1]) + "   ")
        fn.write("{:15.11f}".format(site[2]) + "   ")
        fn.write("\n")


def write_xyz(atoms: typing.List[typing.Union[str, int]], coordinates: np.ndarray, file_name: str, symbols: bool): 
    """
    Writes a xyz file.

    Parameters:
    atoms (list): of atomic symbols/atomic numbers 
    coortinates (ndarray): with xyz coordinates of contained atoms
    file_name (str): name of the output xyz file
    symbols (bool): true if the first column of xyz file should be written as symbol (Au, C, H), false if atomic numbers (79, 6, 1)

    """
    if file_name.__contains__(".xyz"):
        fn = open(file_name,'w')  
    else:
        fn = open(file_name + '.xyz','w')
        
    fn.write(str(len(atoms)) + "\n")
    fn.write("\n")  #write obligatory blank line
    #fn.seek(0)    #rewind
      
    for i, site in enumerate(coordinates):
        if symbols:
            try:
                atom_symbol = ATOM_SYMBOL_DICT[atoms[i]]
            except KeyError:
                atom_symbol = str(atoms[i]) 
        else:
            try:
                atom_symbol = int(atoms[i])
            except KeyError:
                atom_symbol = ATOM_NUMBER_DICT[atoms[i]]
        try:
            if atom_symbol <= 250:
                fn.write(str(atom_symbol)+ "   ")
            else: # in case of surface metallic atoms
                fn.write(str(atom_symbol//10)+ "   ")
        except:
            fn.write(str(atom_symbol)+ "   ")
        fn.write("{:15.11f}".format(site[0]) + "   ")
        fn.write("{:15.11f}".format(site[1]) + "   ")
        fn.write("{:15.11f}".format(site[2]) + "   ")
        fn.write("\n")


def write_xyz(atoms: typing.List[typing.Union[str, int]], coordinates: np.ndarray, file_name: str, symbols: bool): 
    """
    Writes a xyz file.

    Parameters:
    atoms (list): of atomic symbols/atomic numbers 
    coortinates (ndarray): with xyz coordinates of contained atoms
    file_name (str): name of the output xyz file
    symbols (bool): true if the first column of xyz file should be written as symbol (Au, C, H), false if atomic numbers (79, 6, 1)

    """
    if file_name.__contains__(".xyz"):
        fn = open(file_name,'w')  
    else:
        fn = open(file_name + '.xyz','w')
        
    fn.write(str(len(atoms)) + "\n")
    fn.write("\n")  #write obligatory blank line
    #fn.seek(0)    #rewind
      
    for i, site in enumerate(coordinates):
        if symbols:
            try:
                atom_symbol = ATOM_SYMBOL_DICT[atoms[i]]
            except KeyError:
                atom_symbol = str(atoms[i]) 
        else:
            try:
                atom_symbol = int(atoms[i])
            except KeyError:
                atom_symbol = ATOM_NUMBER_DICT[atoms[i]]
        try:
            if atom_symbol <= 250:
                fn.write(str(atom_symbol)+ "   ")
            else: # in case of surface metallic atoms
                fn.write(str(atom_symbol//10)+ "   ")
        except:
            fn.write(str(atom_symbol)+ "   ")
        fn.write("{:15.11f}".format(site[0]) + "   ")
        fn.write("{:15.11f}".format(site[1]) + "   ")
        fn.write("{:15.11f}".format(site[2]) + "   ")
        fn.write("\n")
    fn.close()
        

def write_ani(atoms: typing.List[typing.Union[str, int]], dynamics: np.ndarray, file_name: str, symbols: bool): 
    """
    Writes a ANI file.

    Parameters:
    atoms (list): of atomic symbols/atomic numbers 
    coortinates (ndarray): with xyz coordinates of contained atoms
    file_name (str): name of the output xyz file
    symbols (bool): true if the first column of xyz file should be written as symbol (Au, C, H), false if atomic numbers (79, 6, 1)

    """
    if file_name.__contains__(".ANI"):
        fn = open(file_name,'w')  
    else:
        fn = open(file_name + '.ANI','w')

    for step in dynamics:
        fn.write(str(len(atoms)) + "\n")
        fn.write("\n")  #write obligatory blank line
        #fn.seek(0)    #rewind
        
        for i, site in enumerate(step):
            if symbols:
                try:
                    atom_symbol = ATOM_SYMBOL_DICT[atoms[i]]
                except KeyError:
                    atom_symbol = str(atoms[i]) 
            else:
                try:
                    atom_symbol = int(atoms[i])
                except KeyError:
                    atom_symbol = ATOM_NUMBER_DICT[atoms[i]]
            try:
                if atom_symbol <= 250:
                    fn.write(str(atom_symbol)+ "   ")
                else: # in case of surface metallic atoms
                    fn.write(str(atom_symbol//10)+ "   ")
            except:
                fn.write(str(atom_symbol)+ "   ")
            fn.write("{:15.11f}".format(site[0]) + "   ")
            fn.write("{:15.11f}".format(site[1]) + "   ")
            fn.write("{:15.11f}".format(site[2]) + "   ")
            fn.write("\n")
    fn.close()