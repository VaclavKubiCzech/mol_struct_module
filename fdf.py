'''
FDF MODULE --- DEFINES OPERATIONS FOR FDF FILES
'''

import numpy as np
import typing
from mol_struct_module import constants as const

atomdict = {'Au':' 1',
            'C':' 2',
            'N':' 3',
            'S':' 3',
            'Sn':' 3',
            'P':' 3',
            'H':' 4',
            'Au_surf':' 5',
            'O':' 6',
            'Cl':' 6',
            'F':' 6',
            'Fe':' 6',
            'Ni':' 8',
            'Ag':' 1'} 


def read_fdf_geometry(file_name: str) -> tuple:
    """
    Reads a FDF file and returns the atomic symbols, coordinates, and lattice vectors.

    Parameters:
    file_name (str): Path to the FDF file.

    Returns:
    tuple: A tuple containing a list of atomic symbols, a numpy array of coordinates, and a numpy array of lattice vectors.
    """
    lattice_vectors = []
    coordinates = []
    species_dict = {}
    atoms = []

    with open(file_name, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.lower().startswith('%block latticevectors'):
            i += 1
            while not lines[i].lower().strip().startswith('%endblock'):
                lattice_vectors.append(list(map(float, lines[i].split())))
                i += 1

        elif line.lower().startswith('%block chemicalspecieslabel'):
            i += 1
            while not lines[i].lower().strip().startswith('%endblock'):
                parts = lines[i].split()
                species_index = int(parts[0])
                atomic_number = int(parts[1])
                atomic_name = parts[2]
                species_dict[species_index] = const.ATOM_NUMBER_DICT[atomic_name]
                i += 1

        elif line.lower().startswith('%block atomiccoordinatesandatomicspecies'):
            i += 1
            while not lines[i].lower().strip().startswith('%endblock'):
                parts = lines[i].split()
                coords = list(map(float, parts[:3]))
                species_index = int(parts[3])
                atomic_number = species_dict.get(species_index, 0)
                atoms.append(atomic_number)
                coordinates.append(coords)
                i += 1
        i += 1

    if file_name.find("/"):
        name = file_name[-file_name[::-1].find("/")-1:-4]
    else:
        name = file_name[:-4]

    return atoms, coordinates, lattice_vectors


def write_fdf(atoms: typing.List[typing.Union[str, int]], coordinates: np.ndarray, lattice_vectors: np.ndarray, file_name: str):
    """
    Writes a FDF file.       

    Parameters:
    atoms (list): of atomic symbols/atomic numbers
    coordinates (ndarray): with xyz coordinates of contained atoms
    lattice_vectors (ndarray): with lattice vectors of the structure
    file_name (str): name of the output fdf file
    """
    
    if file_name.__contains__(".fdf"):
        fout = open(file_name,'w')  
    else:
        fout = open(file_name + '.fdf','w')

    #header
    fout.write(f"NumberOfAtoms:         {int(coordinates.size/3)}" + "\n")
    fout.write(f"NumberOfSpecies:        {len(set(atoms))}" + "\n")
    fout.write("\n")
        
    #lattice
    if lattice_vectors.any():
        fout.write("LatticeConstant 1. Ang " + "\n")
        fout.write("%block LatticeVectors" + "\n")
        fout.write(f"{lattice_vectors[0,0]:15.9f}   {lattice_vectors[0,1]:15.9f}   {lattice_vectors[0,2]:15.9f}" + "\n")
        fout.write(f"{lattice_vectors[1,0]:15.9f}   {lattice_vectors[1,1]:15.9f}   {lattice_vectors[1,2]:15.9f}" + "\n")
        fout.write(f"{lattice_vectors[2,0]:15.9f}   {lattice_vectors[2,1]:15.9f}   {lattice_vectors[2,2]:15.9f}" + "\n")
        fout.write("%endblock LatticeVectors" + "\n")
    else:
        fout.write("#No LatticeConstant 1. Ang " + "\n")
        fout.write("#No %block LatticeVectors" + "\n")
        
    fout.write("\n")
    fout.write("\n")
    fout.write("\n")
        
    # chemical specification
    fout.write("%block ChemicalSpeciesLabel  # Species index, atomic number, species label" + "\n")
    for atom in atomdict:
        if const.ATOM_NUMBER_DICT[atom] in atoms:
            # surface metallic atoms are written with number multiplied by 10
            if const.ATOM_NUMBER_DICT[atom] <= 250:
                fout.write(f"    {atomdict[atom]}   {const.ATOM_NUMBER_DICT[atom]}  {atom}                      " + "\n")
            else: # in case of surface metalic atoms
                fout.write(f"    {atomdict[atom]}   {const.ATOM_NUMBER_DICT[atom]//10}  {atom}                      " + "\n")
        elif atom in atoms:
            fout.write(f"    {atomdict[atom]}   {const.ATOM_NUMBER_DICT[atom]}  {atom}                      " + "\n")
    fout.write("%endblock ChemicalSpeciesLabel " + "\n")
        
    fout.write("\n")
    fout.write(f"#  {file_name}" + "\n")
    fout.write("\n")
    fout.write("\n")
        
    # coordinates
    fout.write("AtomicCoordinatesFormat  Ang " + "\n")
    fout.write("%block AtomicCoordinatesAndAtomicSpecies " + "\n") 
    for i, site in enumerate(coordinates):
        #first remove previous number from end of atomic labels to find out atom type
        atomtype = atoms[i]
        if type(atomtype) is str:
            atomnr = atomdict[atomtype]
        if type(atomtype) is int:
            atomnr = atomdict[const.ATOM_SYMBOL_DICT[atomtype]]
      
        line= "{a:15.9f}".format(a=site[0]) + " "
        line+="{a:15.9f}".format(a=site[1]) + " "
        line+="{a:15.9f}".format(a=site[2]) + " "      
        line+= atomnr +'\n'
        fout.write(line)
    fout.write("%endblock AtomicCoordinatesAndAtomicSpecies " + "\n")  
    fout.close()
