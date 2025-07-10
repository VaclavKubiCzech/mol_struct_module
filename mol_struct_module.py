# -*- coding: utf-8 -*-
"""
Molecular structure module

Last version: 23rd June 2025
"""

import numpy as np
from mol_struct_module import geometry, xyz, fdf
from mol_struct_module import constants as const


class MolecularStructure:
    def __init__(self, atoms, coordinates, lattice_vectors=None, name="Unnamed"):
        self.atoms = atoms  # List of atomic numbers or symbols
        self.coordinates = np.array(coordinates)  # Nx3 array of xyz positions
        self.lattice_vectors = lattice_vectors
        if self.lattice_vectors is None:
            self.lattice_vectors = np.zeros((3, 3))  # Default to zero if not provided
        else:
            self.lattice_vectors = np.array(lattice_vectors)  # 3x3 array
        self.name = name # string
          
    def __repr__(self):
        return f"MolecularStructure(atoms={self.atoms},\ncoordinates={self.coordinates},\nlattice_vectors={self.lattice_vectors})"

    def get_formula(self):
        from collections import Counter
        counts = Counter(self.atoms)
        return ''.join(f"{el}-{counts[el]} " for el in sorted(counts))

    def center_of_mass(self):
        return geometry.center_of_mass(self.atoms, self.coordinates)
        
    def write_xyz(self, file_name=None, symbols=False):
        file_name = file_name or (self.name + '.xyz')
        xyz.write_xyz(self.atoms, self.coordinates, file_name=file_name, symbols=symbols)

    def write_fdf(self, file_name=None):
        file_name = file_name or (self.name + '.fdf')
        fdf.write_fdf(self.atoms, self.coordinates, self.lattice_vectors, file_name=file_name)    

    def translate(self, vector):
        # Translate the structure by the given vector.
        translated = geometry.translate_structure(self.coordinates, vector)
        self.coordinates = translated

    def translate_gc(self):
        # Translate to geometrical center.
        translated = geometry.translate_to_geometrical_center(self.coordinates)
        self.coordinates = translated

    def translate_z(self, value):
        # Translate the structure along the z-axis by a given value.
        vector = np.array([0,0,value])
        self.translate(vector)

    def rotate_straight(self, vector):
        # Rotate the structure so that the given vector aligns with the z-axis.
        straight_structure = geometry.rotate_structure_straight(self.coordinates, vector)
        self.coordinates = straight_structure
    
    def rotate_from_to(self, from_vec, to_vec, origin=(0, 0, 0)):
        # Rotate the structure so that from_vec aligns with to_vec, rotating around 'origin'.
        self.coordinates = geometry.rotate_from_to(self.coordinates, from_vec, to_vec, origin=origin)

    def rotate_around_axis(self, from_vec, to_vec, axis, origin=(0, 0, 0)):
        # Rotate the structure around a given axis such that from_vec aligns with to_vec
        self.coordinates = geometry.rotate_around_axis(self.coordinates, from_vec, to_vec, axis, origin=origin)

    def remove_site(self, list_atoms_indices):
        # Remove specified atomic sites from the structure.
        list_atoms_indices = [list_atoms_indices] if type(list_atoms_indices) is not list else list_atoms_indices

        mask = [i not in list_atoms_indices for i in range(len(self.atoms))]
        self.coordinates = self.coordinates[mask]
        self.atoms = [atom for i, atom in enumerate(self.atoms) if mask[i]]

    def add_atom(self, atom_type, position):
        # Add an atom of type atom_type at the specified position.
        self.atoms.append(atom_type if isinstance(atom_type, int) else const.ATOM_NUMBER_DICT[atom_type])
        self.coordinates = np.concatenate((self.coordinates, [np.array(position)]))
        
    def add_au_to_hollow(self, au1, au2, au3, sign):
        # Add an Au atom to a hollow site defined by three other Au atoms. Special case for pyramidal hollow sites.
        pos_au1, pos_au2, pos_au3 = self.coordinates[au1], self.coordinates[au2], self.coordinates[au3]
        auau_dist = 2.970 # angstroms
        triangle_height = np.sqrt(3)*auau_dist/2
        pyramide_height = np.sqrt(2/3)*auau_dist
        
        center = (pos_au1+pos_au2+pos_au3)/3
        normal = np.cross(pos_au1-pos_au2, pos_au1-pos_au3)*sign
        normal /= np.linalg.norm(normal)
        
        new_position = center + normal*pyramide_height
        self.add_atom(79, new_position)
    
    def view(self):
        # Use ase to visualize the structure.
        from ase.io import read
        import ase.visualize
        
        self.write_xyz(file_name="temp.xyz", symbols=True)
        atoms = read("temp.xyz", format="xyz")  # specify format explicitly
        ase.visualize.view(atoms)
        

class MolecularDynamics:
    def __init__(self, atoms, dynamics, lattice_vectors=None, name="Unnamed"):
        self.atoms = atoms  # List of atomic numbers or symbols
        self.dynamics = np.array(dynamics)  # MxNx3 array of positions
        self.lattice_vectors = lattice_vectors
        if self.lattice_vectors is None:
            self.lattice_vectors = np.zeros((3, 3))  # Default to zero if not provided
        else:
            self.lattice_vectors = np.array(lattice_vectors)  # 3x3 array
        self.name = name # string
    
    def __repr__(self):
        return f"MolecularStructure(atoms={self.atoms},\ncoordinates (step 1)={self.dynamics[0]},\nlattice_vectors={self.lattice_vectors})"

    def write_fdf(self, step, file_name=None):
        # Write a FDF file for a specific step in the dynamics.
        file_name = file_name or (self.name + '.fdf')
        fdf.write_fdf(self.atoms, self.dynamics[step], self.lattice_vectors, file_name=file_name)

    def write_cropped_fdf(self, step, file_name=None):
        # Write a cropped (cropped from electrodes, preserves only two Au/metal atoms connected to linkers) FDF file for a specific step in the dynamics.
        file_name = file_name or (self.name + '.fdf')
        new_atoms = []
        new_coordinates = []
        switch = False
        for i, atom in enumerate(self.atoms):
            try:
                if (self.atoms[i+1] != 79) and (self.atoms[i] == 79):
                    switch = True
                elif (self.atoms[i-1] == 79) and (self.atoms[i] == 79):
                    switch = False
                else:
                    switch = True
            except:
                switch = False
            if switch:
                new_atoms.append(atom)
                new_coordinates.append(self.dynamics[step][i])
        new_coordinates = np.array(new_coordinates)
        fdf.write_fdf(new_atoms, new_coordinates, self.lattice_vectors, file_name=file_name)

    def add_step(self, coordinates):
        # Add a new step to the dynamics. NEED TESTING
        self.dynamics = np.concatenate((self.dynamics, [coordinates]), axis=0)

    def write_ani(self, file_name=None, symbols=True):
        # Write the molecular dynamics to an ANI file. NEED TESTING
        file_name = file_name or (self.name + '.ANI')
        xyz.write_ani(self.atoms, self.dynamics, file_name=file_name, symbols=symbols)


def read_fdf_geometry(file_name: str) -> MolecularStructure:
    atoms, coordinates, lattice_vectors = fdf.read_fdf_geometry(file_name)
    return MolecularStructure(atoms, coordinates, lattice_vectors, file_name)


def read_xyz_geometry(file_path: str) -> MolecularStructure:
    atoms, coordinates, lattice_vectors = xyz.read_xyz_geometry(file_path)   
    return MolecularStructure(atoms, coordinates, lattice_vectors, file_path)

def read_ani(file_name: str) -> MolecularDynamics:
    atoms, coordinates, lattice_vectors, name = xyz.read_ani(file_name)
    return MolecularDynamics(atoms, coordinates, lattice_vectors, name)

def lattice_copies(molecule: MolecularStructure, n: int, m: int, l: int, name: str = None) -> MolecularStructure:
    """
    Creates copies of a molecular structure by translating it in a lattice defined by n, m, l.
    n, m, l are the number of copies in each direction.
    The new coordinates are calculated by adding the lattice vectors multiplied by the indices.

    Parameters:
    molecule (MolecularStructure): The original molecular structure.
    n (int): Number of copies in the x-direction.
    m (int): Number of copies in the y-direction.
    l (int): Number of copies in the z-direction.
    name (str, optional): Name for the new molecular structure. If not provided, it will be generated based on the original name.

    Returns:
    MolecularStructure: A new molecular structure containing the copies.
    """
    
    new_coordinates = []
    new_atoms = []
    
    for i in range(n):
        for j in range(m):
            for k in range(l):
                for s, site in enumerate(molecule.coordinates):
                    new_coordinates.append(site + 
                                           i*molecule.lattice_vectors[0] +
                                           j*molecule.lattice_vectors[1] +
                                           k*molecule.lattice_vectors[2])
                    new_atoms.append(molecule.atoms[s])
        
    new_coordinates = np.array(new_coordinates)    
    if not name:
        new_name = molecule.name + "_multiple"
    else:
        new_name = name
        
    return MolecularStructure(new_atoms, new_coordinates, molecule.lattice_vectors, new_name)


def xyz2fdf(file_name: str, out_file=None):
    # Convert XYZ file to FDF format.
    mol = read_xyz_geometry(file_name)
    if out_file:
        mol.name =  out_file
    mol.write_fdf()
    

def fdf2xyz(file_name: str, out_file=None):
    # Convert FDF file to XYZ format.
    mol = read_fdf_geometry(file_name)
    if out_file:
        mol.name =  out_file
    mol.write_xyz()


def connect_two_structures(mol1, mol2):
    """
    Connects two molecular structures by concatenating their coordinates and atoms.
    Assumes that mol1 and mol2 have compatible lattice vectors.
    """
    new_coordinates = np.concatenate((mol1.coordinates, mol2.coordinates))
    new_atoms = []
    for site in mol1.atoms:
        new_atoms.append(site)
    for site in mol2.atoms:
        new_atoms.append(site)
    new_lattice = mol1.lattice_vectors
    name = mol1.name + "_" + mol2.name
    return MolecularStructure(new_atoms, new_coordinates, new_lattice, name)


def connect_molecule_with_cluster(mol, bottom, top):
    """
    Connects a molecule with two electrodes (bottom and top) by aligning the molecule's gold sites with the electrodes' connection sites.
    Assumes that the molecule has two gold/silver sites at the ends.
    The bottom electrode is translated to the bottom gold site of the molecule, and the top electrode is positioned accordingly.
    The bottom electrode's last atom is the binding site for the molecule, and the top electrode's first atom is the binding site for the molecule.

    Parameters:
    mol (MolecularStructure): The molecule to be connected.
    bottom (MolecularStructure): The bottom electrode structure.
    top (MolecularStructure): The top electrode structure.
    Returns:
    MolecularStructure: A new molecular structure that combines the molecule and the electrodes.
    """
    
    # translate molecule so that bottom gold is connected to gold site of bottom electrode
    trans = bottom.coordinates[-1] - mol.coordinates[0]
    mol.translate(trans)
    
    # translate top electrode to proper position so that molecule (with golds) length correspond to connection sites distance
    dist = np.linalg.norm(mol.coordinates[0]-mol.coordinates[-1])
    diff = top.coordinates[0] - bottom.coordinates[-1]
    xysquared = diff[0]**2 + diff[1]**2
    z1 = np.sqrt(dist**2 - xysquared) - diff[2]
    top.translate_z(z1)
    
    # rotate molecule so the other gold matches gold of top connection site molecule
    origin = bottom.coordinates[-1]
    to_vec = top.coordinates[0] - bottom.coordinates[-1]
    from_vec = mol.coordinates[-1] - mol.coordinates[0]
    mol.rotate_from_to(from_vec, to_vec, origin)
    
    # rotate molecule around previous to_vec to optimal position given by au trimer (not optimal for ad atoms!)
    from_vec = mol.coordinates[1] - mol.coordinates[0]
    to_vec = 2*bottom.coordinates[-1] - bottom.coordinates[-2] - bottom.coordinates[-3]
    axis = mol.coordinates[-1]-mol.coordinates[0]
    origin = mol.coordinates[0]
    mol.rotate_around_axis(from_vec, to_vec, axis, origin)

    # remove metal sites from molecule
    mol.remove_site([0,len(mol.atoms)-1])
    
    # create the whole system
    temp_join = connect_two_structures(bottom, mol)
    final_join = connect_two_structures(temp_join, top)
    final_join.lattice_vectors[2][2] = final_join.coordinates[-1][2]
    return final_join


def standard_lattice():
    # returns standard lattice vectors (hexagonal) with third vector being zero
    latt = np.array([[5.939775,  -10.287988,   0.000000],
                     [5.939775,   10.287988,   0.000000],
                     [0.000000,    0.000000,   0.000000]])
    
    return latt


def recognize_lattice_layer(coordinates, tol = 0.05):
    # recognizes the layer of a given set of coordinates based on the average z-coordinate.
    # coordinates is either one vector or a array of vectors.
    def layer(x,y):
        a = 5.939775/4
        b = 10.287988/4
        if (abs(np.array([x,y])/np.array([a,b]) - np.round(np.array([x,y])/np.array([a,b]))) < np.array([tol, tol])).all():
            return 'A'
        elif (abs(np.array([x-a,y-b/3])/np.array([a,b]) - np.round(np.array([x-a,y-b/3])/np.array([a,b]))) < np.array([tol, tol])).all():
            return 'B'
        elif (abs(np.array([x-a,y+b/3])/np.array([a,b]) - np.round(np.array([x-a,y+b/3])/np.array([a,b]))) < np.array([tol, tol])).all():
            return 'C'
        else:
            return None
        
    if np.array(coordinates).ndim == 1:
        return layer(coordinates[0], coordinates[1])
    elif np.array(coordinates).ndim == 2:
        list_layers = []
        for coord in coordinates:
            layer_type = layer(coord[0], coord[1])
            if layer_type:
                list_layers.append(layer_type)
        return list_layers


def return_lattice_positions(layer):
    # returns 16 coordinates of either of A, B, C layers
    a1 = np.array([5.939775, -10.287988])/4
    a2 = np.array([5.939775,  10.287988])/4
    
    lattice = []
    if layer == 'A':
        for i in range(4):
            for j in range(4):
                lattice.append(i * a1 + j * a2)
    elif layer == 'B':
        for i in range(4):
            for j in range(4):
                lattice.append(i * a1 + j * a2 + a1/3 + 2*a2/3)
    elif layer == 'C':
        for i in range(4):
            for j in range(4):
                lattice.append(i * a1 + j * a2 + 2*a1/3 + a2/3)

    return np.array(lattice)


def read_ani(file_name):
    # Reads an ANI file and returns a MolecularDynamics object.
    atoms, dynamics, lattice_vectors, name = xyz.read_ani(file_name)
    return MolecularDynamics(atoms, dynamics, lattice_vectors, name)

# === Example Usage ===

if __name__=="__main__":
    #read fdf file
    mol = read_fdf_geometry("STRUCT.fdf")
    
    #print molecular info
    print("Molecular Structure:")
    print(mol)
    print("\nCenter of Mass:")
    print(mol.center_of_mass())
    
    #create multiple copies by lattice and write its .xyz
    multiple = lattice_copies(mol, n=3, m=3, l=3)
    multiple.write_xyz()