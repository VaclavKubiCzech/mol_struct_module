# -*- coding: utf-8 -*-
"""
Molecular structure module

Last version: 23rd June 2025
"""

import numpy as np


# ========================================================================================
# Definition of global dictionaries
# ========================================================================================
ATOM_NUMBER_DICT = {'Au':79, # dictionary of atomic symbols and their corresponding number
            'C':6,
            'N':7,
            'S':16,
            'Sn':50,
            'P':15,
            'H':1,
            'O':8,
            'Cl':17,
            'F':9,
            'Fe':26,
            'Ni':28,
            'Ag':47,
            'Au_surf':790}

ATOM_SYMBOL_DICT = {79:'Au', # dictionary of atomic symbols and their corresponding number
            6:'C',
            7:'N',
            16:'S',
            50:'Sn',
            15:'P',
            1:'H',
            8:'O',
            17:'Cl',
            9:'F',
            26:'Fe',
            28:'Ni',
            47:'Ag',
            790:'Au_surf'} 
# ========================================================================================


class MolecularStructure:
    def __init__(self, atoms, coordinates, lattice_vectors, name):
        self.atoms = atoms  # List of atomic numbers or symbols
        self.coordinates = np.array(coordinates)  # Nx3 array of positions
        self.lattice_vectors = np.array(lattice_vectors)  # 3x3 array
        self.name = name # string
          
    def __repr__(self):
        return f"MolecularStructure(atoms={self.atoms},\ncoordinates={self.coordinates},\nlattice_vectors={self.lattice_vectors})"

    def get_formula(self):
        from collections import Counter
        counts = Counter(self.atoms)
        return ''.join(f"{el}-{counts[el]} " for el in sorted(counts))

    def center_of_mass(self):
        # calculates center of mass, only if library Fdictable available
        from periodictable import elements
        if type(self.atoms[0]) is str:
            masses = [elements.symbol(el).mass for el in self.atoms]
        else:
            masses = [elements[el].mass for el in self.atoms]
        total_mass = sum(masses)
        com = sum(m * coord for m, coord in zip(masses, self.coordinates)) / total_mass
        return com
    
    def write_xyz(self, filename=None):
      if not filename:
          fn = open(self.name + '.xyz','w')
      else:
          fn = open(filename,'w')
      
      fn.write(str(len(self.atoms)) + "\n")
      fn.write("\n")  #read blank line
      #fn.seek(0)    #rewind
      
      for i, site in enumerate(self.coordinates):
        try:
            if self.atoms[i] <= 250:
                fn.write(str(self.atoms[i])+ "   ")
            else: # in case of surface metallic atoms
                fn.write(str(self.atoms[i]//10)+ "   ")
        except:
            fn.write(str(self.atoms[i])+ "   ")
        fn.write("{:15.11f}".format(site[0]) + "   ")
        fn.write("{:15.11f}".format(site[1]) + "   ")
        fn.write("{:15.11f}".format(site[2]) + "   ")
        fn.write("\n")
        
    def write_xyz_symbol(self, filename=None):
      if not filename:
          fn = open(self.name + '.xyz','w')
      else:
          fn = open(filename,'w')
      
      fn.write(str(len(self.atoms)) + "\n")
      fn.write("\n")  #read blank line
      #fn.seek(0)    #rewind
      
      for i, site in enumerate(self.coordinates):
        fn.write(str(ATOM_SYMBOL_DICT[self.atoms[i]])+ "   ")
        fn.write("{:15.11f}".format(site[0]) + "   ")
        fn.write("{:15.11f}".format(site[1]) + "   ")
        fn.write("{:15.11f}".format(site[2]) + "   ")
        fn.write("\n")

    def write_fdf(self, filename=None):
        #from periodictable import elements
        
        if not filename:
            fout = open(self.name + '.fdf','w')
        else:
            fout = open(filename,'w')
        
        
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
                    'Ag':' 1'}     #CHANGE FOR PROBLEM AT HAND

        
        #header
        fout.write(f"NumberOfAtoms:         {int(self.coordinates.size/3)}" + "\n")
        fout.write(f"NumberOfSpecies:        {len(set(self.atoms))}" + "\n")
        fout.write("\n")
        
        #lattice
        if self.lattice_vectors.any():
            fout.write("LatticeConstant 1. Ang " + "\n")
            fout.write("%block LatticeVectors" + "\n")
            fout.write(f"{self.lattice_vectors[0,0]:15.9f}   {self.lattice_vectors[0,1]:15.9f}   {self.lattice_vectors[0,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[1,0]:15.9f}   {self.lattice_vectors[1,1]:15.9f}   {self.lattice_vectors[1,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[2,0]:15.9f}   {self.lattice_vectors[2,1]:15.9f}   {self.lattice_vectors[2,2]:15.9f}" + "\n")
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
            if ATOM_NUMBER_DICT[atom] in self.atoms:
                if ATOM_NUMBER_DICT[atom] <= 250:
                    fout.write(f"    {atomdict[atom]}   {ATOM_NUMBER_DICT[atom]}  {atom}                      " + "\n")
                else: # in case of surface metalic atoms
                    fout.write(f"    {atomdict[atom]}   {ATOM_NUMBER_DICT[atom]//10}  {atom}                      " + "\n")
        fout.write("%endblock ChemicalSpeciesLabel " + "\n")
        
        fout.write("\n")
        fout.write(f"#  {self.name}" + "\n")
        fout.write("\n")
        fout.write("\n")
        
        # coordinates
        fout.write("AtomicCoordinatesFormat  Ang " + "\n")
        fout.write("%block AtomicCoordinatesAndAtomicSpecies " + "\n") 
        for i, site in enumerate(self.coordinates):
          #first remove previous number from end of atomic labels to find out atom type
          atomtype = self.atoms[i]
          if type(atomtype) is str:
              atomnr = atomdict[atomtype]
          if type(atomtype) is int:
              atomnr = atomdict[ATOM_SYMBOL_DICT[atomtype]]
      
          line= "{a:15.9f}".format(a=site[0]) + " "
          line+="{a:15.9f}".format(a=site[1]) + " "
          line+="{a:15.9f}".format(a=site[2]) + " "      
          line+= atomnr +'\n'
          fout.write(line)
        fout.write("%endblock AtomicCoordinatesAndAtomicSpecies " + "\n")  

    def rotate_straight(self, vector):
        mol = rotate_structure_straight(self, vector)
        self.coordinates = mol.coordinates

    def translate(self, vector):
        mol = translate_structure(self, vector)
        self.coordinates = mol.coordinates

    def translate_gc(self):
        #translate to geometrical center
        mol = translate_to_geometrical_center(self)
        self.coordinates = mol.coordinates

    def rotate_from_to(self, from_vec, to_vec, origin=(0, 0, 0)):
        """
        Rotate the molecule so that from_vec aligns with to_vec, rotating around 'origin'.
        """
        from_vec = np.array(from_vec) / np.linalg.norm(from_vec)
        to_vec = np.array(to_vec) / np.linalg.norm(to_vec)
        origin = np.array(origin)

        # Compute rotation axis and angle
        axis = np.cross(from_vec, to_vec)
        norm_axis = np.linalg.norm(axis)
        
        if norm_axis < 1e-8:
            # Vectors are parallel or anti-parallel
            if np.dot(from_vec, to_vec) > 0:
                return  # No rotation needed
            else:
                # 180-degree rotation around any perpendicular axis
                axis = np.eye(3)[np.argmin(np.abs(from_vec))]  # Pick orthogonal basis
                norm_axis = np.linalg.norm(axis)
        
        axis /= norm_axis
        angle = np.arccos(np.clip(np.dot(from_vec, to_vec), -1.0, 1.0))

        # Rodrigues' rotation formula
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

        # Apply rotation to all coordinates
        self.coordinates = ((self.coordinates - origin) @ R.T) + origin

    def translate_z(self, value):
        vector = np.array([0,0,value])
        self.translate(vector)
        
    def remove_site(self, list_atoms_indices):
        new_coordinates = []
        new_atoms = []
        
        if type(list_atoms_indices) is not list:
            list_atoms_indices = [list_atoms_indices]
        
        for i, site in enumerate(self.atoms):
            if i not in list_atoms_indices:
                new_coordinates.append(self.coordinates[i])
                new_atoms.append(site)
        self.coordinates = np.array(new_coordinates)
        self.atoms = new_atoms
        
    def rotate_around_axis(self, from_vec, to_vec, axis, origin=(0, 0, 0)):
        """
        Rotate the molecule around a given axis such that from_vec aligns with to_vec
        (only via rotation around 'axis'). Rotation is performed around the 'origin'.
        """
        from_vec = np.array(from_vec)
        to_vec = np.array(to_vec)
        axis = np.array(axis)
        origin = np.array(origin)

        # Normalize the axis
        axis = axis / np.linalg.norm(axis)

        # Project vectors onto the plane perpendicular to the axis
        def project_onto_plane(v, axis):
            return v - np.dot(v, axis) * axis

        v1_proj = project_onto_plane(from_vec, axis)
        v2_proj = project_onto_plane(to_vec, axis)

        # Normalize projections
        v1_norm = v1_proj / (np.linalg.norm(v1_proj) + 1e-12)
        v2_norm = v2_proj / (np.linalg.norm(v2_proj) + 1e-12)

        # Compute angle and direction
        angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))
        direction = np.dot(np.cross(v1_norm, v2_norm), axis)
        if direction < 0:
            angle = -angle

        # Rodrigues' rotation matrix around axis
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

        # Rotate coordinates around origin
        self.coordinates = ((self.coordinates - origin) @ R.T) + origin

    def add_atom(self, atom_type, position):
        # adding atom to structure 
        if type(atom_type) is int:
            self.atoms.append(atom_type)
        else:
            self.atoms.append(ATOM_NUMBER_DICT[atom_type])
        
        self.coordinates = np.concatenate((self.coordinates, [np.array(position)]))
                       
    def view(self):
        from ase.io import read
        import ase.visualize
        
        self.write_xyz_symbol(filename="temp.xyz")
        atoms = read("temp.xyz")  # or .xyz, .traj, etc.
        ase.visualize.view(atoms)
        
    def add_au_to_hollow(self, au1, au2, au3, sign):
        # adding au atom into hollow site of au layer, sign = +-1
        
        pos_au1 = self.coordinates[au1]
        pos_au2 = self.coordinates[au2]
        pos_au3 = self.coordinates[au3]
        
        auau_dist = 2.970 # angstroms
        triangle_height = np.sqrt(3)*auau_dist/2
        pyramide_height = np.sqrt(2/3)*auau_dist
        
        center = (pos_au1+pos_au2+pos_au3)/3
        normal = np.cross(pos_au1-pos_au2, pos_au1-pos_au3)*sign
        normal /= np.linalg.norm(normal)
        
        new_position = center + normal*pyramide_height
        
        self.add_atom(79, new_position)
        
        
        
        
        

class MolecularDynamics:
    def __init__(self, atoms, dynamics, lattice_vectors, name):
        self.atoms = atoms  # List of atomic numbers or symbols
        self.dynamics = np.array(dynamics)  # MxNx3 array of positions
        self.lattice_vectors = np.array(lattice_vectors)  # 3x3 array
        self.name = name # string
        
    
    def __repr__(self):
        return f"MolecularStructure(atoms={self.atoms},\ncoordinates (step 1)={self.dynamics[0]},\nlattice_vectors={self.lattice_vectors})"

    def write_fdf(self, step, filename=None):
        #from periodictable import elements
        
        if not filename:
            fout = open(self.name + '.fdf','w')
        else:
            fout = open(filename,'w')
        
        
        atomdict = {'Au':' 1',
                    'C':' 2',
                    'N':' 3',
                    'S':' 3',
                    'Sn':' 3',
                    'P':' 3',
                    'H':' 4',
                    'O':' 6',
                    'Cl':' 6',
                    'F':' 6',
                    'Fe':' 6',
                    'Ni':' 8',
                    'Ag':' 1'}     #CHANGE FOR PROBLEM AT HAND

        
        #header
        fout.write(f"NumberOfAtoms:         {int(self.dynamics[step].size/3)}" + "\n")
        fout.write(f"NumberOfSpecies:        {len(set(self.atoms))}" + "\n")
        fout.write("\n")
        
        #lattice
        if self.lattice_vectors.any():
            fout.write("LatticeConstant 1. Ang " + "\n")
            fout.write("%block LatticeVectors" + "\n")
            fout.write(f"{self.lattice_vectors[0,0]:15.9f}   {self.lattice_vectors[0,1]:15.9f}   {self.lattice_vectors[0,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[1,0]:15.9f}   {self.lattice_vectors[1,1]:15.9f}   {self.lattice_vectors[1,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[2,0]:15.9f}   {self.lattice_vectors[2,1]:15.9f}   {self.lattice_vectors[2,2]:15.9f}" + "\n")
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
            if ATOM_NUMBER_DICT[atom] in self.atoms:
                fout.write(f"    {atomdict[atom]}   {ATOM_NUMBER_DICT[atom]}  {atom}                      " + "\n")
        fout.write("%endblock ChemicalSpeciesLabel " + "\n")
        
        fout.write("\n")
        fout.write(f"#  {self.name}" + "\n")
        fout.write("\n")
        fout.write("\n")
        
        # coordinates
        fout.write("AtomicCoordinatesFormat  Ang " + "\n")
        fout.write("%block AtomicCoordinatesAndAtomicSpecies " + "\n") 
        for i, site in enumerate(self.dynamics[step]):
          #first remove previous number from end of atomic labels to find out atom type
          atomtype = self.atoms[i]
          if type(atomtype) is str:
              atomnr = atomdict[atomtype]
          if type(atomtype) is int:
              atomnr = atomdict[ATOM_SYMBOL_DICT[atomtype]]
      
          line= "{a:15.10f}".format(a=site[0]) + " "
          line+="{a:15.10f}".format(a=site[1]) + " "
          line+="{a:15.10f}".format(a=site[2]) + " "      
          line+= atomnr +'\n'
          fout.write(line)
        fout.write("%endblock AtomicCoordinatesAndAtomicSpecies " + "\n")  
        
        
    def write_cropped_fdf(self, step, filename=None):
        # TBD
        #from periodictable import elements
        
        if not filename:
            fout = open(self.name + '.fdf','w')
        else:
            fout = open(filename,'w')
        
        
        atomdict = {'Au':' 1',
                    'C':' 2',
                    'N':' 3',
                    'S':' 3',
                    'Sn':' 3',
                    'P':' 3',
                    'H':' 4',
                    'O':' 6',
                    'Cl':' 6',
                    'F':' 6',
                    'Fe':' 6',
                    'Ni':' 8',
                    'Ag':' 1'}     #CHANGE FOR PROBLEM AT HAND
        
        #header
        number_of_atoms = 2
        for atom in self.atoms:
            if (atom != 79):
                number_of_atoms += 1
        fout.write(f"NumberOfAtoms:         {number_of_atoms}" + "\n")
        fout.write(f"NumberOfSpecies:        {len(set(self.atoms))}" + "\n")
        fout.write("\n")
        
        #lattice
        if self.lattice_vectors.any():
            fout.write("LatticeConstant 1. Ang " + "\n")
            fout.write("%block LatticeVectors" + "\n")
            fout.write(f"{self.lattice_vectors[0,0]:15.9f}   {self.lattice_vectors[0,1]:15.9f}   {self.lattice_vectors[0,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[1,0]:15.9f}   {self.lattice_vectors[1,1]:15.9f}   {self.lattice_vectors[1,2]:15.9f}" + "\n")
            fout.write(f"{self.lattice_vectors[2,0]:15.9f}   {self.lattice_vectors[2,1]:15.9f}   {self.lattice_vectors[2,2]:15.9f}" + "\n")
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
            if ATOM_NUMBER_DICT[atom] in self.atoms:
                fout.write(f"    {atomdict[atom]}   {ATOM_NUMBER_DICT[atom]}  {atom}                      " + "\n")
        fout.write("%endblock ChemicalSpeciesLabel " + "\n")
        
        fout.write("\n")
        fout.write(f"#  {self.name}" + "\n")
        fout.write("\n")
        fout.write("\n")
            
        # coordinates
        fout.write("AtomicCoordinatesFormat  Ang " + "\n")
        fout.write("%block AtomicCoordinatesAndAtomicSpecies " + "\n")
        switch = False
        for i, site in enumerate(self.dynamics[step]):
            #first remove previous number from end of atomic labels to find out atom type
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
                atomtype = self.atoms[i]
                if type(atomtype) is str:
                    atomnr = atomdict[atomtype]
                if type(atomtype) is int:
                    atomnr = atomdict[ATOM_SYMBOL_DICT[atomtype]]
              
                line= "{a:15.10f}".format(a=site[0]) + " "
                line+="{a:15.10f}".format(a=site[1]) + " "
                line+="{a:15.10f}".format(a=site[2]) + " "      
                line+= atomnr +'\n'
                fout.write(line)
        fout.write("%endblock AtomicCoordinatesAndAtomicSpecies " + "\n")  


    def add_step(self, coordinates):
        new_dyn = []
        for dyn in self.dynamics:
            new_dyn.append(dyn)
        new_dyn.append(coordinates)
        new_dyn = np.array(new_dyn)
        self.dynamics = new_dyn
        
    
    def write_ani(self, filename=None):
        if not filename:
            fn = open(self.name + '.ANI','w')
        else:
            fn = open(filename,'w')
            
        for step in self.dynamics:
                    
            fn.write(str(len(self.atoms)) + "\n")
            fn.write("\n")  #read blank line
            #fn.seek(0)    #rewind
            
            for i, site in enumerate(step):
              fn.write(str(atom_symbol_dict[self.atoms[i]])+ "   ")
              fn.write("{:15.11f}".format(site[0]) + "   ")
              fn.write("{:15.11f}".format(site[1]) + "   ")
              fn.write("{:15.11f}".format(site[2]) + "   ")
              fn.write("\n")
        


def read_fdf_geometry(filename):
    lattice_vectors = []
    atomic_coords = []
    species_dict = {}
    atoms = []

    with open(filename, 'r') as file:
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
                species_dict[species_index] = ATOM_NUMBER_DICT[atomic_name]
                i += 1

        elif line.lower().startswith('%block atomiccoordinatesandatomicspecies'):
            i += 1
            while not lines[i].lower().strip().startswith('%endblock'):
                parts = lines[i].split()
                coords = list(map(float, parts[:3]))
                species_index = int(parts[3])
                atomic_number = species_dict.get(species_index, 0)
                atoms.append(atomic_number)
                atomic_coords.append(coords)
                i += 1
        i += 1
    
    if filename.find("/"):
        name = filename[-filename[::-1].find("/")-1:-4]
    else:
        name = filename[:-4]
    
    structure = MolecularStructure(atoms, atomic_coords, lattice_vectors, name)
    return structure


def read_xyz_geometry(filename):
    #from periodictable import elements
    
    atoms = []
    coordinates = []
    lattice_vectors = np.zeros((3, 3))  # No lattice vectors in XYZ

    with open(filename, 'r') as f:
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
    
    if filename.find("/"):
        name = filename[-filename[::-1].find("/")-1:-4]
    else:
        name = filename[:-4]

    return MolecularStructure(atoms, coordinates, lattice_vectors, name)


def create_copies_by_lattice(molecule, n, m, l, name=None):
    # creates same structure, n vinthin first lattice vector and so on
    
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


def rotational_matrix(axis, angle):
    if axis=="x":
        case = 0
    elif axis=="y":
        case = 1
    elif axis=="z":
        case = 2
        
    matrix = np.zeros((3,3))
    cos = np.cos(angle)
    sin = np.sin(angle)
    
    for i in range(3):
        for j in range(i+1):
            if (i == j) and (i == case):
                matrix[i,j] = 1
            elif (i == j) and (i != case):
                matrix[i,j] = cos
            elif (i != j) and (i != case) and (j != case):
                if case == 1:
                    matrix[i,j] = -sin
                    matrix[j,i] = -matrix[i,j]
                else:
                    matrix[i,j] = sin
                    matrix[j,i] = -matrix[i,j]
    print(matrix)
    return matrix
                

def rotate_structure_straight(molecule, vector):
    # transforms coordinates so that the vector is in z-axis after transformation
    vector = vector/np.linalg.norm(vector)
    x, y, z = vector[0], vector[1], vector[2]
    phi = np.sign(y)*np.arccos(x/np.sqrt(x**2+y**2))
    theta = np.arccos(z/np.sqrt(x**2+y**2+z**2))
    
    rotmat_z = rotational_matrix("z", -phi)
    rotmat_y = rotational_matrix("y", -theta)
    
    #CHECK
    eps = 1e-4
    unit_test = np.matmul(rotmat_y, np.matmul(rotmat_z, vector))
    if (np.abs(unit_test[0]) >= eps) or (np.abs(unit_test[1]) >= eps):
        print("Error: rotate_structure_strainght finction does not work properly")
        print(unit_test, np.linalg.norm(unit_test))
    
    new_coordinates = molecule.coordinates
    for i, site in enumerate(molecule.coordinates):
        site_new = np.matmul(rotmat_y, np.matmul(rotmat_z, site))
        new_coordinates[i] = site_new
    
    return MolecularStructure(molecule.atoms, new_coordinates, molecule.lattice_vectors, molecule.name)
    

def translate_structure(molecule, vector):
    # translates all sites by same vector
    
    new_coordinates = molecule.coordinates
    for i, site in enumerate(molecule.coordinates):
        new_coordinates[i] = site + vector
    
    return MolecularStructure(molecule.atoms, new_coordinates, molecule.lattice_vectors, molecule.name)
     

def translate_to_geometrical_center(molecule):
    # translates structure to its center of mass   
    
    x_aux, y_aux, z_aux = 0, 0, 0
    for i, site in enumerate(molecule.coordinates):
        x_aux += site[0]
        y_aux += site[1]
        z_aux += site[2]
    
    x_aux /= i+1
    y_aux /= i+1
    z_aux /= i+1
    
    new_mol = translate_structure(molecule, -np.array([x_aux,y_aux,z_aux]))
    
    return new_mol
           

def read_ani(filename, name="trajectory"):
    from periodictable import elements
    atoms = []
    dynamics = []
    lattice_vectors = np.zeros((3, 3))

    with open(filename, 'r') as f:
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
                    atoms.append(elements.symbol(symbol).number)
                i += 1
            step_coords = np.array(step_coords)
        dynamics.append(step_coords)
    dynamics = np.array(dynamics)
    return MolecularDynamics(atoms, dynamics, lattice_vectors, name)


def xyz2fdf(filename, outfile=None):
    mol = read_xyz_geometry(filename)
    if outfile:
        mol.name =  outfile
    mol.write_fdf()
    

def fdf2xyz(filename, outfile=None):
    mol = read_fdf_geometry(filename)
    if outfile:
        mol.name =  outfile
    mol.write_xyz()


def add_two_structures(mol1, mol2):
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
    # beware: molecule must have two gold/silver sites at the ends, functions for trimer tips
    
    # translate molecule so that it bottom gold is connected to gold site of bottom electrode
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
    
    # rotate molecule around previous to_vec to optimal position given by au trimer
    from_vec = mol.coordinates[1] - mol.coordinates[0]
    to_vec = 2*bottom.coordinates[-1] - bottom.coordinates[-2] - bottom.coordinates[-3]
    axis = mol.coordinates[-1]-mol.coordinates[0]
    origin = mol.coordinates[0]
    mol.rotate_around_axis(from_vec, to_vec, axis, origin)

    # remove metal sites from molecule
    mol.remove_site([0,len(mol.atoms)-1])
    
    
    # create the whole system
    temp_join = add_two_structures(bottom, mol)
    final_join = add_two_structures(temp_join, top)
    final_join.lattice_vectors[2][2] = final_join.coordinates[-1][2]
    return final_join

def standard_lattice():
    # returns standard lattice vectors (hexagonal) with third vector being zero
    
    latt = np.array([[5.939775,  -10.287988,   0.000000],
                     [5.939775,   10.287988,   0.000000],
                     [0.000000,    0.000000,   0.000000]])
    
    return latt

#function for returning A B C lattice vectors 4x4
''
def return_lattice_positions(layer):
    # returns 16 coordinates of either of A, B, C layers
    
    # standard lattice
    a1 = np.array([5.939775, -10.287988])/4
    a2 = np.array([5.939775,  10.287988])/4
    
    lattice = []
    
    match layer:
        case 'A':
            for i in range(4):
                for j in range(4):
                    lattice.append(i * a1 + j * a2)
            
        case 'B':
            for i in range(4):
                for j in range(4):
                    lattice.append(i * a1 + j * a2 + a1/3 + 2*a2/3)
            
        case 'C':
            for i in range(4):
                for j in range(4):
                    lattice.append(i * a1 + j * a2 + 2*a1/3 + a2/3)

    return np.array(lattice)
''

def read_ani(filename):
    # function enabling reading .ANI file for further manipulation

        
    atoms = []
    dynamics = []
    lattice_vectors = np.zeros((3, 3))  # No lattice vectors in XYZ

    # CHANGE (copied from read_xyz())
    with open(filename, 'r') as f:
        lines = f.readlines()

    num_atoms = int(lines[0].strip())

    not_last_step = True
    j = 0
    while not_last_step:
        # iterate over all atoms and write 
        coordinates = []
        for i in range(num_atoms):
            parts = lines[j*(num_atoms+2) + 2 + i].strip().split()
            if j == 0:
                symbol = parts[0]
                try:
                    atoms.append(ATOM_NUMBER_DICT[symbol])
                except:
                    atoms.append(int(symbol))
            coords = list(map(float, parts[1:4]))
            coordinates.append(coords)
        dynamics.append(coordinates)        
        j += 1
        try:
            t = lines[j*(num_atoms+2) + num_atoms] 
        except:
            not_last_step = False
    
    if filename.find("/"):
        name = filename[-filename[::-1].find("/")-1:-4]
    else:
        name = filename[:-4]

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
    multiple = create_copies_by_lattice(mol, n=3, m=3, l=3)
    multiple.write_xyz()