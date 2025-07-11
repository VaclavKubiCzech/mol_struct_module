
from mol_struct_module import *
import numpy as np

atoms = ['H', 'H']
coordinates = np.array([[0.0, 0.0, 0.0],
                         [1.0, 0.0, 0.0]])  

hydrogen = MolecularStructure(atoms, coordinates)
hydrogen.write_xyz("hydrogen_test.xyz", symbols=True)
hydrogen.write_fdf("hydrogen_test.fdf")


#hydrogen = mol_struct_module.read_xyz_geometry("hydrogen_test.xyz")
h2 = read_fdf_geometry("hydrogen_test.fdf")
#h2.remove_site([0])
print(h2)
h2.rotate_from_to(np.array([1, 0, 0]), np.array([0, 0, 1]))

h2_ani = MS2MD(h2)

h2_ani.add_step(np.array([[0.0, 0.0, 0.0],
                         [1.0, 0.0, 1.0]]))

h2_ani.write_ani(file_name="test.ANI")
