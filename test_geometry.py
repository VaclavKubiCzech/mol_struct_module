import unittest
import numpy as np
from mol_struct_module import geometry

class TestCenterOfMass(unittest.TestCase):
    def test_center_of_mass_h2(self):
        # H2 molecule: two hydrogen atoms 1 Å apart on x-axis
        atoms = ['H', 'H']
        coordinates = np.array([[0.0, 0.0, 0.0],
                             [1.0, 0.0, 0.0]])
        
        # Expected center of mass is at 0.5 Å on x-axis
        expected_com = np.array([0.5, 0.0, 0.0])
        
        com = geometry.center_of_mass(atoms, coordinates)
        
        np.testing.assert_allclose(com, expected_com, atol=1e-6)

    def test_center_of_mass_with_atomic_numbers(self):
        # Water molecule approximation: O at origin, Hs at ±1 along x
        atoms = [8, 1, 1]
        coordinates = np.array([[0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [-1.0, 0.0, 0.0]])

        com = geometry.center_of_mass(atoms, coordinates)

        # Expected COM is closer to O due to higher mass
        self.assertTrue(np.allclose(com[0], 0.0, atol=0.1))  # roughly centered

if __name__ == '__main__':
    unittest.main()
