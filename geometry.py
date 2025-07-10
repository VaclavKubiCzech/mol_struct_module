'''
GEOMETRY MODULE --- DEFINES GEOMETRIC OPERAITIONS WITH ATOMIC STRUCTURES 
'''

import typing
import numpy as np


def center_of_mass(atoms: typing.List[typing.Union[str, int]], coordinates: np.ndarray) -> np.ndarray:
    """
    Calculates the center of mass of a molecular structure given atomic symbols and their coordinates.
    """
    from periodictable import elements
    if type(atoms[0]) is str:
        masses = [elements.symbol(el).mass for el in atoms]
    else:
        masses = [elements[el].mass for el in atoms]
    total_mass = sum(masses)
    com = sum(m * coord for m, coord in zip(masses, coordinates)) / total_mass
    return com


# TRANSLATION FUNCTIONS
def translate_structure(coordinates: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """
    Translates the coordinates of a molecular structure by a given vector.
    """
    new_coordinates = coordinates.copy()
    for i, site in enumerate(coordinates):
        new_coordinates[i] = site + vector
    return new_coordinates


def translate_to_geometrical_center(coordinates: np.ndarray) -> np.ndarray:
    """
    Translates the coordinates of a molecular structure to its geometrical center.
    """
    x_aux, y_aux, z_aux = 0, 0, 0
    for i, site in enumerate(coordinates):
        x_aux += site[0]
        y_aux += site[1]
        z_aux += site[2]
    x_aux /= i+1
    y_aux /= i+1
    z_aux /= i+1
    
    new_coordinates = translate_structure(coordinates, -np.array([x_aux,y_aux,z_aux]))
    return new_coordinates


# ROTATION FUNCTIONS
def rotate_from_to(coordinates: np.ndarray, from_vec: np.ndarray, to_vec: np.ndarray, origin=(0, 0, 0)):
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
    new_coordinates = coordinates.copy()
    new_coordinates = ((new_coordinates - origin) @ R.T) + origin
    return new_coordinates


def rotate_around_axis(coordinates: np.ndarray, from_vec: np.ndarray, to_vec: np.ndarray, axis: np.ndarray, origin=(0, 0, 0)):
    """
    Rotate the molecule around a given axis such that from_vec aligns with to_vec
    (only via rotation around 'axis'). Axis intersects the 'origin'.
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
    new_coordinates = coordinates.copy()
    new_coordinates = ((new_coordinates - origin) @ R.T) + origin
    return new_coordinates


def rotational_matrix(axis: str, angle: float) -> np.ndarray:
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
    #print(matrix)
    return matrix
                

def rotate_structure_straight(coordinates: np.ndarray, vector: np.ndarray):
    # transforms coordinates so that the vector is in z-axis after transformation
    vector = vector/np.linalg.norm(vector)
    x, y, z = vector[0], vector[1], vector[2]
    phi = np.sign(y)*np.arccos(x/np.sqrt(x**2+y**2))
    theta = np.arccos(z/np.sqrt(x**2+y**2+z**2))
    
    rotmat_z = rotational_matrix("z", -phi)
    rotmat_y = rotational_matrix("y", -theta)
    
    eps = 1e-4
    unit_test = np.matmul(rotmat_y, np.matmul(rotmat_z, vector))
    if (np.abs(unit_test[0]) >= eps) or (np.abs(unit_test[1]) >= eps):
        print("Error: rotate_structure_strainght finction does not work properly")
        print(unit_test, np.linalg.norm(unit_test))

    new_coordinates = coordinates.copy()
    for i, site in enumerate(coordinates):
        site_new = np.matmul(rotmat_y, np.matmul(rotmat_z, site))
        new_coordinates[i] = site_new
    return new_coordinates
