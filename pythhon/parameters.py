from enum import Enum, auto
import numpy as np
import pathlib
import shutil
import os

# real = np.float64
real = float
# real = np.float32
# intg = np.uint8
intg = int
size_type = np.uint8

debug_mode = -1

class BoundaryType(Enum):
    DISPLACEMENT = "DISPLACEMENT"
    PRESSURE = "PRESSURE"
    SLIDE = "SLIDE"


class ShapeType(Enum):
    POINT = "POINT"
    SEGMENT = "SEGMENT"
    TRIANGLE = "TRIANGLE"
    QUADRANGLE = "QUADRANGLE"
    POLYGON = "POLYGON"
    TETRAHEDRON = "TETRAHEDRON"
    PRISM = "PRISM"
    PYRAMID = "PYRAMID"
    HEXAHEDRON = "HEXAHEDRON"


class QuadratureType(Enum):
    GAUSS = "GAUSS"


class BasisType(Enum):
    MONOMIAL = "MONOMIAL"


class ElementType(Enum):
    HDG_HIGH = "HDG_HIGH"
    HHO_HIGH = "HHO_HIGH"
    HDG_EQUAL = "HDG_EQUAL"
    HHO_EQUAL = "HHO_EQUAL"
    HDG_LOW = "HDG_LOW"
    HHO_LOW = "HHO_LOW"


class FieldType(Enum):
    SCALAR = "SCALAR"
    VECTOR = "VECTOR"
    # DISPLACEMENT_SMALL_STRAIN = "DISPLACEMENT_SMALL_STRAIN"
    # DISPLACEMENT_FINITE_STRAIN = "DISPLACEMENT_FINITE_STRAIN"
    # DISPLACEMENT_SMALL_STRAIN_PLANE_STRAIN = "DISPLACEMENT_SMALL_STRAIN_PLANE_STRAIN"
    # DISPLACEMENT_SMALL_STRAIN_PLANE_STRESS = "DISPLACEMENT_SMALL_STRAIN_PLANE_STRESS"
    # DISPLACEMENT_FINITE_STRAIN_PLANE_STRAIN = "DISPLACEMENT_FINITE_STRAIN_PLANE_STRAIN"
    # DISPLACEMENT_FINITE_STRAIN_PLANE_STRESS = "DISPLACEMENT_FINITE_STRAIN_PLANE_STRESS"
    DISPLACEMENT = "DISPLACEMENT"
    DISPLACEMENT_PLANE_STRAIN = "DISPLACEMENT_PLANE_STRAIN"
    DISPLACEMENT_PLANE_STRESS = "DISPLACEMENT_PLANE_STRAIN"


class StrainTempType(Enum):
    SMALL_STRAIN = ("SMALL_STRAIN",)
    FINITE_STRAIN = ("FINITE_STRAIN",)

class StressType(Enum):
    PIOLA_KIRCHOFF_1 = "PK1"
    CAUCHY = "CAUCHY"

class StrainType(Enum):
    DISPLACEMENT_TRANSFORMATION_GRADIENT = "F"
    DISPLACEMENT_SYMMETRIC_GRADIENT = "GRAD_SYM_U"


class DerivationType(Enum):
    SYMMETRIC = "SYMMETRIC"
    FULL = "FULL"


# def get_strain_type(field_type: FieldType):
#     if field_type in [
#         FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRESS,
#         FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRAIN,
#         FieldType.DISPLACEMENT_FINITE_STRAIN,
#     ]:
#         return StrainTempType.FINITE_STRAIN
#     else:
#         return StrainTempType.SMALL_STRAIN

def get_project_path():
    return pathlib.Path(__file__).parent.parent.absolute()

def get_res_file_path(res_file_name: str, suffix: str):
    project_path = get_project_path()
    return os.path.join(project_path, "res/{}_{}.txt".format(res_file_name, suffix))



# voigt_vector_plane = {(0, 0): 0, (1, 1): 1, (0, 1): 2, (1, 0): 3}
# # voigt_vector_plane = {(0, 0): 0, (1, 1): 3, (0, 1): 1, (1, 0): 2}
# voigt_scalar_plane = {(0, 0): 0, (0, 1): 1}
