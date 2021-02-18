from pythhon.parameters import *


class Field:
    label: str
    field_type: FieldType
    stress_type: StressType
    strain_type: StrainType
    derivation_type: DerivationType

    def __init__(
        self,
        label: str,
        euclidean_dimension:int,
        field_type: FieldType,
        stress_type: StressType,
        strain_type: StrainType,
        derivation_type: DerivationType = DerivationType.FULL,
    ):
        """

        Args:
            label:
            euclidean_dimension:
            field_type:
            stress_type:
            strain_type:
            derivation_type:
        """
        self.label = label
        self.field_type = field_type
        self.stress_type = stress_type
        self.strain_type = strain_type
        self.euclidean_dimension = euclidean_dimension
        if self.field_type == FieldType.SCALAR:
            self.derivation_type = derivation_type
            self.field_dimension = 1
            self.gradient_dimension = euclidean_dimension
            self.voigt_indices = {
                (0, 0): 0,
                (0, 1): 1
            }
            self.voigt_coefficients = {
                (0, 0): 1.0,
                (0, 1): 1.0
            }
        elif self.field_type in [FieldType.VECTOR, FieldType.DISPLACEMENT]:
            self.derivation_type = derivation_type
            if self.derivation_type == DerivationType.FULL:
                self.field_dimension = euclidean_dimension
                self.gradient_dimension = euclidean_dimension * euclidean_dimension
                if euclidean_dimension == 2:
                    self.voigt_indices = {
                        (0, 0): 0,
                        (1, 1): 1,
                        (0, 1): 2,
                        (1, 0): 3
                    }
                    self.voigt_coefficients = {
                        (0, 0): 1.0,
                        (1, 1): 1.0,
                        (0, 1): 1.0,
                        (1, 0): 1.0,
                    }
                elif euclidean_dimension == 3:
                    self.voigt_indices = {
                        (0, 0): 0,
                        (1, 1): 1,
                        (2, 2): 2,
                        (0, 1): 3,
                        (1, 0): 4,
                        (0, 2): 5,
                        (2, 0): 6,
                        (1, 2): 7,
                        (2, 1): 8,
                    }
                    self.voigt_coefficients = {
                        (0, 0): 1.0,
                        (1, 1): 1.0,
                        (2, 2): 1.0,
                        (0, 1): 1.0,
                        (1, 0): 1.0,
                        (0, 2): 1.0,
                        (2, 0): 1.0,
                        (1, 2): 1.0,
                        (2, 1): 1.0,
                    }
            elif self.derivation_type == DerivationType.SYMMETRIC:
                self.field_dimension = euclidean_dimension
                self.gradient_dimension = euclidean_dimension + (euclidean_dimension ** 2 - euclidean_dimension)/2
                if euclidean_dimension == 2:
                    self.voigt_indices = {
                        (0, 0): 0,
                        (1, 1): 1,
                        (0, 1): 2
                    }
                    self.voigt_coefficients = {
                        (0, 0): 1.0,
                        (1, 1): 1.0,
                        (0, 1): np.sqrt(2.0)
                    }
                elif euclidean_dimension == 3:
                    self.voigt_indices = {
                        (0, 0): 0,
                        (1, 1): 1,
                        (2, 2): 2,
                        (0, 1): 3,
                        (0, 2): 4,
                        (1, 2): 5,
                    }
                    self.voigt_coefficients = {
                        (0, 0): 1.0,
                        (1, 1): 1.0,
                        (2, 2): 1.0,
                        (0, 1): np.sqrt(2.0),
                        (0, 2): np.sqrt(2.0),
                        (1, 2): np.sqrt(2.0),
                    }
        elif self.field_type in [FieldType.DISPLACEMENT_PLANE_STRAIN, FieldType.DISPLACEMENT_PLANE_STRESS]:
            if euclidean_dimension != 2:
                raise ValueError("The euclidean dimension must be 2 !")
            self.derivation_type = derivation_type
            if self.derivation_type == DerivationType.SYMMETRIC:
                self.field_dimension = 2
                self.gradient_dimension = 4
                self.voigt_indices = {
                    (0, 0): 0,
                    (1, 1): 1,
                    (0, 1): 3
                }
                self.voigt_coefficients = {
                    (0, 0): 1.0,
                    (1, 1): 1.0,
                    (0, 1): np.sqrt(2.0)
                }
            elif self.derivation_type == DerivationType.FULL:
                self.field_dimension = 2
                self.gradient_dimension = 5
                self.voigt_indices = {
                    (0, 0): 0,
                    (1, 1): 1,
                    (0, 1): 3,
                    (1, 0): 4
                }
                self.voigt_coefficients = {
                    (0, 0): 1.0,
                    (1, 1): 1.0,
                    (0, 1): 1.0,
                    (1, 0): 1.0,
                }




#     self.derivation_type = derivation_type
#     if self.derivation_type == DerivationType.FULL:
#         self.field_dimension = euclidean_dimension
#         self.gradient_dimension = euclidean_dimension * euclidean_dimension + 1
#         if euclidean_dimension == 2:
#             self.voigt_indices = {
#                 (0, 0): 0,
#                 (1, 1): 1,
#                 (0, 1): 2,
#                 (1, 0): 3
#             }
#             self.voigt_coefficients = {
#                 (0, 0): 1.0,
#                 (1, 1): 1.0,
#                 (0, 1): 1.0,
#                 (1, 0): 1.0,
#             }
#         elif euclidean_dimension == 3:
#             self.voigt_indices = {
#                 (0, 0): 0,
#                 (1, 1): 1,
#                 (2, 2): 2,
#                 (0, 1): 3,
#                 (1, 0): 4,
#                 (0, 2): 5,
#                 (2, 0): 6,
#                 (1, 2): 7,
#                 (2, 1): 8,
#             }
#             self.voigt_coefficients = {
#                 (0, 0): 1.0,
#                 (1, 1): 1.0,
#                 (2, 2): 1.0,
#                 (0, 1): 1.0,
#                 (1, 0): 1.0,
#                 (0, 2): 1.0,
#                 (2, 0): 1.0,
#                 (1, 2): 1.0,
#                 (2, 1): 1.0,
#             }
#     elif self.derivation_type == DerivationType.SYMMETRIC:
#         self.field_dimension = euclidean_dimension
#         self.gradient_dimension = euclidean_dimension * euclidean_dimension
#         if euclidean_dimension == 2:
#             self.voigt_indices = {
#                 (0, 0): 0,
#                 (1, 1): 1,
#                 (0, 1): 2
#             }
#             self.voigt_coefficients = {
#                 (0, 0): 1.0,
#                 (1, 1): 1.0,
#                 (0, 1): np.sqrt(2.0)
#             }
#         elif euclidean_dimension == 3:
#             self.voigt_indices = {
#                 (0, 0): 0,
#                 (1, 1): 1,
#                 (2, 2): 2,
#                 (0, 1): 3,
#                 (0, 2): 4,
#                 (1, 2): 5,
#             }
#             self.voigt_coefficients = {
#                 (0, 0): 1.0,
#                 (1, 1): 1.0,
#                 (2, 2): 1.0,
#                 (0, 1): np.sqrt(2.0),
#                 (0, 2): np.sqrt(2.0),
#                 (1, 2): np.sqrt(2.0),
#             }
# elif self.field_type == FieldType.DISPLACEMENT_FINITE_STRAIN:
#     self.derivation_type = DerivationType.FULL
#     self.field_dimension = euclidean_dimension
#     self.gradient_dimension = euclidean_dimension * euclidean_dimension
#     if euclidean_dimension == 2:
#         self.voigt_indices = {
#             (0, 0): 0,
#             (1, 1): 1,
#             (0, 1): 2,
#             (1, 0): 3
#         }
#         self.voigt_coefficients = {
#             (0, 0): 1.0,
#             (1, 1): 1.0,
#             (0, 1): 1.0,
#             (1, 0): 1.0,
#         }
#     elif euclidean_dimension == 3:
#         self.voigt_indices = {
#             (0, 0): 0,
#             (1, 1): 1,
#             (2, 2): 2,
#             (0, 1): 3,
#             (1, 0): 4,
#             (0, 2): 5,
#             (2, 0): 6,
#             (1, 2): 7,
#             (2, 1): 8,
#         }
#         self.voigt_coefficients = {
#             (0, 0): 1.0,
#             (1, 1): 1.0,
#             (2, 2): 1.0,
#             (0, 1): 1.0,
#             (1, 0): 1.0,
#             (0, 2): 1.0,
#             (2, 0): 1.0,
#             (1, 2): 1.0,
#             (2, 1): 1.0,
#         }
# elif self.field_type == FieldType.DISPLACEMENT_SMALL_STRAIN:
#     self.derivation_type = DerivationType.SYMMETRIC
#     self.field_dimension = euclidean_dimension
#     self.gradient_dimension = euclidean_dimension * euclidean_dimension
#     if euclidean_dimension == 2:
#         self.voigt_indices = {
#             (0, 0): 0,
#             (1, 1): 1,
#             (0, 1): 2
#         }
#         self.voigt_coefficients = {
#             (0, 0): 1.0,
#             (1, 1): 1.0,
#             (0, 1): np.sqrt(2.0)
#         }
#     elif euclidean_dimension == 3:
#         self.voigt_indices = {
#             (0, 0): 0,
#             (1, 1): 1,
#             (2, 2): 2,
#             (0, 1): 3,
#             (0, 2): 4,
#             (1, 2): 5,
#         }
#         self.voigt_coefficients = {
#             (0, 0): 1.0,
#             (1, 1): 1.0,
#             (2, 2): 1.0,
#             (0, 1): np.sqrt(2.0),
#             (0, 2): np.sqrt(2.0),
#             (1, 2): np.sqrt(2.0),
#         }
# elif (
#     self.field_type == FieldType.DISPLACEMENT_SMALL_STRAIN_PLANE_STRAIN
#     or self.field_type == FieldType.DISPLACEMENT_SMALL_STRAIN_PLANE_STRESS
# ):
#     self.derivation_type = DerivationType.SYMMETRIC
#     if self.euclidean_dimension != 2:
#         raise ValueError("The euclidean dimension must be 2 !")
#     self.field_dimension = 2
#     self.gradient_dimension = 4
#     self.voigt_indices = {
#         (0, 0): 0,
#         (1, 1): 1,
#         (0, 1): 3
#     }
#     self.voigt_coefficients = {
#         (0, 0): 1.0,
#         (1, 1): 1.0,
#         (0, 1): np.sqrt(2.0)
#     }
# elif (
#     self.field_type == FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRAIN
#     or self.field_type == FieldType.DISPLACEMENT_FINITE_STRAIN_PLANE_STRESS
# ):
#     self.derivation_type = DerivationType.FULL
#     if self.euclidean_dimension != 2:
#         raise ValueError("The euclidean dimension must be 2 !")
#     self.field_dimension = 2
#     self.gradient_dimension = 5
#     self.voigt_indices = {
#         (0, 0): 0,
#         (1, 1): 1,
#         (0, 1): 3,
#         (1, 0): 4
#     }
#     self.voigt_coefficients = {
#         (0, 0): 1.0,
#         (1, 1): 1.0,
#         (0, 1): 1.0,
#         (1, 0): 1.0,
#     }
# self.derivation_type = derivation_type