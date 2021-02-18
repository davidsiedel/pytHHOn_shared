from typing import Callable
from pythhon.parameters import *


class BoundaryCondition:
    boundary_name: str
    function: Callable
    boundary_type: BoundaryType
    direction: int

    def __init__(self, boundary_name: str, function: Callable, boundary_type: BoundaryType, direction: int):
        """

        Args:
            boundary_name:
            function:
            boundary_type:
            direction:
        """
        self.__check_function(function)
        self.boundary_name = boundary_name
        self.function = function
        self.boundary_type = boundary_type
        self.direction = direction
        return

    def __check_function(self, function: Callable):
        """

        Args:
            function:
        """
        if hasattr(function, "__call__"):
            try:
                function(0.0, np.array([0.0]))
            except:
                raise ValueError
