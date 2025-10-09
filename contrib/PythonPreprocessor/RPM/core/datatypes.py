from typing import List, Optional, Union
import sys
import os

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from ..components.base import RPMEntity, imported_pydantic

# --- Helper classes ---

class PhysicalQuantity(RPMEntity):
    """Represents a physical quantity with a value and a unit."""
    unit: str
    value: Union[float, l.MBVar, List[Union[float, l.MBVar]]]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(value={self.value}, unit='{self.unit}')"

    def to_SI(self):
        # TODO: implement conversion to SI units
        pass

class ReferenceSystem(RPMEntity):
    """Defines a coordinate system in the MBDyn model."""
    name: str
    mbdyn_label: Optional[l.MBVar] = None
    component_axis: Optional[str] = 'x'
    base_reference: str
    position_wrt_base: l.Position
    orientation_wrt_base: l.Position
    velocity_wrt_base: l.Position
    angular_velocity_wrt_base: l.Position
    mirror: bool = False

    def __init__(self, **kwargs):
        # A custom __init__ is used to set default values for mutable types like l.Position.
        # If we set a mutable default directly on the class (e.g., `position_wrt_base: l.Position = l.Position(...)`),
        # all instances of ReferenceSystem would share the same Position object, leading to unintended side effects.
        # This method ensures that a new, unique Position object is created for each instance that doesn't provide one.
        # It works regardless of whether Pydantic is installed.
        kwargs.setdefault('position_wrt_base', l.Position(relative_position=[0., 0., 0.], reference=''))
        kwargs.setdefault('orientation_wrt_base', l.Position(relative_position=[1., 0., 0., 0., 1., 0., 0., 0., 1.], reference=''))
        kwargs.setdefault('velocity_wrt_base', l.Position(relative_position=[0., 0., 0.], reference=''))
        kwargs.setdefault('angular_velocity_wrt_base', l.Position(relative_position=[0., 0., 0.], reference=''))
        super().__init__(**kwargs)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}', base_reference='{self.base_reference}')"


class LumpedMass(RPMEntity):
    """Represents a lumped mass with inertial properties."""
    name: str
    mbdyn_label: Optional[l.MBVar] = None
    reference: Optional[str] = None
    mass: PhysicalQuantity
    inertia_tensor: PhysicalQuantity

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}', mass={self.mass.value})"
        