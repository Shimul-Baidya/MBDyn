from typing import List, Optional, Union
import sys
import os
from pydantic import Field

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from ..components.base import RPMEntity

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
    position_wrt_base: l.Position = Field(
        default_factory=lambda: l.Position(relative_position=[l.null()], reference='')
    )
    orientation_wrt_base: l.Position = Field(
        default_factory=lambda: l.Position(relative_position=[l.eye()], reference='')
    )
    velocity_wrt_base: l.Position = Field(
        default_factory=lambda: l.Position(relative_position=[l.null()], reference='')
    )
    angular_velocity_wrt_base: l.Position = Field(
        default_factory=lambda: l.Position(relative_position=[l.null()], reference='')
    )
    mirror: bool = False

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
        