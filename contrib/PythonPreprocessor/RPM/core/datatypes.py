from typing import List, Optional, Union
import sys
import os
from pydantic import model_validator

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
    position_wrt_base: Optional[l.Position] = None
    orientation_wrt_base: Optional[l.Position] = None
    velocity_wrt_base: Optional[l.Position] = None
    angular_velocity_wrt_base: Optional[l.Position] = None
    mirror: bool = False

    @model_validator(mode='after')
    def set_default_positions(self):
        """Set default Positions using base_reference after construction."""
        if self.position_wrt_base is None:
            self.position_wrt_base = l.Position(
                relative_position=l.null(), 
                reference=self.base_reference
            )
        if self.orientation_wrt_base is None:
            self.orientation_wrt_base = l.Position(
                relative_position=l.eye(), 
                reference=self.base_reference
            )
        if self.velocity_wrt_base is None:
            self.velocity_wrt_base = l.Position(
                relative_position=l.null(), 
                reference=self.base_reference
            )
        if self.angular_velocity_wrt_base is None:
            self.angular_velocity_wrt_base = l.Position(
                relative_position=l.null(), 
                reference=self.base_reference
            )
        return self

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
        