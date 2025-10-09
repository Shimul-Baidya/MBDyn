from typing import List
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from .base import RotorcraftComponent
from ..core.datatypes import ReferenceSystem

class Airframe(RotorcraftComponent):
    """Represents the main body of the aircraft."""
    reference_system: ReferenceSystem

    def _create_references(self) -> List[l.Reference]:
        ref = l.Reference(idx=20000, # Using the hardcoded value for now
                         position=self.reference_system.position_wrt_base, 
                         orientation=self.reference_system.orientation_wrt_base, 
                         velocity=self.reference_system.velocity_wrt_base, 
                         angular_velocity=self.reference_system.angular_velocity_wrt_base)
        return [ref]

    def _create_nodes(self) -> List[l.Node]:
        node = l.DynamicNode(
            idx=20000, # Using the hardcoded value for now
            position=l.Position(relative_position=[0., 0., 0.], reference=self.reference_system.base_reference),
            orientation=l.Position(relative_position=[l.eye()], reference=self.reference_system.base_reference),
            velocity=l.Position(relative_position=[l.null()], reference=self.reference_system.base_reference),
            angular_velocity=l.Position(relative_position=[l.null()], reference=self.reference_system.base_reference)
        )
        return [node]

    def _create_elements(self) -> List[l.Element]:
        return []
