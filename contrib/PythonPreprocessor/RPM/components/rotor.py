from typing import List, TYPE_CHECKING
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from .base import RotorcraftComponent
from ..core.datatypes import ReferenceSystem

if TYPE_CHECKING:
    from ..core.processor import ModelProcessor

# This is not final, it's only used for testing ModelProcessor
class Rotor(RotorcraftComponent):
    """Represents a rotor assembly."""
    reference_system: ReferenceSystem
    n_blades: int

    def _create_references(self, processor: 'ModelProcessor') -> List[l.Reference]:
        ref = l.Reference(idx=self.reference_system.mbdyn_label,
                         position=self.reference_system.position_wrt_base, 
                         orientation=self.reference_system.orientation_wrt_base, 
                         velocity=self.reference_system.velocity_wrt_base, 
                         angular_velocity=self.reference_system.angular_velocity_wrt_base)
        return [ref]

    def _create_nodes(self, processor: 'ModelProcessor') -> List[l.Node]:
        node_idx = processor.create_internal_label(owner=self, entity_type='node', counter=1)
        node = l.DynamicNode(
            idx=node_idx,
            position=l.Position(relative_position=[0., 0., 0.], reference=self.reference_system.base_reference),
            orientation=l.Position(relative_position=[l.eye()], reference=self.reference_system.base_reference),
            velocity=l.Position(relative_position=[l.null()], reference=self.reference_system.base_reference),
            angular_velocity=l.Position(relative_position=[l.null()], reference=self.reference_system.base_reference)
        )
        return [node]

    def _create_elements(self, processor: 'ModelProcessor') -> List[l.Element]:
        return []
