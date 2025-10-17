from typing import List, TYPE_CHECKING
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from .base import RotorcraftComponent
from ..core.datatypes import ReferenceSystem

if TYPE_CHECKING:
    from ..core.processor import ModelProcessor

class Hub(RotorcraftComponent):
    """
    Represents the central rotating component of the rotor assembly.
    
    The Hub class provides the rotating reference frame to which the blades are attached.
    It connects to the static Rotor reference frame via a joint that is driven to rotate 
    at the specified rotor speed.
    
    According to the design specification:
    - Generates 1 Reference object (rotating frame)
    - Generates 1 Node object (DynamicNode at the geometric center)
    - Generates 0 Element objects by default (massless kinematic entity)
    - Mass can be added via lumped_masses attribute
    
    Attributes:
        reference_system: The coordinate system definition for this hub.
    
    Example:
        >>> hub = Hub(
        ...     reference_system=ReferenceSystem(
        ...         name='hub_r1',
        ...         base_reference='parent',
        ...         position_wrt_base=Position([0., 0., 0.], 'ROTOR_1'),
        ...         angular_velocity_wrt_base=Position([0., 0., 0.], 'ROTOR_1')
        ...     )
        ... )
    """
    reference_system: ReferenceSystem

    def _create_references(self, processor: 'ModelProcessor') -> List[l.Reference]:
        """
        Creates a single rotating MBDyn Reference for the hub.
        """
        ref = l.Reference(
            idx=self.reference_system.mbdyn_label,
            position=self.reference_system.position_wrt_base, 
            orientation=self.reference_system.orientation_wrt_base, 
            velocity=self.reference_system.velocity_wrt_base, 
            angular_velocity=self.reference_system.angular_velocity_wrt_base
        )
        return [ref]

    def _create_nodes(self, processor: 'ModelProcessor') -> List[l.Node]:
        """
        Creates a single MBDyn DynamicNode at the hub center.
        
        This node represents the geometric center of the hub and serves as the 
        attachment point for blade roots, yokes, and the mast.
        """
        # Create internal node label using processor
        # Counter is 1 for the single hub node (station 0)
        node_idx = processor.create_internal_label(
            owner=self, 
            entity_type='NODE', 
            counter=1
        )
        
        node = l.DynamicNode(
            idx=node_idx,
            position=l.Position(
                relative_position=[0., 0., 0.], 
                reference=self.reference_system.base_reference
            ),
            orientation=l.Position(
                relative_position=l.eye(), 
                reference=self.reference_system.base_reference
            ),
            velocity=l.Position(
                relative_position=l.null(), 
                reference=self.reference_system.base_reference
            ),
            angular_velocity=l.Position(
                relative_position=l.null(), 
                reference=self.reference_system.base_reference
            )
        )
        return [node]

    def _create_elements(self, processor: 'ModelProcessor') -> List[l.Element]:
        """Creates MBDyn elements for the hub (currently none by default)."""

        return []