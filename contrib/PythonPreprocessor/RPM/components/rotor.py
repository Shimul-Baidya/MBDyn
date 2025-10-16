from typing import List, TYPE_CHECKING
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from .base import RotorcraftComponent
from ..core.datatypes import ReferenceSystem, PhysicalQuantity

if TYPE_CHECKING:
    from ..core.processor import ModelProcessor

class Rotor(RotorcraftComponent):
    """
    Represents a rotor assembly as a non-physical container.
    
    The Rotor class serves as a static (non-rotating) reference frame anchor point
    for the entire rotor system. It holds key design parameters (n_blades, radius, 
    precone) that are used by its child components (Hub, Mast, Blades, etc.).
    
    According to the design specification:
    - Generates 1 Reference object (static, non-rotating)
    - Generates 0 Node objects (returns empty list)
    - Generates 0 Element objects (returns empty list)
    
    Attributes:
        reference_system: The coordinate system definition for this rotor
        n_blades: Number of blades in the rotor
        radius: Full radius of the rotor
        precone: Built-in angle of the blades (coning angle)
        precone_start: Radial position where precone begins
    """
    reference_system: ReferenceSystem
    n_blades: int
    radius: PhysicalQuantity
    precone: PhysicalQuantity
    precone_start: PhysicalQuantity

    def _create_references(self, processor: 'ModelProcessor') -> List[l.Reference]:
        """
        Creates a single static (non-rotating) MBDyn Reference for the rotor.
        
        This reference frame is positioned relative to its parent (typically the Airframe)
        and serves as the base reference for all rotating components (Hub, Mast, etc.).
        
        Returns:
            List containing one Reference object
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
        Returns an empty list as Rotor is a non-physical container.
        
        The Rotor does not represent a physical component, so it has no nodes.
        Nodes are created by its child components (Hub, Mast, Blades, etc.).
        
        Returns:
            Empty list
        """
        return []

    def _create_elements(self, processor: 'ModelProcessor') -> List[l.Element]:
        """
        Returns an empty list as Rotor is a non-physical container.
        
        The Rotor does not represent a physical component, so it has no elements.
        Elements are created by its child components.
        
        Returns:
            Empty list
        """
        return []
