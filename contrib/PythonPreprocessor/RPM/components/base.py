from abc import ABC, abstractmethod
from typing import List
import sys
import os

# Adjust path for new location
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

class RPMEntity(BaseModel):
    """Base class for all RPM entities, providing Pydantic validation."""
    model_config = ConfigDict(extra='forbid', use_attribute_docstrings=True)

def _collect_entities_from_components(components: List['RotorcraftComponent']) -> tuple[List[l.Reference], List[l.Node], List[l.Element]]:
    """Helper function to collect entities from a list of components."""
    all_refs, all_nodes, all_elements = [], [], []
    for component in components:
        refs, nodes, elements = component.collect_entities()
        all_refs.extend(refs)
        all_nodes.extend(nodes)
        all_elements.extend(elements)
    return all_refs, all_nodes, all_elements

class RotorcraftComponent(RPMEntity, ABC):
    """Abstract base class for all rotorcraft components (e.g., Airframe, Rotor, Blade)."""
    # reference_system will be defined in subclasses that need it
    
    # Pydantic safely handles mutable defaults with = []
    sub_components: List['RotorcraftComponent'] = []
    _references: List[l.Reference] = []
    _nodes: List[l.Node] = []
    _elements: List[l.Element] = []

    def add_sub_component(self, component: 'RotorcraftComponent'):
        self.sub_components.append(component)

    def generate_mbdyn_entities(self, processor: 'ModelProcessor'):
        """
        Orchestrates the creation of all MBDyn entities for this component.
        This method is called AFTER labels have been processed.
        
        Args:
            processor: The ModelProcessor instance that assigned labels
        """
        self._references = self._create_references(processor)
        self._nodes = self._create_nodes(processor)
        self._elements = self._create_elements(processor)

        # Recursively generate for sub-components
        for sub in self.sub_components:
            sub.generate_mbdyn_entities(processor)

    @abstractmethod
    def _create_references(self, processor: 'ModelProcessor') -> List[l.Reference]:
        """Creates MBDyn Reference objects for this component."""
        pass

    @abstractmethod
    def _create_nodes(self, processor: 'ModelProcessor') -> List[l.Node]:
        """Creates MBDyn Node objects for this component."""
        pass

    @abstractmethod
    def _create_elements(self, processor: 'ModelProcessor') -> List[l.Element]:
        """Creates MBDyn Element objects for this component."""
        pass

    def collect_entities(self) -> tuple[List[l.Reference], List[l.Node], List[l.Element]]:
        """Recursively collects all MBDyn entities from this component and its sub-components."""
        all_refs = list(self._references)
        all_nodes = list(self._nodes)
        all_elements = list(self._elements)

        sub_refs, sub_nodes, sub_elements = _collect_entities_from_components(self.sub_components)
        all_refs.extend(sub_refs)
        all_nodes.extend(sub_nodes)
        all_elements.extend(sub_elements)
        
        return all_refs, all_nodes, all_elements

RotorcraftComponent.model_rebuild()

class Rotorcraft(RPMEntity):
    """Represents the entire rotorcraft, managing all major components."""
    name: str
    root_components: List[RotorcraftComponent] = []

    def add_root_component(self, component: RotorcraftComponent):
        self.root_components.append(component)

    def generate_all_mbdyn_entities(self, processor: 'ModelProcessor'):
        """Generates all MBDyn entities for the entire rotorcraft model."""
        for component in self.root_components:
            component.generate_mbdyn_entities(processor)

    def collect_all_entities(self) -> tuple[List[l.Reference], List[l.Node], List[l.Element]]:
        """Collects all MBDyn entities from all root components in the model."""
        return _collect_entities_from_components(self.root_components)

    def write_mbdyn_files(self, output_dir: str):
        """Writes all generated MBDyn entities to their respective files."""
        # This method would collect all _nodes, _elements, etc. from the entire
        # component tree and write them to files.
        pass


# Public API declaration
__all__ = [
    # Classes that other modules will use
    'RPMEntity',
    'RotorcraftComponent', 
    'Rotorcraft',
]