from abc import ABC, abstractmethod
from typing import List
import sys
import os

# Adjust path for new location
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l

# --- Optional Pydantic Boilerplate ---
imported_pydantic = False
try:
    from pydantic import BaseModel, ConfigDict, field_validator, model_validator
    imported_pydantic = True
    class _EntityBase(BaseModel):
        """Configuration for Entity with pydantic available"""
        model_config = ConfigDict(extra='forbid', use_attribute_docstrings=True)

except ImportError:
    class _EntityBasePlaceholder:
        """Placeholder with minimal functionality for running a correct model when some libraries aren't available"""

        def __init__(self, *args, **kwargs):
            if len(args) > 0:
                raise TypeError(
                    'MBDyn entities cannot be initialized using positional arguments')
            for key, value in kwargs.items():
                setattr(self, key, value)

    def placeholder(*args, **kwargs):
        """Ignores all arguments"""
        return None

    # HACK: This forces code analysis to always use the definition with pydantic
    exec('_EntityBase = _EntityBasePlaceholder')
    exec('ConfigDict = placeholder')

    def identity_decorator(*args, **kwargs):
        """Ignores all decorator arguments and returns the wrapped function unchanged"""
        def identity(func):
            return func
        return identity

    field_validator = identity_decorator
    model_validator = identity_decorator
# --- End Boilerplate ---

class RPMEntity(_EntityBase):
    """Base class for all RPM entities, providing optional Pydantic validation."""
    pass

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
    
    # Remove the default values here - they'll be set in __init__
    sub_components: List['RotorcraftComponent'] = []
    _references: List[l.Reference] = []
    _nodes: List[l.Node] = []
    _elements: List[l.Element] = []

    def add_sub_component(self, component: 'RotorcraftComponent'):
        self.sub_components.append(component)

    def generate_mbdyn_entities(self):
        """
        Orchestrates the creation of all MBDyn entities for this component.
        This method is called AFTER labels have been processed.
        """
        # This check will be specific to components that have a reference system
        # if self.reference_system.mbdyn_label is None:
        #     raise RuntimeError(f"Component '{self.reference_system.name}' has not been processed.")

        self._references = self._create_references()
        self._nodes = self._create_nodes()
        self._elements = self._create_elements()

        # Recursively generate for sub-components
        for sub in self.sub_components:
            sub.generate_mbdyn_entities()

    @abstractmethod
    def _create_references(self) -> List[l.Reference]:
        """Creates MBDyn Reference objects for this component."""
        pass

    @abstractmethod
    def _create_nodes(self) -> List[l.Node]:
        """Creates MBDyn Node objects for this component."""
        pass

    @abstractmethod
    def _create_elements(self) -> List[l.Element]:
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

if imported_pydantic:
    RotorcraftComponent.model_rebuild()

class Rotorcraft(RPMEntity):
    """Represents the entire rotorcraft, managing all major components."""
    name: str
    root_components: List[RotorcraftComponent] = []

    def add_root_component(self, component: RotorcraftComponent):
        self.root_components.append(component)

    def generate_all_mbdyn_entities(self):
        """Generates all MBDyn entities for the entire rotorcraft model."""
        for component in self.root_components:
            component.generate_mbdyn_entities()

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
    
    # Decorators/functions that other modules might need
    'field_validator',
    'model_validator',
    'ConfigDict',
    
    # Flags that other modules might check
    'imported_pydantic'
]