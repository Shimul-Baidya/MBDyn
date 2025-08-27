from abc import ABC, abstractmethod
from typing import List, Optional, Union
from MBDynLib import *

imported_pydantic = False
try:
    from pydantic import BaseModel, ConfigDict, field_validator, FieldValidationInfo, model_validator
    imported_pydantic = True
    class _EntityBase(BaseModel):
        """Configuration for Entity with pydantic available"""
        model_config = ConfigDict(extra='forbid',
                                  use_attribute_docstrings=True)

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
    validate_call = identity_decorator


class MBEntity(_EntityBase, ABC):
    """Base class for every 'thing' to put in MBDyn file, other than numbers"""

    @abstractmethod
    def __str__(self) -> str:
        """Has to be overridden to output the MBDyn syntax"""
        pass

# --- Helper classes ---

class PhysicalQuantity(MBEntity):
    """Represents a physical quantity with a value and a unit."""
    unit: str
    value: Union[float, MBVar, List[Union[float, MBVar]]]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(value={self.value}, unit='{self.unit}')"

    def to_SI(self):
        # TODO: implement conversion to SI units
        pass

class ReferenceSystem(MBEntity):
    """Defines a coordinate system in the MBDyn model."""
    label: str
    component_axis: Optional[str] = 'x'
    base_reference: str
    position_wrt_base: Position
    orientation_wrt_base: Position
    velocity_wrt_base: Position
    angular_velocity_wrt_base: Position
    mirror: Optional[bool] = False

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(label='{self.label}', base_reference='{self.base_reference}')"

class LumpedMass(MBEntity):
    """Represents a lumped mass with inertial properties."""
    label: str
    reference: Optional[str] = None  # if None, default to the component's main reference frame
                                    # logic will be in the parent RotorcraftComponent that owns the mass
    mass: PhysicalQuantity
    relative_center_of_mass: PhysicalQuantity
    diag_inertia_matrix: PhysicalQuantity

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(label='{self.label}', mass={self.mass.value})"

# --- Mixin for collecting MBDyn entities ---

class ComponentCollectorMixin:
    """
    Mixin class that provides recursive methods to collect all MBDyn
    entities from a tree of components.
    """
    @property
    def _components_to_collect(self) -> List['RotorcraftComponent']:
        """Abstract property that child classes must override to specify the list of components to iterate over."""
        raise NotImplementedError

    def get_all_references(self) -> List[Reference]:
        """Gathers all MBDyn reference nodes from this component and its sub-components."""
        all_refs = self._create_references()
        for component in self._components_to_collect:
            all_refs.extend(component.get_all_references())
        return all_refs

    def get_all_nodes(self) -> List[Node]:
        """Gathers all MBDyn structural nodes from this component and its sub-components."""
        all_nodes = self._create_nodes()
        for component in self._components_to_collect:
            all_nodes.extend(component.get_all_nodes())
        return all_nodes

    def get_all_elements(self) -> List[Element]:
        """Gathers all MBDyn elements from this component and its sub-components."""
        all_elements = self._create_elements()
        for component in self._components_to_collect:
            all_elements.extend(component.get_all_elements())
        return all_elements

# --- Main Base Classes ---

class RotorcraftComponent(MBEntity, ComponentCollectorMixin, ABC):
    """Abstract base class for any physical component of the rotorcraft"""
    reference_system: ReferenceSystem
    lumped_masses: Optional[List[LumpedMass]] = []
    sub_components: List['RotorcraftComponent'] = []

    @property
    def _components_to_collect(self) -> List['RotorcraftComponent']:
        return self.sub_components

    def add_sub_component(self, component: 'RotorcraftComponent'):
        self.sub_components.append(component)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(reference_system_label='{self.reference_system.label}', sub_components={len(self.sub_components)})"

    @abstractmethod
    def _create_references(self) -> List[Reference]:
        pass

    @abstractmethod
    def _create_nodes(self) -> List[Node]:
        pass

    @abstractmethod
    def _create_elements(self) -> List[Element]:
        pass

RotorcraftComponent.model_rebuild() # Resolve self reference
