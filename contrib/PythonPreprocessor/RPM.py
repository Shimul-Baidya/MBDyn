from abc import ABC, abstractmethod
from typing import List, Optional, Union
import MBDynLib as l


imported_pydantic = False
try:
    from pydantic import BaseModel, ConfigDict, field_validator, model_validator
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


class RPMEntity(_EntityBase):
    """Base class for all RPM entities, providing Pydantic validation."""
    pass

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
    label: str
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
        if 'position_wrt_base' not in kwargs:
            kwargs['position_wrt_base'] = l.Position(relative_position=[0., 0., 0.], reference='')
        if 'orientation_wrt_base' not in kwargs:
            kwargs['orientation_wrt_base'] = l.Position(relative_position=[1., 0., 0., 0., 1., 0., 0., 0., 1.], reference='')
        if 'velocity_wrt_base' not in kwargs:
            kwargs['velocity_wrt_base'] = l.Position(relative_position=[0., 0., 0.], reference='')
        if 'angular_velocity_wrt_base' not in kwargs:
            kwargs['angular_velocity_wrt_base'] = l.Position(relative_position=[0., 0., 0.], reference='')
        super().__init__(**kwargs)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(label='{self.label}', base_reference='{self.base_reference}')"


class LumpedMass(RPMEntity):
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

    def get_all_references(self) -> List[l.Reference]:
        """Gathers all MBDyn reference nodes from this component and its sub-components."""
        all_refs = self._create_references()
        for component in self._components_to_collect:
            all_refs.extend(component.get_all_references())
        return all_refs

    def get_all_nodes(self) -> List[l.Node]:
        """Gathers all MBDyn structural nodes from this component and its sub-components."""
        all_nodes = self._create_nodes()
        for component in self._components_to_collect:
            all_nodes.extend(component.get_all_nodes())
        return all_nodes

    def get_all_elements(self) -> List[l.Element]:
        """Gathers all MBDyn elements from this component and its sub-components."""
        all_elements = self._create_elements()
        for component in self._components_to_collect:
            all_elements.extend(component.get_all_elements())
        return all_elements

# --- Main Base Classes ---

class RotorcraftComponent(RPMEntity, ComponentCollectorMixin, ABC):
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
    def _create_references(self) -> List[l.Reference]:
        pass

    @abstractmethod
    def _create_nodes(self) -> List[l.Node]:
        pass

    @abstractmethod
    def _create_elements(self) -> List[l.Element]:
        pass

RotorcraftComponent.model_rebuild() # Resolve self reference


class Rotorcraft(RPMEntity, ComponentCollectorMixin):
    """The main container for the root components and entire model"""
    name: str
    root_components: List[RotorcraftComponent] = []

    @property
    def _components_to_collect(self) -> List['RotorcraftComponent']:
        return self.root_components
    
    def add_root_component(self, component: RotorcraftComponent):
        self.root_components.append(component)

    # CONCRETE IMPLEMENTATIONS FOR THE TOP-LEVEL CONTAINER
    def _create_references(self) -> List[l.Reference]: return []
    def _create_nodes(self) -> List[l.Node]: return []
    def _create_elements(self) -> List[l.Element]: return []
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}', root_components={len(self.root_components)})"


# --- Component Classes ---

class Airframe(RotorcraftComponent):
    # The only attribute needed here to be input by the user is the 
    # reference_system, which will be inherited from the 'RotorcraftComponent' class

    def _create_references(self) -> List[l.Reference]:
        # The airframe's primary reference system is the only one it creates.
        ref = l.Reference(idx=20000, # self.reference_system.label # Using the hardcoded value for now
                        position=self.reference_system.position_wrt_base, 
                        orientation=self.reference_system.orientation_wrt_base, 
                        velocity=self.reference_system.velocity_wrt_base, 
                        angular_velocity=self.reference_system.angular_velocity_wrt_base)
        return [ref]

    def _create_nodes(self) -> List[l.Node]:
        # The airframe has a single dynamic node at its reference point.
        node = l.DynamicNode(
            idx=10000, # self.reference_system.label # Using the hardcoded value for now
            position=self.reference_system.position_wrt_base,
            orientation=self.reference_system.orientation_wrt_base,
            velocity=self.reference_system.velocity_wrt_base,
            angular_velocity=self.reference_system.angular_velocity_wrt_base
        )
        return [node]

    def _create_elements(self) -> List[l.Element]:
        # A simple airframe does not create any elements itself.
        return []