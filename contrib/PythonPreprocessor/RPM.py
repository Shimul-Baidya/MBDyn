from abc import ABC, abstractmethod
from tkinter import NO
from turtle import position
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
    mirror: Optional[bool] = NO

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


class RotorcraftComponent(MBEntity):
    """The base class for all rotorcraft components.

    This class provides the fundamental attributes and methods for creating a hierarchical
    rotorcraft model. It is designed to be subclassed by specific components like
    Rotor, Blade, Hub, etc.
    """

    label: str
    reference: Union[Reference, 'RotorcraftComponent', str]
    position: Position
    orientation: Position
    children: List['RotorcraftComponent'] = []

    @abstractmethod
    def _mbdyn_str(self) -> str:
        """
        Generates the MBDyn input string for this component ONLY.
        This method is intended to be overridden by concrete subclasses.
        """
        raise NotImplementedError(
            f"The component '{self.label}' of type '{self.__class__.__name__}' must implement the _mbdyn_str method."
        )

    def __str__(self) -> str:
        """Recursively generates the MBDyn input string for this component and all its children.

        It orchestrates the generation of the MBDyn output by combining the 
        component's own string with the strings of all its children.
        """
        parent_str = self._mbdyn_str()
        children_str = "\n".join(str(child) for child in self.children)
        return "\n".join(filter(None, [parent_str, children_str]))

    # --- Hierarchy Management Methods ---

    def add_child(self, component: 'RotorcraftComponent') -> 'RotorcraftComponent':
        """Adds a child component and returns the child for method chaining."""
        self.children.append(component)
        return component

    def get_all_components(self) -> List['RotorcraftComponent']:
        """Returns a flat list of this component and all its descendants."""
        components = [self]
        for child in self.children:
            components.extend(child.get_all_components())
        return components

    def get_components_by_type(self, component_type: type) -> List['RotorcraftComponent']:
        """
        Finds all components of a specific type (e.g., Blade) in the hierarchy.
        """
        found_components = []
        # First, check if the current component itself is the type we're looking for.
        if isinstance(self, component_type):
            found_components.append(self)
        # Then, recursively ask all children to do the same search and add their findings.
        for child in self.children:
            found_components.extend(child.get_components_by_type(component_type))
        return found_components
