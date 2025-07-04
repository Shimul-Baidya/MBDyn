#
#MBDyn (C) is a multibody analysis code.
#http://www.mbdyn.org
#
#Copyright (C) 1996-2023
#
#Pierangelo Masarati	<pierangelo.masarati@polimi.it>
#Paolo Mantegazza	<paolo.mantegazza@polimi.it>
#
#Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
#via La Masa, 34 - 20156 Milano, Italy
#http://www.aero.polimi.it
#
#Changing this copyright notice is forbidden.
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation (version 2 of the License).
#
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#COPYRIGHT (C) 2016
#
#Marco Morandini <marco.morandini@polimi.it>
#Mattia Alioli   <mattia.alioli@polimi.it>
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.


from abc import ABC, abstractmethod
import builtins
from enum import Enum
from numbers import Number, Integral
import sys
from typing import Optional, Tuple, Union, List, Literal, Any, ClassVar, Tuple
from typing_extensions import Self
import warnings


assert sys.version_info >= (3,6), 'Syntax for variable annotations (PEP 526) was introduced in Python 3.6'

declared_ConstMBVars = {}
declared_IfndefMBVars = {}
declared_MBVars = {}

MBDynLib_simplify = True

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


def errprint(*args, **kwargs):
    print(*args, file = sys.stderr, **kwargs)

def get_value(x):
    if isinstance(x, expression):
        return x.__get__()
    else:
        return x

def simplify_null_element_multiplication(l, r):
    if MBDynLib_simplify:
        if l == 0 or r == 0:
            return True
        else:
            return False
    else:
        return False

def simplify_null_element_division(l, r):
    assert get_value(r) != 0, (
        'Error, division by zero: \'' + str(l) + ' / ' + str(r) + 
        '\'\n')
    if MBDynLib_simplify:
        if l == 0:
            return True
        else:
            return False
    else:
        return False

def simplify_neutral_element(l, r, op, ne):
    if MBDynLib_simplify:
        #if get_value(l) == ne:
        if l == ne:
            return r
        #elif get_value(r) == ne:
        elif r == ne:
            return l
        else:
            return op(l, r)
    else:
        return op(l, r)

class expression:
    def __init__(self):
        pass
    def __float__(self):
        return float(self.__get__())
    def __int__(self):
        return int(self.__get__())
    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return self.__get__() == other
        return NotImplemented
    def __neg__(self):
        return negative(self)
    def __add__(self, other):
            return simplify_neutral_element(self, other, addition, 0) #addition(self, other)
    def __sub__(self, other):
            return simplify_neutral_element(self, other, subtraction, 0) #subtraction(self, other)
    def __pow__(self, other):
            return power(self, other)
    def __mul__(self, other):
            if simplify_null_element_multiplication(self, other):
                return 0
            else:
                return simplify_neutral_element(self, other, multiplication, 1) #multiplication(self, other)
    def __truediv__(self, other):
            if simplify_null_element_division(self, other):
                return 0
            else:
                return division(self, other)
    def __radd__(self, other):
            return simplify_neutral_element(other, self, addition, 0) #addition(other, self)
    def __rsub__(self, other):
            return simplify_neutral_element(other, self, subtraction, 0) #subtraction(other, self)
    def __rmul__(self, other):
            if simplify_null_element_multiplication(other, self):
                return 0
            else:
                return simplify_neutral_element(other, self, multiplication, 1) #multiplication(other, self)
    def __rtruediv__(self, other):
            if simplify_null_element_division(other, self):
                return 0
            else:
                return division(other, self)

class negative(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return -get_value(self.left)
    def __str__(self):
        ls = str(self.left)
        if isinstance(self.left, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.right, expression):
            ls = '(' + str(self.left) +')'
        return '-' + ls

import math
class sin(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return math.sin(get_value(self.left))
    def __str__(self):
        ls = str(self.left)
        return 'sin(' + ls + ')'

class cos(expression):
    def __init__(self, left):
        expression.__init__(self)
        self.left = left
    def __get__(self):
        return math.cos(get_value(self.left))
    def __str__(self):
        ls = str(self.left)
        return 'cos(' + ls + ')'

class tan(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.tan(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'tan(' + ls + ')'

class asin(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.asin(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'asin(' + ls + ')'

class acos(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.acos(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'acos(' + ls + ')'

class sqrt(expression):
	def __init__(self, left):
		expression.__init__(self)
		self.left = left
	def __get__(self):
		return math.sqrt(get_value(self.left))
	def __str__(self):
		ls = str(self.left)
		return 'sqrt(' + ls + ')'

class terminal_expression(expression):
    def __init__(self, value):
        expression.__init__(self)
        self.value = value
    def __get___(self):
        return self.value
    def __str__(self):
        return str(self.value)

class binary_expression(expression):
    def __init__(self, left, right):
        expression.__init__(self)
        self.left = left
        self.right = right
    def __trunc__(self):
        y = self.__get__()
        assert isinstance(y, int), (
                'Error, __trunc__  required for expression \n\'' + 
                str(self) + 
                '\'\nof type ' + str(type(y)) +
                ' \n')
        return y
    def __index__(self):
        return self.__trunc__()

class atan2(binary_expression):
	def __init__(self, left, right):
		binary_expression.__init__(self, left, right)
	def __get__(self):
		return math.atan2(get_value(self.left), get_value(self.right))
	def __str__(self):
		ls = str(self.left)
		rs = str(self.right)
		return 'atan2(' + ls + ', ' + rs + ')'

class addition(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) + get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        return ls + ' + ' + rs
            
class subtraction(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) - get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        return ls + ' - ' + rs
            
class multiplication(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) * get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, addition) or isinstance(self.left, subtraction):
            ls = '(' + ls + ')'
        if isinstance(self.right, addition) or isinstance(self.right, subtraction):
            rs = '(' + rs + ')'
        return ls + ' * ' + rs
            
class division(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return get_value(self.left) / get_value(self.right)
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, addition) or isinstance(self.left, subtraction):
            ls = '(' + ls + ')'
        if isinstance(self.right, terminal_expression) or isinstance(self.right, MBVar) or isinstance(self.right, power):
            pass
        elif isinstance(self.right, expression):
            rs = '(' + rs + ')'
        return ls + ' / ' + rs

class power(binary_expression):
    def __init__(self, left, right):
        binary_expression.__init__(self, left, right)
    def __get__(self):
        return pow(get_value(self.left), get_value(self.right))
    def __str__(self):
        ls = str(self.left)
        rs = str(self.right)
        if isinstance(self.left, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.left, expression):
            ls = '(' + ls + ')'
        if isinstance(self.right, terminal_expression) or isinstance(self.right, MBVar):
            pass
        elif isinstance(self.right, expression):
            rs = '(' + rs + ')'
        return ls + ' ^ ' + rs


# Enum entries need annotated type to find descriptions for documentation
class MBVarType(str, Enum):
    """Built-in types in math parser"""
    BOOL: str = 'bool'
    INTEGER: str = 'integer'
    REAL: str = 'real'
    STRING: str = 'string'

class MBVarModifiers(str, Enum):
    CONST: str = 'const'
    DEFINE: str = 'ifndef const'

class MBVar(MBEntity, terminal_expression):
    name: str
    var_type: str
    expression: Any

    var_types: ClassVar[Tuple[str]] = tuple([t.value for t in MBVarType] +\
                                  [f'{MBVarModifiers.CONST} {t.value}' for t in MBVarType] +\
                                  [f'{MBVarModifiers.DEFINE} {t.value}' for t in MBVarType])

    @field_validator('var_type')
    def validate_var_type(cls, v):
        assert v.strip('const ') in cls.var_types, (
            f'\n-------------------\nERROR: MBVar: unknown variable type {v}\n\t' +
            '\n-------------------\n'
        )
        return v

    @model_validator(mode='after')
    def validate_declarations(self):
        assert self.name not in declared_ConstMBVars, (
            '\n-------------------\nERROR: re-defining an already declared const variable:\n\t' +
            f'{self.var_type} {self.name}\n-------------------\n'
        )
        
        assert self.name not in declared_IfndefMBVars, (
            '\n-------------------\nERROR: re-defining an already declared ifndef variable:\n\t' +
            f'{self.var_type} {self.name}\n-------------------\n'
        )
        return self

    def __init__(self, name: str, var_type: str, expression: Any):
        super().__init__(name=name, var_type=var_type, expression=expression)
        self.declare()

    def __get__(self):
        return get_value(self.expression)

    def __trunc__(self):
        y = self.__get__()
        assert isinstance(y, int), (
            f'Error, __trunc__ required for expression \n\'{self}\'\nof type {type(y)}\n'
        )
        return self.expression.__trunc__()

    def __index__(self):
        return self.expression.__trunc__()

    def __str__(self):
        return str(self.name)
    
    def __repr__(self):
        # This ensures that when the object is used in lists or other contexts,
        # it still outputs just the name
        return self.__str__()

    def __lt__(self, other):
        return self.__get__() < other

    def __gt__(self, other):
        return self.__get__() > other

    def __eq__(self, other):
        return self.__get__() == other

    def __le__(self, other):
        return self.__get__() <= other

    def __ge__(self, other):
        return self.__get__() >= other

    def declare(self):
        if self.name in declared_MBVars:
            assert declared_MBVars[self.name].var_type == self.var_type, (
                '\n-------------------\nERROR: re-defining an already declared variable of type ' +
                f'{declared_MBVars[self.name].var_type}\nwith different type {self.var_type}\n' +
                '\n-------------------\n'
            )
            if 'string' in self.var_type:
                print(f'set: {self.name} = "{str(self.expression)}";')
            else:
                print(f'set: {self.name} = {str(self.expression)};')
        else:
            declared_MBVars[self.name] = self
            if 'string' in self.var_type:
                print(f'set: {self.var_type} {self.name} = "{str(self.expression)}";')
            else:
                print(f'set: {self.var_type} {self.name} = {str(self.expression)};')
        setattr(builtins, self.name, self)

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

class ConstMBVar(MBVar):
    def __init__(self, name: str, var_type: str, value: Any):
        super().__init__(name=name, var_type=f'const {var_type}', expression=value)

    def declare(self):
        super().declare()
        declared_ConstMBVars[self.name] = self

class IfndefMBVar(MBVar):
    def __init__(self, name: str, var_type: str, value: Any):
        if name not in declared_MBVars:
            super().__init__(name=name, var_type=f'ifndef {var_type}', expression=value)


class null(MBEntity):
    def __str__(self):
        return 'null'

class eye(MBEntity):
    def __str__(self):
        return 'eye'

class Position(MBEntity):
    """Position definition for MBDyn elements"""
        
    relative_position: List[Union[float, MBVar, null, eye]]
    reference: Union['Reference', Literal['global', 'node', 'other node', '']] # TODO: Make reference an optional field (remove '')

    @field_validator('relative_position', mode='before')
    def ensure_list(cls, v):
        if not isinstance(v, list):
            return [v]
        return v

    @field_validator('reference')
    def validate_reference(cls, v):
        if isinstance(v, str):
            if v not in {'global', 'node', 'other node', ''}:
                raise ValueError("Invalid literal for reference")
        elif not isinstance(v, Reference):
            raise ValueError("reference must be either a Reference instance or one of the specified strings")
        return v

    def __str__(self):
        s = ''
        if self.reference != '':
            s = 'reference, ' + str(self.reference) + ', '
        s = s + ', '.join(str(i) for i in self.relative_position)
        return s

    def isnull(self) -> bool:
        return (self.reference == '') and isinstance(self.relative_position[0], null)

    def iseye(self) -> bool:
        return (self.reference == '') and isinstance(self.relative_position[0], eye)
    
class Reference(MBEntity):
    idx: Union[int, MBVar]
    position: Position
    orientation: Position
    velocity: Position
    angular_velocity: Position    
    def __str__(self):
        s = 'reference: '
        s = s + str(self.idx) + ', \n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity) + ';\n'
        return s

if imported_pydantic:
    Position.model_rebuild()

class Node(MBEntity):
    """This class isn't directly used to create instances, but it's child classes are."""

    idx: Union[int, MBVar]
    position: Position
    orientation: Position
    velocity: Position
    angular_velocity: Position
    node_type: Literal['dynamic', 'static', 'modal'] = 'dynamic'
    scale: Optional[Union[Literal['default'], float, MBVar]] = 'default'
    output: Optional[Union[Literal['yes', 'no'], int, bool]] = 'yes'
    def __str__(self):
        s = f"structural: {self.idx}, {self.node_type},\n"
        s += f"\t{self.position},\n"
        s += f"\t{self.orientation},\n"
        s += f"\t{self.velocity},\n"
        s += f"\t{self.angular_velocity}"
        if self.scale != 'default':
            s += f",\n\tscale, {self.scale}"
        if self.output != 'yes':
            s += f",\n\toutput, {self.output}"
        return s
    
class DynamicNode(Node):
    accelerations: Optional[Union[Literal['yes', 'no'], int, bool]] = None
    def __init__(self, idx, pos, orient, vel, angular_vel, accelerations=None):
        super().__init__(idx=idx, position=pos, orientation=orient, velocity=vel, angular_velocity=angular_vel, node_type='dynamic')
        self.accelerations = accelerations
    def __str__(self):
        s = super().__str__()
        if self.accelerations is not None:
            s += f",\n\taccelerations, {self.accelerations}"
        s += ';\n'
        return s

class StaticNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        super().__init__(idx=idx, position=pos, orientation=orient, velocity=vel, angular_velocity=angular_vel, node_type='static')
    def __str__(self):
        return super().__str__() + ';\n'

class ModalNode(Node):
    """
    The modal node is basically a regular dynamic node that must be used to describe the rigid reference
    motion of a modal joint.
    """

    accelerations: Optional[Union[Literal['yes', 'no'], int, bool]] = None
    def __init__(self, idx, pos, orient, vel, angular_vel, accelerations=None):
        super().__init__(idx=idx, position=pos, orientation=orient, velocity=vel, angular_velocity=angular_vel, node_type='modal')
        self.accelerations = accelerations
    def __str__(self):
        s = super().__str__()
        if self.accelerations is not None:
            s += f",\n\taccelerations, {self.accelerations}"
        s += ';\n'
        return s

class DisplacementNode(MBEntity):
    idx: Union[int, MBVar]
    position: Position
    velocity: Position
    node_type: Literal['dynamic', 'static'] = 'dynamic'
    scale: Optional[Union[Literal['default'], float, MBVar]] = 'default'
    output: Optional[Union[Literal['yes', 'no'], int, bool]] = 'yes'
    def __str__(self):
        s = f"structural: {self.idx}, {self.node_type} displacement,\n"
        s += f"\t{self.position},\n"
        s += f"\t{self.velocity}"
        if self.scale != 'default':
            s += f",\n\tscale, {self.scale}"
        if self.output != 'yes':
            s += f",\n\toutput, {self.output}"
        return s

class DynamicDisplacementNode(DisplacementNode):
    accelerations: Optional[Union[Literal['yes', 'no'], int, bool]] = None
    def __init__(self, idx, pos, vel, accelerations=None):
        super().__init__(idx=idx, position=pos, velocity=vel, node_type='dynamic')
        self.accelerations = accelerations
    def __str__(self):
        s = super().__str__()
        if self.accelerations is not None:
            s += f",\n\taccelerations, {self.accelerations}"
        return s + ';\n'

class StaticDisplacementNode(DisplacementNode):
    def __init__(self, idx, pos, vel):
        super().__init__(idx=idx, position=pos, velocity=vel, node_type='static')
    def __str__(self):
        return super().__str__() + ';\n'

class PointMass(MBEntity):
    idx: Union[int, MBVar]
    node: Node
    mass: Union[float, MBVar]
    output: Optional[Union[Literal['yes', 'no'], int, bool]] = 'yes'
    
    def __str__(self):
        s = f"body: {self.idx}, {self.node}, {self.mass}"
        if self.output != 'yes':
            s += f", output, {self.output}"
        s += ";\n"
        return s
    
class Element(MBEntity):
    """
    Abstract base class for all elements
    """

    idx: Union[MBVar, int]
    output: Optional[Union[Literal['yes', 'no'], int, bool]] = 'yes'

    @abstractmethod
    def element_type(self) -> str:
        """Every element class must define this to return its MBDyn syntax name"""
        raise NotImplementedError("called elelemt_type of abstract Element")

    def element_header(self) -> str:
        """common syntax for start of any element"""
        return f'{self.element_type()}: {self.idx}'
    
    def element_footer(self) -> str:
        s = ''
        if self.output != 'yes':
            s = s + f''',\n\toutput, {self.output}'''
        s = s + ';\n'
        return s

    @staticmethod
    def check_unit_vector3(value: List[Union[float, MBVar]]):
        if not len(value) == 3:
            raise ValueError("relative_direction must be a 3-dimensional vector")

        magnitude = sum(v**2 for v in value)
        if not (0.999 <= magnitude <= 1.001):  # Allowing some tolerance for floating-point precision
            raise ValueError("relative_direction must be a unit vector (magnitude = 1)")

class Body(Element):
    node: Node
    mass: Union[float, MBVar]
    position: Position
    inertial_matrix: Position 
    inertial: Optional[Position] = None

    def element_type(self):
        return 'body'
    
    def __str__(self):
        s = f'{self.element_header()}, {self.node.idx}'
        s += f',\n\t{self.mass}'
        s += f',\n\t{self.position}'
        s += f',\n\t{self.inertial_matrix}'
        if self.inertial is not None:
            s += f',\n\t{self.inertial}'
        s += self.element_footer()
        return s

# Force Elements
class StructuralForce(Element):
    node: Node
    ftype: Literal['absolute', 'follower', 'total']
    position: Optional[Position] = None
    force_drive: Optional[List] = None # TODO: Needs TplDriveCaller
    force_orientation: Optional[Position] = None
    moment_orientation: Optional[Position] = None
    moment_drive: Optional[List] = None # TODO: Needs TplDriveCaller
    
    @model_validator(mode='after')
    def validate_fields_based_on_ftype(self):
        if self.ftype in ['absolute', 'follower']:
            # For 'absolute' or 'follower' types, position and force_drive are required
            if self.position is None:
                raise ValueError(f"{self.__class__.__name__}: position is required when ftype is {self.ftype}")
            if self.force_drive is None:
                raise ValueError(f"{self.__class__.__name__}: force_drive is required when ftype is {self.ftype}")
        return self    

    def element_type(self):
        return 'force'
    
    def __str__(self):
        s = f'{self.element_header()}, {self.ftype}'
        s += f',\n\t{self.node.idx}'
        if self.ftype == 'absolute' or self.ftype == 'follower':
            s += f',\n\t\tposition, {self.position}'
            s += f',\n\t\t'
            s += ', '.join(str(i) for i in self.force_drive)
        elif self.ftype == 'total':
            if self.position is not None:
                s += f',\n\t\tposition, {self.position}'
            if self.force_orientation is not None:
                s += f',\n\t\tforce orientation, {self.force_orientation}'
            if self.moment_orientation is not None:
                s += f',\n\t\tmoment orientation, {self.moment_orientation}'
            if self.force_drive is not None:
                s += f',\n\t\tforce, '
                s += ', '.join(str(i) for i in self.force_drive)
            if self.moment_drive is not None:
                s += f',\n\t\tmoment, '
                s += ', '.join(str(i) for i in self.moment_drive)
        s += self.element_footer()
        return s

class StructuralInternalForce(Element):
    nodes: List[Node]
    ftype: Literal['absolute', 'follower', 'total']
    positions: Optional[List[Position]] = None
    force_drive: Optional[List] = None  # TODO: Needs TplDriveCaller
    force_orientation: Optional[List[Position]] = None
    moment_orientation: Optional[List[Position]] = None
    moment_drive: Optional[List] = None  # TODO: Needs TplDriveCaller

    @model_validator(mode='after')
    def validate_fields_based_on_ftype(self):
        if self.ftype in ['absolute', 'follower']:
            if self.positions is None:
                raise ValueError(f"{self.__class__.__name__}: positions is required when ftype is {self.ftype}")
            if self.force_drive is None:
                raise ValueError(f"{self.__class__.__name__}: force_drive is required when ftype is {self.ftype}")
        elif self.ftype == 'total':
            if self.force_orientation is None or self.moment_orientation is None:
                raise ValueError(f"{self.__class__.__name__}: force_orientation and moment_orientation are required when ftype is total")
        return self

    @model_validator(mode='after')
    def validate_length_of_lists(self):
        if len(self.nodes) != 2:
            raise ValueError(f"{self.__class__.__name__}: nodes must have length 2")
        if self.positions is not None and len(self.positions) != 2:
            raise ValueError(f"{self.__class__.__name__}: positions must have length 2")
        if self.force_orientation is not None and len(self.force_orientation) != 2:
            raise ValueError(f"{self.__class__.__name__}: force_orientation must have length 2")
        if self.moment_orientation is not None and len(self.moment_orientation) != 2:
            raise ValueError(f"{self.__class__.__name__}: moment_orientation must have length 2")
        return self

    def element_type(self):
        return 'force'
    
    def __str__(self):
        s = f'{self.element_header()}, {self.ftype} internal'
        s += f',\n\t{self.nodes[0].idx}'
        if self.positions:
            s += f',\n\t\tposition, {self.positions[0]}'
        if self.ftype == 'total':
            if self.force_orientation:
                s += f',\n\t\tforce orientation, {self.force_orientation[0]}'
            if self.moment_orientation:
                s += f',\n\t\tmoment orientation, {self.moment_orientation[0]}'
        s += f',\n\t{self.nodes[1].idx}'
        if self.positions:
            s += f',\n\t\tposition, {self.positions[1]}'
        if self.ftype == 'total':
            if self.force_orientation:
                s += f',\n\t\tforce orientation, {self.force_orientation[1]}'
            if self.moment_orientation:
                s += f',\n\t\tmoment orientation, {self.moment_orientation[1]}'
            if self.force_drive:
                s += f',\n\t\tforce, '
                s += ', '.join(str(i) for i in self.force_drive)
            if self.moment_drive:
                s += f',\n\t\tmoment, '
                s += ', '.join(str(i) for i in self.moment_drive)
        else:  # ftype = { absolute|follower }
            s += f',\n\t\t'
            s += ', '.join(str(i) for i in self.force_drive)
        s += self.element_footer()
        return s

class StructuralCouple(Element):
    node: Node
    ctype: Literal['absolute', 'follower']
    position: Optional[Position] = None
    couple_drive: List # TODO: Needs TplDriveCaller

    def element_type(self):
        return 'couple'

    def __str__(self):
        s = f'{self.element_header()}, {self.ctype}'
        s += f',\n\t{self.node.idx}'
        if self.position:
            s += f',\n\t\tposition, {self.position}'
        s += f',\n\t\t'
        s += ', '.join(str(i) for i in self.couple_drive)
        s += self.element_footer()
        return s

class StructuralInternalCouple(Element):
    nodes: List[Node]
    ctype: Literal['absolute', 'follower']
    positions: Optional[List[Position]] = None
    couple_drive: List # TODO: Needs TplDriveCaller

    def element_type(self):
        return 'couple'
    
    @model_validator(mode='after')
    def validate_length_of_lists(self):
        if len(self.nodes) != 2:
            raise ValueError(f"{self.__class__.__name__}: nodes must have length 2")
        if self.positions is not None and len(self.positions) != 2:
            raise ValueError(f"{self.__class__.__name__}: positions must have length 2")
        return self
    
    def __str__(self):
        s = f'{self.element_header()}, {self.ctype} internal'
        s += f',\n\t{self.nodes[0].idx}'
        if self.positions:
            s += f',\n\t\tposition, {self.positions[0]}'
        s += f',\n\t{self.nodes[1].idx}'
        if self.positions:
            s += f',\n\t\tposition, {self.positions[1]}'
        s += f',\n\t\t'
        s += ', '.join(str(i) for i in self.couple_drive)
        s += self.element_footer()
        return s

# Joint Elements
class AngularAcceleration(Element):
    """
    This joint imposes the absolute angular acceleration of a node about a given axis.
    """

    node_label: Union[int, MBVar] # TODO: Take input as Node and use it's idx
    relative_direction: List[Union[float, MBVar]]
    acceleration: 'DriveCaller'

    def element_type(self):
        return 'joint'
    
    @field_validator('relative_direction')
    def validate_relative_direction(cls, v):
        Element.check_unit_vector3(v)
        return v
    
    def __str__(self):
        s = f'''{self.element_header()}, angular acceleration'''
        s += f''',\n\t{self.node_label}, {self.relative_direction}'''
        s += f''',\n\t{self.acceleration}'''
        s += self.element_footer()
        return s

class AngularVelocity(Element):
    """
    Represents a joint imposing the absolute angular velocity of a node about a given axis.
    """

    node_label: Union[int, MBVar]
    relative_direction: List[Union[float, MBVar]]
    velocity: 'DriveCaller'

    def element_type(self):
        return 'joint'

    @field_validator('relative_direction')
    def validate_relative_direction(cls, v):
        Element.check_unit_vector3(v)
        return v
    
    def __str__(self):
        s = f'''{self.element_header()}, angular velocity'''
        s += f''',\n\t{self.node_label}, '''
        s += ', '.join(str(i) for i in self.relative_direction)
        s += f''',\n\t{self.velocity}'''
        s += self.element_footer()
        return s
    
class AxialRotation(Element):
    """
    This joint is equivalent to a revolute hinge, but the angular velocity about axis 3 is imposed by means of the driver.
    """

    node_1_label: Union[int, MBVar]
    position_1: Position
    orientation_mat_1: Position
    node_2_label: Union[int, MBVar]
    position_2: Position
    orientation_mat_2: Position
    angular_velocity: 'DriveCaller'

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, axial rotation'
        s += f''',\n\t{self.node_1_label}'''
        s += f''',\n\t\tposition, {self.position_1}'''
        s += f''',\n\t\torientation, {self.orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}'''
        s += f''',\n\t\tposition, {self.position_2}'''
        s += f''',\n\t\torientation, {self.orientation_mat_2}'''
        s += f''',\n\t{self.angular_velocity}'''
        s += self.element_footer()
        return s
    
class Beam(Element):
    nodes: List[Node]
    positions: List[Position]
    orientations: List[Position]
    const_laws_orientations: List[Union[Position, Literal['same']]]
    const_laws: List[Union['ConstitutiveLaw', 'NamedConstitutiveLaw', Literal['same']]]
    custom_output: Optional[List] = None # TODO: Add custom output class

    @model_validator(mode='after')
    def validate_const_laws(self):
        assert self.const_laws_orientations[0] != 'same', (
            f'\n-------------------\nERROR:' + 
            f'{self.__class__.__name__}: the first constitutive law orientation must not be "same";\n' +
            f'\n-------------------\n')
        assert self.const_laws[0] != 'same', (
            f'\n-------------------\nERROR:' + 
            f'{self.__class__.__name__}: the first constitutive law must not be "same";\n' +
            f'\n-------------------\n')
        return self

    @model_validator(mode='after')
    def validate_lengths(self):
        assert len(self.nodes) == 3 or len(self.nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a beam with ' + str(len(self.nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(self.nodes) == len(self.positions), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(self.nodes)) +
            ' nodes and ' + str(len(self.positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert len(self.nodes) == len(self.orientations), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(self.nodes)) +
            ' nodes and ' + str(len(self.orientations)) + ' relative orientations;\n' +
            '\n-------------------\n')
        assert len(self.const_laws) == len(self.const_laws_orientations), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(self.const_laws)) +
            ' constitutive laws and ' + str(len(self.const_laws_orientations)) + ' constitutive law orientations;\n' +
            '\n-------------------\n')
        return self
    
    @model_validator(mode='before')
    def adjust_const_laws_for_two_nodes(cls, values):
        if len(values['nodes']) == 2:
            # Convert single instance to list for two-node beams
            if 'const_laws' in values and not isinstance(values['const_laws'], list):
                values['const_laws'] = [values['const_laws']]
            if 'const_laws_orientations' in values and not isinstance(values['const_laws_orientations'], list):
                values['const_laws_orientations'] = [values['const_laws_orientations']]
        return values
    
    def element_type(self):
        if len(self.nodes) == 3:
            return 'beam3'
        else:
            return 'beam2'
        
    def __str__(self):
        s = f'{self.element_header()}'
        for (node, position, orientation) in zip(self.nodes, self.positions, self.orientations):
            s += f',\n\t{node.idx}'
            s += f',\n\t\tposition, {position}'
            s += f',\n\t\torientation, {orientation}'
        for (cl_or, cl) in zip(self.const_laws_orientations, self.const_laws):
            s += f',\n\t{cl_or}'
            s += f',\n\t{cl}'
        if self.custom_output is not None:
            s += f',\n\tcustom output, ' 
            s += f', '.join(str(i) for i in self.custom_output)
        s += self.element_footer()
        return s
        
class BeamSlider(Element):
    """
    This joint implements a slider, e.g. it constrains a structural node on a string of three-node beams.
    """

    slider_node_label: Union[int, MBVar]
    position: Position
    orientation: Optional[Position]
    slider_type: Optional[str] = None  # should be one of 'spherical', 'classic', or 'spline'
    beam_number: Union[int, MBVar]
    three_node_beam: 'Beam'
    first_node_offset: Union[str, Position]
    first_node_orientation: Optional[Union[str, Position]]
    mid_node_offset: Position
    mid_node_orientation: Optional[Union[str, Position]]
    end_node_offset: Position
    end_node_orientation: Optional[Union[str, Position]]
    initial_beam: Optional[Beam]
    initial_node: Optional[Node]
    smearing_factor: Optional[Union[float, MBVar, int]]

    def element_type(self):
        return 'joint'

    @field_validator('slider_type')
    def validate_slider_type(cls, v):
        allowed_types = {'spherical', 'classic', 'spline'}
        if v is not None and v not in allowed_types:
            raise ValueError(f"slider_type must be one of {allowed_types}")
        return v
    
    def __str__(self):
        s = f'{self.element_header()}, kinematic'
        s += f''',\n\t{self.slider_node_label}'''
        s += f''',\n\t\t{self.position}'''
        if self.orientation is not None:
            s += f''',\n\t\thinge, {self.orientation}'''
        if self.slider_type is not None:
            s += f''',\n\ttype, {self.slider_type}'''
        s += f''',\n\t{self.beam_number}'''
        s += f''',\n\t\t{self.three_node_beam}'''
        if s.endswith(';\n'):
            s = s[:-2] # Remove the last two characters
        if isinstance(self.first_node_offset, str) and self.first_node_offset == 'same':
            s += f''',\n\t\t\tsame'''
        else:
            s += f''',\n\t\t\t{self.first_node_offset}'''
        if self.first_node_orientation is not None:
            s += f''',\n\t\thinge, {self.first_node_orientation}'''
        s += f''',\n\t\t\t{self.mid_node_offset}'''
        if self.mid_node_orientation is not None:
            s += f''',\n\t\thinge, {self.mid_node_orientation}'''
        s += f''',\n\t\t\t{self.end_node_offset}'''
        if self.end_node_orientation is not None:
            s += f''',\n\t\thinge, {self.end_node_orientation}'''
        if self.initial_beam is not None:
            s += f''',\n\tinitial beam, {self.initial_beam}'''
            if s.endswith(';\n'):
                s = s[:-2]  # Remove the last two characters
        if self.initial_node is not None:
            s += f''',\n\tinitial node, {self.initial_node}'''
            if s.endswith(';\n'):
                s = s[:-2]  # Remove the last two characters
        if self.smearing_factor is not None:
            s += f''',\n\tsmearing, {self.smearing_factor}'''
        s += self.element_footer()
        return s

class Brake(Element):
    """
    This element models a wheel brake, i.e., a constraint that applies a frictional internal torque between two
    nodes about an axis. The frictional torque depends on the normal force that is applied as an external
    input by means of the same friction models implemented for regular joints.
    """

    node_1_label: Union[int, MBVar]
    position_1: Position
    orientation_mat_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    position_2: Position
    orientation_mat_2: Optional[Position] = None
    average_radius: Union[float, MBVar]
    preload: Optional[Union[float, MBVar, int]] = None
    friction_model: str  # TODO: Implement FrictionModel class
    shape_function: str  # TODO: Implement ShapeFunction class
    normal_force: 'DriveCaller'

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, brake'
        s += f''',\n\t{self.node_1_label}, {self.position_1}'''
        if self.orientation_mat_1 is not None:
            s += f''',\n\thinge, {self.orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}, {self.position_2}'''
        if self.orientation_mat_2 is not None:
            s += f''',\n\thinge, {self.orientation_mat_2}'''
        s += f''',\n\tfriction, {self.average_radius}'''
        if self.preload is not None:
            s += f''',\n\t\tpreload, {self.preload}'''
        s += f''',\n\t\t{self.friction_model}'''
        s += f''',\n\t\t{self.shape_function}'''
        s += f''',\n\t{self.normal_force}'''
        s += self.element_footer()
        return s
        
class CardanoHinge(Element):
    '''
    This joint implements a Cardano's joint, also known as Hooke's joint or Universal joint, which is made
    of a sequence of two revolute hinges orthogonal to each other, one about relative axis 2 and one about
    relative axis 3 of the reference systems defined by the two orientation statements. 
    
    In other words, this joint constrains the relative axis 3 of node 1 to be always orthogonal to the 
    relative axis 2 of node 2. As a result, torque is transmitted about axis 1 of both nodes. The relative 
    position is constrained as well.

    Note: This joint does not represent a constant velocity joint, so, when a steady deflection between 
    the two nodes is present, a constant velocity about axis 1 of one node results in an oscillating 
    velocity about axis 1 for the other node.
    '''

    node_1: Node
    position_1: Position
    orientation_mat_1: Optional[Position] = None
    node_2: Node
    position_2: Position
    orientation_mat_2: Optional[Position] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, cardano hinge'
        s += f',\n\t{self.node_1.idx}'
        s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2.idx}'
        s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += self.element_footer()
        return s
    
class CardanoPin(Element):
    """
    This joint implements a 'Cardano' joint between a node and the ground.
    The absolute position is also constrained.
    """

    node_label: Union[int, MBVar]
    position: Position
    orientation_mat: Optional[Position] = None
    absolute_pin_position: Position
    absolute_pin_orientation_mat: Optional[Position] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, cardano pin'
        s += f',\n\t{self.node_label},'
        s += f'\n\t\tposition, {self.position}'
        if self.orientation_mat is not None:
            s += f',\n\t\torientation, {self.orientation_mat}'
        s += f',\n\tposition, {self.absolute_pin_position}'
        if self.absolute_pin_orientation_mat is not None:
            s += f',\n\torientation, {self.absolute_pin_orientation_mat}'
        s += self.element_footer()
        return s

class CardanoRotation(Element):
    """
    This joint implements a 'Cardano' joint, which is made of a sequence of two orthogonal revolute hinges.
    The relative position is not constrained.
    """

    node_1_label: Union[int, MBVar]
    orientation_mat_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    orientation_mat_2: Optional[Position] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, cardano rotation'
        s += f',\n\t{self.node_1_label}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2_label}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += self.element_footer()
        return s

class DeformableAxial(Element):
    """
    This joint implements a configuration dependent moment that is exchanged between two nodes about
    an axis rigidly attached to the first node. 
    """

    node_1_label: Union[int, MBVar]
    position_1: Optional[Position] = None
    orientation_mat_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    position_2: Optional[Position] = None
    orientation_mat_2: Optional[Position] = None
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, deformable axial'
        s += f',\n\t{self.node_1_label}'
        if self.position_1 is not None:
            s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2_label}'
        if self.position_2 is not None:
            s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s

class DeformableHinge(Element):
    """
    This joint implements a configuration dependent moment that is exchanged between two nodes. The
    moment may depend, by way of a generic 3D constitutive law, on the relative orientation and angular
    velocity of the two nodes, expressed in the reference frame of node 1.
    """

    node_1: Node
    position_1: Optional[Position] = None
    orientation_mat_1: Optional[Position] = None
    node_2: Node
    position_2: Optional[Position] = None
    orientation_mat_2: Optional[Position] = None
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']
    orientation_desc: Optional[Literal['euler123', 'euler313', 'euler321', 'orientation vector', 'orientation matrix']] = None

    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 3D constitutive law with law_type 'D3_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, deformable hinge'
        s += f',\n\t{self.node_1.idx}'
        if self.position_1 is not None:
            s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2.idx}'
        if self.position_2 is not None:
            s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += f',\n\t{self.const_law}'
        if self.orientation_desc is not None:
            s += f''',\n\torientation description, {self.orientation_desc}'''
        s += self.element_footer()
        return s

class Distance(Element):
    """
    This joint forces the distance between two points, each relative to a node, to assume the value indicated
    by the drive. If no offset is given, the points are coincident with the node themselves.
    """

    node_1_label: Union[int, MBVar]
    position_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    position_2: Optional[Position] = None
    distance: Union['DriveCaller', str]

    @field_validator('distance')
    def validate_distance(cls, value):
        if isinstance(value, str):
            if value != "from nodes":
                raise ValueError('Invalid value for distance. It must be "from nodes" if a string.')
        return value
    
    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, distance'
        s += f''',\n\t{self.node_1_label}'''
        if self.position_1 is not None:
            s += f''', position, {self.position_1}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.position_2 is not None:
            s += f''', position, {self.position_2}'''
        s += f''',\n\t{self.distance}'''
        s += self.element_footer()
        return s

class DriveDisplacement(Element):
    '''
    This joint imposes the relative position between two points optionally offset from two structural nodes,
    in the form of a vector that expresses the direction of the displacement in the reference frame of node 1,
    whose amplitude is defined by a drive.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position
    node_2_label: Union[int, MBVar]
    position_2: Position
    relative_position: 'List'

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, drive displacement'
        s += f''',\n\t{self.node_1_label}, {self.position_1}'''
        s += f''',\n\t{self.node_2_label}, {self.position_2}'''
        s += f''',\n\t{", ".join(str(i) for i in self.relative_position)}'''
        s += self.element_footer()
        return s
    
class DriveDisplacementPin(Element):
    '''
    This joint imposes the relative position between two points optionally offset from two structural nodes,
    in the form of a vector that expresses the direction of the displacement in the reference frame of node 1,
    whose amplitude is defined by a drive.
    '''

    node_label: Union[int, MBVar]
    node_offset: Position
    offset: Position
    position: 'List'

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, drive displacement pin'
        s += f''',\n\t{self.node_label}, {self.node_offset}'''
        s += f''',\n\t{self.offset}'''
        s += f''',\n\t{", ".join(str(i) for i in self.position)}'''
        s += self.element_footer()
        return s

class DriveHinge(Element):
    '''
    This joint imposes the relative orientation between two nodes, in the form of a rotation about an axis
    whose amplitude is defined by a drive.
    '''

    node_1_label: Union[int, MBVar]
    relative_orientation_mat_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    relative_orientation_mat_2: Optional[Position] = None
    hinge_orientation: 'List'

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, drive hinge'
        s += f''',\n\t{self.node_1_label}'''
        if self.relative_orientation_mat_1 is not None:
            s += f''', orientation, {self.relative_orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.relative_orientation_mat_2 is not None:
            s += f''', orientation, {self.relative_orientation_mat_2}'''
        s += f''',\n\t{", ".join(str(i) for i in self.hinge_orientation)}'''
        s += self.element_footer()
        return s
    
class GimbalRotation(Element):
    '''
    A homokinetic joint without position constraints; this joint, in conjunction with a spherical hinge 
    joint, should be used to implement an ideal tiltrotor gimbal instead of a cardano
    rotation. It is equivalent to a series of two Cardano's joints (the cardano hinge) rotated 90 degrees
    apart, each accounting for half the relative rotation between axis 3 of each side of the joint.
    '''

    node_1_label: Union[int, MBVar]
    relative_orientation_mat_1: Optional[Union[Position]] = None
    node_2_label: Union[int, MBVar]
    relative_orientation_mat_2: Optional[Union[Position]] = None
    orientation_desc: Optional[Literal['euler123', 'euler313', 'euler321', 'orientation vector', 'orientation matrix']] = None
    """The type of orientation description"""

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, gimbal rotation'
        s += f''',\n\t{self.node_1_label}'''
        if self.relative_orientation_mat_1 is not None:
            s += f''', orientation, {self.relative_orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.relative_orientation_mat_2 is not None:
            s += f''', orientation, {self.relative_orientation_mat_2}'''
        if self.orientation_desc is not None:
            s += f''',\n\torientation description, {self.orientation_desc}'''
        s += self.element_footer()
        return s

class ImposedDisplacement(Element):
    '''
    This joint imposes the relative position between two points, optionally offset from two structural nodes,
    along a given direction that is rigidly attached to the first node. The amplitude of the displacement is
    defined by a drive.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position
    node_2_label: Union[int, MBVar]
    position_2: Position
    direction: List[Union[float, MBVar]]
    relative_position: 'DriveCaller'

    @field_validator('direction')
    def validate_direction(cls, v):
        Element.check_unit_vector3(v)
        return v

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, imposed displacement'
        s += f''',\n\t{self.node_1_label}, {self.position_1}'''
        s += f''',\n\t{self.node_2_label}, {self.position_2}'''
        s += f''',\n\t{self.direction}'''
        s += f''',\n\t{self.relative_position}'''
        s += self.element_footer()
        return s

class ImposedDisplacementPin(Element):
    '''
    This joint imposes the absolute displacement of a point optionally offset from a structural node, along
    a direction defined in the absolute reference frame. The amplitude of the displacement is defined by a
    drive.
    '''

    node_label: Union[int, MBVar]
    node_offset: Position
    offset: Position
    direction: List[Union[float, MBVar]]
    position: 'DriveCaller'

    @field_validator('direction')
    def validate_direction(cls, v):
        Element.check_unit_vector3(v)
        return v

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, imposed displacement pin'
        s += f''',\n\t{self.node_label}, {self.node_offset}'''
        s += f''',\n\t{self.offset}'''
        s += f''',\n\t{self.direction}'''
        if self.position.idx is not None and self.position.idx >= 0:
            s += f''',\n\treference, {self.position.idx}'''
        else:
            s += f''',\n\t{self.position}'''
        s += self.element_footer()
        return s
    
class InLine(Element):
    '''
    This joint forces a point relative to the second node to move along a line attached to the first node.
    '''
    
    node_1_label: Union[int, MBVar]
    position: Optional[Position] = None
    orientation: Optional[Union[Position, List]] = None
    node_2_label: Union[int, MBVar]
    offset: Optional[Position] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, in line'
        s += f''',\n\t{self.node_1_label}'''
        if self.position is not None:
            s += f''', position, {self.position}'''
        if self.orientation is not None:
            s += f'''\n\t, orientation, {self.orientation}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.offset is not None:
            s += f''', offset, {self.offset}'''
        s += self.element_footer()
        return s
    
class InPlane(Element):
    '''
    This joint forces a point relative to the second node to move in a plane attached to the first node.
    '''
    
    node_1_label: Union[int, MBVar]
    position: Optional[Position] = None
    relative_direction: List[Union[float, MBVar]]
    node_2_label: Union[int, MBVar]
    offset: Optional[Position] = None

    @field_validator('relative_direction')
    def validate_relative_direction(cls, v):
        Element.check_unit_vector3(v)
        return v

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, in plane'
        s += f''',\n\t{self.node_1_label}'''
        if self.position is not None:
            s += f''', position, {self.position}'''
        s += f''',\n\t{self.relative_direction}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.offset is not None:
            s += f''', offset, {self.offset}'''
        s += self.element_footer()
        return s
    
class LinearAcceleration(Element):
    '''
    This joint imposes the absolute linear acceleration of a node along a given axis.
    '''
    
    node_label: Union[int, MBVar]
    relative_direction: List[Union[float, MBVar]]
    acceleration: 'DriveCaller'

    @field_validator('relative_direction')
    def validate_relative_direction(cls, v):
        Element.check_unit_vector3(v)
        return v

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, linear acceleration'
        s += f''',\n\t{self.node_label}'''
        s += f''',\n\t {self.relative_direction}'''
        if self.acceleration.idx is not None and self.acceleration.idx >= 0:
            s += f''',\n\treference, {self.acceleration.idx}'''
        else:
            s += f''',\n\t{self.acceleration}'''
        s += self.element_footer()
        return s
    
class LinearVelocity(Element):
    '''
    This joint imposes the absolute linear velocity of a node along a given axis.
    '''
    
    node_label: Union[int, MBVar]
    relative_direction: List[Union[float, MBVar]]
    velocity: 'DriveCaller'

    @field_validator('relative_direction')
    def validate_relative_direction(cls, v):
        Element.check_unit_vector3(v)
        return v

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, linear velocity'
        s += f''',\n\t{self.node_label}'''
        s += f''',\n\t {self.relative_direction}'''
        if self.velocity.idx is not None and self.velocity.idx >= 0:
            s += f''',\n\treference, {self.velocity.idx}'''
        else:
            s += f''',\n\t{self.velocity}'''
        s += self.element_footer()
        return s
    
class Modal(Element):
    pass

class PlaneDisplacement(Element):
    '''
    This joint allows two nodes to move in the common relative 1-2 plane and to rotate about the common
    relative axis 3.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position
    orientation_mat_1: Optional[Union[Position, List]] = None
    node_2_label: Union[int, MBVar]
    position_2: Position
    orientation_mat_2: Optional[Union[Position, List]] = None

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, plane displacement'
        s += f''',\n\t{self.node_1_label}, position, {self.position_1}'''
        if self.orientation_mat_1 is not None:
            s += f''',\n\torientation, {self.orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}, position, {self.position_2}'''
        if self.orientation_mat_2 is not None:
            s += f''',\n\torientation, {self.orientation_mat_2}'''
        s += self.element_footer()
        return s
    
class PlaneDisplacementPin(Element):
    '''
    This joint allows a node to move in the relative 1–2 plane and to rotate about the relative axis 3 with
    respect to an absolute point and plane.
    '''

    node_label: Union[int, MBVar]
    relative_offset: Position
    relative_orientation_mat: Optional[Union[Position, List]] = None
    absolute_pin_position: Position
    absolute_pin_orientation_mat: Optional[Union[Position, List]] = None

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, plane displacement pin'
        s += f''',\n\t{self.node_label}'''
        s += f''',\n\t\tposition, {self.relative_offset}'''
        if self.relative_orientation_mat is not None:
            s += f''',\n\t\torientation, {self.relative_orientation_mat}'''
        s += f''',\n\tposition, {self.absolute_pin_position}'''
        if self.absolute_pin_orientation_mat is not None:
            s += f''',\n\torientation, {self.absolute_pin_orientation_mat}'''
        s += self.element_footer()
        return s

class Prismatic(Element):
    '''
    This joints constrains the relative orientation of two nodes, so that their orientations remain parallel.
    The relative position is not constrained. The initial orientation of the joint must be compatible: use the
    orientation keyword to assign the joint initial orientation.
    '''

    node_1_label: Union[int, MBVar]
    relative_orientation_mat_1: Optional[Union[Position, List]] = None
    node_2_label: Union[int, MBVar]
    relative_orientation_mat_2: Optional[Union[Position, List]] = None

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, prismatic'
        s += f''',\n\t{self.node_1_label}'''
        if self.relative_orientation_mat_1 is not None:
            s += f''', orientation, {self.relative_orientation_mat_1}'''
        s += f''',\n\t{self.node_2_label}'''
        if self.relative_orientation_mat_2 is not None:
            s += f''', orientation, {self.relative_orientation_mat_2}'''
        s += self.element_footer()
        return s

class RevoluteHinge(Element):
    '''
    This joint only allows the relative rotation of two nodes about a given axis, which is axis 3 in the reference
    systems defined by the two orientation statements.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position
    orientation_mat_1: Optional[Position] = None
    node_2_label: Union[int, MBVar]
    position_2: Position
    orientation_mat_2: Optional[Position] = None
    initial_theta: Optional[Union[float, MBVar]] = None
    friction: Optional[Union[float, MBVar]] = None
    preload: Optional[Union[float, MBVar]] = None
    friction_model: Optional[str] = None # TODO: Define FrictionModel
    shape_function: Optional[str] = None # TODO: Define ShapeFunction

    def element_type(self):
        return 'joint'
    
    @model_validator(mode='after')
    def check_friction_parameters(self):
        if self.friction is not None:
            if self.friction_model is None or self.shape_function is None:
                raise ValueError("When 'friction' is specified, 'friction_model' and 'shape_function' must also be specified.")
            # 'preload' is optional when 'friction' is specified
        else:
            # If 'friction' is not specified, none of the friction-related parameters should be specified
            if any(param is not None for param in [self.preload, self.friction_model, self.shape_function]):
                raise ValueError("If 'friction' is not specified, 'preload', 'friction_model', and 'shape_function' should not be specified.")
        return self

    def __str__(self):
        s = f'{self.element_header()}, revolute hinge'
        s += f',\n\t{self.node_1_label}'
        s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2_label}'
        s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        if self.initial_theta is not None:
            s += f',\n\tinitial theta, {self.initial_theta}'
        if self.friction is not None:
            s += f',\n\tfriction, {self.friction}'
            if self.preload is not None:
                s += f',\n\t\tpreload, {self.preload}'
            s += f',\n\t\t{self.friction_model}'
            s += f',\n\t\t{self.shape_function}'
        s += self.element_footer()
        return s

class RevolutePin(Element):
    """
    This joint only allows the absolute rotation of a node about a given axis, which is axis 3 in the reference
    systems defined by the two orientation statements.
    """

    node_label: Union[int, MBVar]
    relative_offset: Position
    relative_orientation_mat: Optional[Union[Position, list]] = None
    absolute_pin_position: Position
    absolute_pin_orientation_mat: Optional[Union[Position, list]] = None
    initial_theta: Optional[Union[float, MBVar]] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, revolute pin'
        s += f',\n\t{self.node_label}'
        s += f',\n\t\tposition, {self.relative_offset}'
        if self.relative_orientation_mat is not None:
            s += f',\n\t\torientation, {self.relative_orientation_mat}'
        s += f',\n\tposition, {self.absolute_pin_position}'
        if self.absolute_pin_orientation_mat is not None:
            s += f',\n\torientation, {self.absolute_pin_orientation_mat}'
        if self.initial_theta is not None:
            s += f',\n\tinitial theta, {self.initial_theta}'
        s += self.element_footer()
        return s
    
class RevoluteRotation(Element):
    '''
    This joint allows the relative rotation of two nodes about a given axis, which is axis 3 in the reference
    systems defined by the two orientation statements. The relative position is not constrained.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Optional[Position] = None
    orientation_mat_1: Optional[Union[Position, list]] = None
    node_2_label: Union[int, MBVar]
    position_2: Optional[Position] = None
    orientation_mat_2: Optional[Union[Position, list]] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, revolute rotation'
        s += f',\n\t{self.node_1_label}'
        if self.position_1 is not None:
            s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2_label}'
        if self.position_2 is not None:
            s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += self.element_footer()
        return s

class Rod(Element):
    '''
    The rod element represents a force between two nodes that depends on the relative position and velocity
    of two points, each rigidly attached to a structural node. The direction of the force is also based on
    the relative position of the points: it is the line that passes through them. If no offset is defined, the
    points are the nodes themselves.
    '''

    node_1: Node
    position_1: Optional[Position] = None
    node_2: Node
    position_2: Optional[Position] = None
    rod_length: Union[float, MBVar, Literal['from nodes']]  # Can be a float or 'from nodes'
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']

    def element_type(self):
        return 'joint'
        
    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 1D constitutive law with law_type 'SCALAR_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")

    def __str__(self):
        s = f'{self.element_header()}, rod'
        s += f',\n\t{self.node_1.idx}'
        if self.position_1 is not None:
            s += f',\n\t\tposition, {self.position_1}'
        s += f',\n\t{self.node_2.idx}'
        if self.position_2 is not None:
            s += f',\n\t\tposition, {self.position_2}'
        s += f',\n\t{self.rod_length}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s
    
class RodWithOffset(Element):
    '''
    Analogous to the rod joint with the optional offsets.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position  # Required
    node_2_label: Union[int, MBVar]
    position_2: Position  # Required
    rod_length: Union[float, MBVar, str]  # Can be a float, MBVar or 'from nodes'
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']  # Should be a 1D constitutive law

    def element_type(self):
        return 'joint'

    @field_validator('rod_length')
    def validate_rod_length(cls, v):
        if isinstance(v, str):
            if v.lower() != 'from nodes':
                raise ValueError("rod_length must be a float or the string 'from nodes'")
            return v.lower()
        else:
            return v

    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 1D constitutive law with law_type 'SCALAR_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")

    def __str__(self):
        s = f'{self.element_header()}, rod with offset'
        s += f',\n\t{self.node_1_label}'
        s += f',\n\t\t{self.position_1}'
        s += f',\n\t{self.node_2_label}'
        s += f',\n\t\t{self.position_2}'
        if isinstance(self.rod_length, str) and self.rod_length == 'from nodes':
            s += f',\n\tfrom nodes'
        else:
            s += f',\n\t{self.rod_length}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s

class RodBezier(Element):
    '''
    This joint, in analogy with the rod joint, represents a force that acts between two points each rigidly
    attached to a structural node.

    The force on node 1 acts along the line connecting the insertion point, defined in the reference frame of
    the node by <relative_offset_1> and a first intermediate point defined, also in the reference frame of
    the node, by <relative_offset_2>. In the same way, the force on node 2 is applied at the insertion
    point of the element, defined in the reference frame of node 2 by <relative_offset_4> and acts along
    the line connecting a second intermediate point defined by <relative_offset_3>.

    The element internally is represented as a Bézier spline of order 3, which length and instantaneous
    lengthening velocity are calculated using Gauss-Legendre quadrature.

    The absolute value of the force depends on the strain and strain rate of the curve as in the standard rod
    element as determined by the ConstitutiveLaw<1D> <const_law>.
    '''

    node_1_label: Union[int, MBVar]
    position_1: Position
    position_2: Position
    node_2_label: Union[int, MBVar]
    position_3: Position
    position_4: Position
    rod_length: Union[float, MBVar, str]  # Can be a float or 'from nodes'
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']  # Should be a 1D constitutive law
    integration_order: int = 2  # Defaults to 2
    integration_segments: int = 3  # Defaults to 3

    def element_type(self):
        return 'joint'

    @field_validator('rod_length')
    def validate_rod_length(cls, v):
        if isinstance(v, str):
            if v.lower() != 'from nodes':
                raise ValueError("rod_length must be a float or the string 'from nodes'")
            return v.lower()
        else:
            return v

    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 1D constitutive law with law_type 'SCALAR_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")

    @field_validator('integration_order')
    def validate_integration_order(cls, v):
        if not isinstance(v, int):
            raise TypeError("integration_order must be an integer")
        if not (1 <= v <= 10):
            raise ValueError("integration_order must be between 1 and 10")
        return v

    @field_validator('integration_segments')
    def validate_integration_segments(cls, v):
        if not isinstance(v, int):
            raise TypeError("integration_segments must be an integer")
        if v <= 0:
            raise ValueError("integration_segments must be a positive integer")
        return v

    def __str__(self):
        s = f'{self.element_header()}, rod bezier'
        s += f',\n\t{self.node_1_label}'
        s += f',\n\t\t{self.position_1}'
        s += f',\n\t\t{self.position_2}'
        s += f',\n\t{self.node_2_label}'
        s += f',\n\t\t{self.position_3}'
        s += f',\n\t\t{self.position_4}'
        if isinstance(self.rod_length, str) and self.rod_length == 'from nodes':
            s += f',\n\tfrom nodes'
        else:
            s += f',\n\t{self.rod_length}'
        # Include integration order only if it differs from the default value of 2
        if self.integration_order != 2:
            s += f',\n\tintegration order, {self.integration_order}'
        # Include integration segments only if it differs from the default value of 3
        if self.integration_segments != 3:
            s += f',\n\tintegration segments, {self.integration_segments}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s
    
class SphericalPin(Element):
    '''
    This joint constrains the absolute position of a node; the relative orientation is not constrained.
    **Note**: This joint is equivalent to a spherical hinge when one node is grounded.
    '''

    node_label: Union[int, 'MBVar']
    position: Optional[Position] = None
    orientation_mat: Optional[Union[Position, list]] = None
    absolute_pin_position: Position
    absolute_orientation_mat: Optional[Union[Position]] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, spherical pin'
        s += f',\n\t{self.node_label}'
        if self.position is not None:
            s += f',\n\t\tposition, {self.position}'
        if self.orientation_mat is not None:
            s += f',\n\t\torientation, {self.orientation_mat}'
        s += f',\n\tposition, {self.absolute_pin_position}'
        if self.absolute_orientation_mat is not None:
            s += f',\n\torientation, {self.absolute_orientation_mat}'
        s += self.element_footer()
        return s
    
class ViscousBody(Element):
    '''
    This element defines a force and a moment that depend on the absolute linear and angular velocity of
    a body, projected in the reference frame of the node itself. The force and moment are defined as a 6D
    viscous constitutive law.
    '''

    node_label: Union[int, MBVar]
    position: Optional[Position] = None
    orientation_mat: Optional[Union[Position, list]] = None
    const_law: Union['ConstitutiveLaw','NamedConstitutiveLaw']  # Should be a 6D constitutive law

    def element_type(self):
        return 'joint'
    
    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 6D constitutive law with law_type 'D6_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")
    
    def __str__(self):
        s = f'{self.element_header()}, viscous body'
        s += f',\n\t{self.node_label}'
        if self.position is not None:
            s += f',\n\t\tposition, {self.position}'
        if self.orientation_mat is not None:
            s += f',\n\t\torientation, {self.orientation_mat}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s

class Clamp(Element):
    node: Node
    position: Union[Position, Literal['node']]
    orientation_mat: Union[List, Literal['node']]

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, clamp, {self.node.idx}'
        s += f',\n\tposition, {self.position}'
        s += f',\n\torientation, {self.orientation_mat}'
        s += self.element_footer()
        return s

class TotalJoint(Element):
    nodes: List[Node]
    positions: Optional[List[Position]] = None
    position_orientations: Optional[List[Position]] = None
    rotation_orientations: Optional[List[Position]] = None
    position_status: Optional[List[Union[Literal['active', 'inactive', 'position', 'velocity'], bool]]] = None
    orientation_status: Optional[List[Union[Literal['active', 'inactive', 'rotation', 'angular velocity'], bool]]] = None
    position_drive: Optional[List] = None # TODO: Needs tpl drive
    orientation_drive: Optional[List] = None # TODO: Needs tpl drive

    @model_validator(mode='after')
    def validate_total_joint(self):
        assert len(self.nodes) == 2, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(self.nodes)) +
            ' nodes;\n' +
            '\n-------------------\n')
        
        if self.positions is not None:
            assert len(self.nodes) == len(self.positions), (  # either both or none must be given
                '\n-------------------\nERROR:' +
                ' defining a total joint with ' + str(len(self.nodes)) +
                ' nodes and ' + str(len(self.positions)) + ' relative positions;\n' +
                '\n-------------------\n')
        else:
            assert self.position_orientations is None and self.rotation_orientations is None, (
                '\n-------------------\nERROR:' +
                ' position orientations and rotation orientations cannot be given if positions is not given' +
                '\n-------------------\n')
            
        if self.position_orientations is not None:
            assert len(self.nodes) == len(self.position_orientations), (  # either both or none must be given
                '\n-------------------\nERROR:' +
                ' defining a total joint with ' + str(len(self.nodes)) +
                ' nodes and ' + str(len(self.position_orientations)) + ' position orientations;\n' +
                '\n-------------------\n')
        if self.rotation_orientations is not None:
            assert len(self.nodes) == len(self.rotation_orientations), (  # either both or none must be given
                '\n-------------------\nERROR:' +
                ' defining a total joint with ' + str(len(self.nodes)) +
                ' nodes and ' + str(len(self.rotation_orientations)) + ' rotation orientations;\n' +
                '\n-------------------\n')
        if self.position_status is not None:
            assert len(self.position_status) == 3, (
                '\n-------------------\nERROR:' +
                ' defining a total joint with ' + str(len(self.position_status)) +
                ' position status;\n' +
                '\n-------------------\n')
        if self.orientation_status is not None:
            assert len(self.orientation_status) == 3, (
                '\n-------------------\nERROR:' +
                ' defining a total joint with ' + str(len(self.orientation_status)) +
                ' orientation status;\n' +
                '\n-------------------\n')
        return self

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, total joint'
        s += f',\n\t{self.nodes[0].idx}'
        if self.positions is not None:
            s += f',\n\t\tposition, {self.positions[0]}'
        if self.position_orientations is not None:
            s += f',\n\t\tposition orientation, {self.position_orientations[0]}'
        if self.rotation_orientations is not None:
            s += f',\n\t\trotation orientation, {self.rotation_orientations[0]}'
        s += f',\n\t{self.nodes[1].idx}'
        if self.positions is not None:
            s += f',\n\t\tposition, {self.positions[1]}'
        if self.position_orientations is not None:
            s += f',\n\t\tposition orientation, {self.position_orientations[1]}'
        if self.rotation_orientations is not None:
            s += f',\n\t\trotation orientation, {self.rotation_orientations[1]}'
        if self.position_status is not None:
            s += f',\n\tposition constraint, '
            s += ', '.join(str(ps) for ps in self.position_status)
            if self.position_drive is not None:
                s += f',\n\t\t'
                s += ', '.join(str(i) for i in self.position_drive)
        if self.orientation_status is not None:
            s += f',\n\torientation constraint, '
            s += ', '.join(str(os) for os in self.orientation_status)
            if self.orientation_drive is not None:
                s += f',\n\t\t'
                s += ', '.join(str(i) for i in self.orientation_drive)
        s += self.element_footer()
        return s

class TotalPinJoint(Element):
    node: Node
    rel_position: Optional[Position] = None
    rel_position_orientation: Optional[Position] = None
    rel_rotation_orientation: Optional[Position] = None
    abs_position: Optional[Position] = None
    abs_position_orientation: Optional[Position] = None
    abs_rotation_orientation: Optional[Position] = None
    position_status: Optional[List[bool]] = None
    orientation_constraints: Optional[List[bool]] = None
    position_status: Optional[List[Union[Literal['active', 'inactive', 'position', 'velocity'], bool]]] = None
    orientation_status: Optional[List[Union[Literal['active', 'inactive', 'rotation', 'angular velocity'], bool]]] = None
    position_drive: Optional[List] = None # TODO: Needs tpl drive
    orientation_drive: Optional[List] = None # TODO: Needs tpl drive

    @model_validator(mode='after')
    def validate_total_pin_joint(self):
        if self.rel_position is None:
            assert self.rel_position_orientation is None and self.rel_rotation_orientation is None, (
                f'\n-------------------\nERROR:' +
                f' {self.__class__.__name__}: relative position orientation and rotation ' +
                f'orientation cannot be given if relative position is not given' +
                f'\n-------------------\n')
        if self.abs_position is None:
            assert self.abs_position_orientation is None and self.abs_rotation_orientation is None, (
                f'\n-------------------\nERROR:' +
                f' {self.__class__.__name__}: absolute position orientation and rotation ' +
                f'orientation cannot be given if absolute position is not given' +
                f'\n-------------------\n')
        assert len(self.position_status) == 3, (
            f'\n-------------------\nERROR:' +
            f' {self.__class__.__name__}: position status must be given as a list of 3 statuses' +
            f'\n-------------------\n')
        assert len(self.orientation_status) == 3, (
            f'\n-------------------\nERROR:' +
            f' {self.__class__.__name__}: orientation status must be given as a list of 3 statuses' +
            f'\n-------------------\n')
        return self

    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, total pin joint'
        s += f',\n\t{self.node.idx}'
        if self.rel_position is not None:
            s += f',\n\t\tposition, {self.rel_position}'
        if self.rel_position_orientation is not None:
            s += f',\n\t\tposition orientation, {self.rel_position_orientation}'
        if self.rel_rotation_orientation is not None:
            s += f',\n\t\trotation orientation, {self.rel_rotation_orientation}'
        if self.abs_position is not None:
            s += f',\n\t# GROUND'
            s += f',\n\tposition, {self.abs_position}'
        if self.abs_position_orientation is not None:
            s += f',\n\tposition orientation, {self.abs_position_orientation}'
        if self.abs_rotation_orientation is not None:
            s += f',\n\trotation orientation, {self.abs_rotation_orientation}'
        if self.position_status is not None:
            s += f',\n\tposition constraint, '
            s += ', '.join(str(ps) for ps in self.position_status)
            if self.position_drive is not None:
                s += f',\n\t\t'
                s += ', '.join(str(i) for i in self.position_drive)
        if self.orientation_status is not None:
            s += f',\n\torientation constraint, '
            s += ', '.join(str(os) for os in self.orientation_status)
            if self.orientation_drive is not None:
                s += f',\n\t\t'
                s += ', '.join(str(i) for i in self.orientation_drive)
        s += self.element_footer()
        return s

class JointRegularization(Element):
    coefficients: Union[List[float], List[MBVar], float, MBVar]

    def element_type(self):
        return 'joint regularization'

    def __str__(self):
        s = f'{self.element_header()}, tikhonov'
        if isinstance(self.coefficients, list):
            s += f',\n\tlist, '
            s += ', '.join(str(co) for co in self.coefficients)
        else:
            s += f',\n\t{self.coefficients}'
        s += self.element_footer()
        return s
    
class DeformableDisplacement(Element):
    node_1: Node
    position_1: Position
    orientation_mat_1: Optional[Position] = None
    node_2: Node
    position_2: Position
    orientation_mat_2: Optional[Position] = None
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']

    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 3D constitutive law with law_type 'D3_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")
                
    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, deformable displacement' # According to the manual, 'deformable displacement joint' but in the previous implementation it was given as 'deformable displacement'
        s += f',\n\t{self.node_1.idx}'
        s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2.idx}'
        s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += f',\n\t{self.const_law}'
        s += self.element_footer()
        return s

class DeformableJoint(Element):
    node_1: Node
    position_1: Position
    orientation_mat_1: Optional[Position] = None
    node_2: Node
    position_2: Position
    orientation_mat_2: Optional[Position] = None
    const_law: Union['ConstitutiveLaw', 'NamedConstitutiveLaw']
    orientation_desc: Optional[Literal['euler123', 'euler313', 'euler321', 'orientation vector', 'orientation matrix']] = None

    @field_validator('const_law')
    def validate_const_law(cls, v):
        if isinstance(v, ConstitutiveLaw):
            if v.law_type != ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW:
                raise ValueError("const_law must be a 6D constitutive law with law_type 'D6_ISOTROPIC_LAW'")
            return v
        elif isinstance(v, NamedConstitutiveLaw):
            return v
        else:
            raise TypeError("const_law must be an instance of ConstitutiveLaw or NamedConstitutiveLaw")
                
    def element_type(self):
        return 'joint'
    
    def __str__(self):
        s = f'{self.element_header()}, deformable joint'
        s += f',\n\t{self.node_1.idx}'
        s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2.idx}'
        s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += f',\n\t{self.const_law}'
        if self.orientation_desc is not None:
            s += f',\n\t{self.orientation_desc}'
        s += self.element_footer()
        return s
        
class SphericalHinge(Element):
    '''
    This joint constrains the relative position of two nodes; the relative orientation is not constrained.
    '''

    node_1: Node
    position_1: Optional[Position] = None
    orientation_mat_1: Optional[Position] = None
    node_2: Node
    position_2: Optional[Position] = None
    orientation_mat_2: Optional[Position] = None

    def element_type(self):
        return 'joint'

    def __str__(self):
        s = f'{self.element_header()}, spherical hinge'
        s += f',\n\t{self.node_1.idx}'
        if self.position_1 is not None:
            s += f',\n\t\tposition, {self.position_1}'
        if self.orientation_mat_1 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_1}'
        s += f',\n\t{self.node_2.idx}'
        if self.position_2 is not None:
            s += f',\n\t\tposition, {self.position_2}'
        if self.orientation_mat_2 is not None:
            s += f',\n\t\torientation, {self.orientation_mat_2}'
        s += self.element_footer()
        return s

class Shell(Element):
    shell_type: Literal['shell4eas', 'shell4easans']
    nodes: List[Node]
    const_law_data: List

    @field_validator('const_law_data', mode='before')
    def validate_const_law(cls, v):
        if isinstance(v, list):
            return v
        return [v]
    
    @field_validator('nodes')
    def validate_nodes_count(cls, v):
        if len(v) == 4:
            return v
        raise ValueError(f'{cls.__name__}: must have 4 nodes')
    
    def element_type(self):
        return self.shell_type
    
    def __str__(self):
        s = f'{self.element_header()}'
        s += ',\n\t'
        s += ', '.join(str(i.idx) for i in self.nodes)
        s += ',\n\t'
        s += ', '.join(str(i) for i in self.const_law_data)
        s += self.element_footer()
        return s

class AerodynamicBody(Element):
    node: Node
    position: Position
    orientation: Position
    span: Union[float, MBVar]
    chord: List 
    aero_center: List
    b_c_point: List
    twist: List
    integration_points: Union[int, MBVar]
    induced_velocity: Optional[Union[int, MBVar]] = None
    tip_loss: Optional[List] = None
    control: Optional['DriveCaller'] = None
    airfoil_data: Optional[List] = []
    unsteady: Optional[Literal['bielawa']] = None
    jacobian: Optional[Union[Literal['yes', 'no'], bool]] = 'no'
    custom_output: Optional[List] = None
    
    @model_validator(mode='after')
    def validate_aerodynamic_body(self):
        # Validate span
        if not (isinstance(self.span, float) or 
                (isinstance(self.span, MBVar) and self.span.var_type in ('real', 'const real'))):
            raise ValueError(f'{self.__class__.__name__}: Surface span must be numeric or a real MBVar')
        
        # Validate integration_points
        if not ((isinstance(self.integration_points, int) and self.integration_points > 0) or 
                (isinstance(self.integration_points, MBVar) and self.integration_points.var_type in ('integer', 'const integer'))):
            raise ValueError(f'{self.__class__.__name__}: Integration points must be a positive integer or an integer MBVar')
                                
        return self
    
    def element_type(self):
        return "aerodynamic body"
    
    def __str__(self):
        s = f"{self.element_header()}\n\t{self.node.idx}"
        if self.induced_velocity:
            s += f",\n\t\tinduced velocity, {self.induced_velocity}"
        s += f",\n\t\t{self.position}"
        s += f",\n\t\t{self.orientation}"
        s += f",\n\t\t{self.span}"
        s += f",\n\t\t{', '.join(str(i) for i in self.chord)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.aero_center)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.b_c_point)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.twist)}"
        if self.tip_loss:
            s += f",\n\t\ttip loss, {', '.join(str(i) for i in self.tip_loss)}"
        s += f"\n\t\t{self.integration_points}"
        if self.control:
            s += f",\n\t\tcontrol, {self.control}"
        if self.airfoil_data:
            s += f",\n\t\t{', '.join(str(i) for i in self.airfoil_data)}"
        if self.unsteady:
            s += f",\n\t\tunsteady, {self.unsteady}"
        s += f",\n\t\tjacobian, {self.jacobian}"
        if self.output != 'yes':
            s += f",\n\toutput, {self.output}"
        if self.custom_output:
            s += f",\n\tcustom output, {', '.join(str(i) for i in self.custom_output)}"
        s += f";\n"
        return s


class AerodynamicBeam(Element):   
    beam: Beam
    positions: List[Position]
    orientations: List[Position]
    chord: List 
    aero_center: List
    b_c_point: List
    twist: List
    integration_points: Union[int, MBVar]
    induced_velocity: Optional[Union[int, MBVar]] = None
    tip_loss: Optional[List] = None
    control: Optional['DriveCaller'] = None
    airfoil_data: Optional[List] = []
    unsteady: Optional[Literal['bielawa']] = None
    jacobian: Optional[Union[Literal['yes', 'no'], bool]] = 'no'
    custom_output: Optional[List] = None # TODO: Add custom output class
    
    @model_validator(mode='after')
    def validate_aerodynamic_beam(self):
        # Validate positions and orientations
        if len(self.positions) not in {2, 3}:
            raise ValueError(f'{self.__class__.__name__}: must have 2 or 3 relative surface offsets')
        
        if len(self.orientations) not in {2, 3}:
            raise ValueError(f'{self.__class__.__name__}: must have 2 or 3 relative surface orientations')
        
        if len(self.positions) != len(self.orientations):
            raise ValueError(f'{self.__class__.__name__}: Number of positions ({len(self.positions)}) must match number of orientations ({len(self.orientations)})')
        
        # Validate integration_points
        if not ((isinstance(self.integration_points, int) and self.integration_points > 0) or 
                (isinstance(self.integration_points, MBVar) and self.integration_points.var_type in ('integer', 'const integer'))):
            raise ValueError(f'{self.__class__.__name__}: Integration points must be a positive integer or an integer MBVar')
                
        return self

    def element_type(self):
        return f"aerodynamic beam{len(self.positions)}"
    
    def __str__(self):
        s = f"{self.element_header()}\n\t {self.beam.idx}"
        if self.induced_velocity:
            s += f",\n\t\tinduced velocity {self.induced_velocity}"
        for pos, ori in zip(self.positions, self.orientations):
            s += f",\n\t\t{pos}"
            s += f",\n\t\t{ori}"
        s += f",\n\t\t{', '.join(str(i) for i in self.chord)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.aero_center)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.b_c_point)}"
        s += f",\n\t\t{', '.join(str(i) for i in self.twist)}"
        if self.tip_loss:
            s += f",\n\t\ttip loss, {', '.join(str(i) for i in self.tip_loss)}"
        s += f",\n\t\t{self.integration_points}"
        if self.control:
            s += f",\n\t\tcontrol, {', '.join(str(i) for i in self.control)}"
        if self.airfoil_data:
            s += f",\n\t\t{', '.join(str(i) for i in self.airfoil_data)}"
        if self.unsteady:
            s += f",\n\t\tunsteady, {self.unsteady[0]}"
        s += f",\n\t\tjacobian, {self.jacobian}"
        if self.output != 'yes':
            s += f",\n\toutput, {self.output}"
        if self.custom_output:
            s += f",\n\tcustom output, {', '.join(str(i) for i in self.custom_output)}"
        s += f";\n"
        return s

# General stuff
class NodeDof(MBEntity):
    """
    A node in MBDyn is an entity that owns public degrees of freedom and instantiates the corresponding
    public equations. It can lend them to other entities, called elements, to let them write contributions to
    public equations, possibly depending on the value of the public degrees of freedom.
    Usually elements access nodal degrees of freedom through well-deﬁned interfaces, at a high level. But
    in a few cases, nodal degrees of freedom must be accessed at a very low level, with the bare knowledge of
    the node label, the node type, the internal number of the degree of freedom, and the order of that degree
    of freedom (algebraic or diﬀerential, if any). The data that allows an entity to track a nodal degree of
    freedom is called NodeDof
    """

    node_label: Union[int, MBVar]    
    node_type: Literal['abstract', 'electric', 'hydraulic', 'parameter', 'structural', 'thermal']
    """refers to a non-scalar node type"""
    
    dof_number: Optional[Union[int, MBVar]] = None
    """required to indicate the requested degree of freedom"""
    
    dof_order: Optional[Literal['algebraic', 'differential']] = None

    def __str__(self) -> str:
        s = f'{self.node_label}, {self.node_type}'
        if self.dof_number is not None:
            s += f', {self.dof_number}'
        if self.dof_order is not None:
            s += f', {self.dof_order}'
        return s
    

# Drives
class DriveCaller(MBEntity):
    """
    Abstract class for C++ type `DriveCaller`. Every time some entity can be driven, i.e. a value can be expressed
    as dependent on some external input, an object of the class  `DriveCaller` is used.

    The `drive` essentially represents a scalar function, whose value can change over time or,
    through some more sophisticated means, can depend on the state of the analysis.
    Usually, the dependence over time is implicitly assumed, unless otherwise specified.

    For example, the amplitude of the force applied by a  `force` element is defined by means of a `drive`;
    as such, the value of the `drive` is implicitly calculated as a function of the time.
    However, a  `dof drive` uses a subordinate `drive` to compute its value based on the value of
    a degree of freedom of the analysis; as a consequence, the value of the `dof drive` is represented
    by the value of the subordinate `drive` when evaluated as a function of that specific degree of freedom
    at the desired time (function of function).

    The family of the `DriveCaller` object is very large.
    """

    # model_config = ConfigDict(arbitrary_types_allowed=True)

    idx: Optional[Union[MBVar, int]] = None
    """Index of this drive to reuse with references"""

    @abstractmethod
    def drive_type(self) -> str:
        """Every drive class must define this to return its MBDyn syntax name"""
        raise NotImplementedError("called drive_type of abstract DriveCaller")

    def drive_header(self) -> str:
        """common syntax for start of any drive caller"""
        # it's not just `__str__` to still require overriding it in specific drives
        if self.idx is not None and self.idx >= 0:
            # The idx possibly being None communicates the intent more clearly than checking if it's >=0
            return f'drive caller: {self.idx}, {self.drive_type()}'
        else:
            return self.drive_type()

class ArrayDriveCaller(DriveCaller):
    '''
    this is simply a front-end for the linear combination of <len(drives)> normal drives. <len(drives)> must be
    at least 1, in which case a simple drive caller is created, otherwise an array of drive callers is created and
    at every call their value is added to give the ﬁnal value of the array drive
    '''

    drives: List[DriveCaller]
    """List of drive callers to be used in the array"""

    @field_validator('drives')
    def validate_drives_not_empty(cls, v):
        """Validate that drives contains at least one drive caller"""
        if len(v) < 1:
            raise ValueError("array drive must contain at least one drive caller")
        return v

    def drive_type(self) -> str:
        return 'array'
    
    def __str__(self):
        return self._str_with_indent()
        
    def _str_with_indent(self, indent=""):
        """Helper method to handle indentation for nested arrays"""
        # Base string with appropriate indentation for current level
        s = self.drive_header()
        s += f", {len(self.drives)}"
        # Add each drive with appropriate indentation
        for drive in self.drives:
            if hasattr(drive, 'idx') and drive.idx is not None and drive.idx >= 0:
                s += f",\n{indent}\treference, {drive.idx}"
            elif isinstance(drive, ArrayDriveCaller):
                # For nested ArrayDriveCaller, increase indentation
                nested_str = drive._str_with_indent(indent + "\t")
                # Remove the header part before including
                if drive.idx is not None and drive.idx >= 0:
                    nested_str = nested_str.replace(f"drive caller: {drive.idx}, ", "")
                nested_str = nested_str.replace(f"{drive.drive_type()}", f"{indent}\t{drive.drive_type()}")
                s += f",\n{nested_str}"
            else:
                s += f",\n{indent}\t{drive}"
        return s

class BistopDriveCaller(DriveCaller):
    '''
    This drive caller returns 1.0 (TRUE) when its status is active and 0.0 (FALSE) when it is inactive.
    When in inactive status, it turns to active if the activation_condition is TRUE. When in active
    status, it turns to inactive if the deactivation_condition is TRUE.
    This drive caller is useful to implement a "robust" and irreversible status change
    '''
    
    initial_status: Optional[Literal['active', 'inactive']] = 'active'
    activation_condition: DriveCaller 
    deactivation_condition: DriveCaller
    
    def drive_type(self) -> str:
        return 'bistop'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\tinitial status, {self.initial_status},'
        s += '\n\t# activation condition drive'
        if self.activation_condition.idx is not None and self.activation_condition.idx >= 0:
            s += f'\n\treference, {self.activation_condition.idx}'
        else:
            s += f'\n\t{self.activation_condition}'
        s += '\n\t# deactivation condition drive'
        if self.deactivation_condition.idx is not None and self.deactivation_condition.idx >= 0:
            s += f'\n\treference, {self.deactivation_condition.idx}'
        else:
            s += f'\n\t{self.deactivation_condition}'
        return s

class ConstDriveCaller(DriveCaller):
    """An example of `DriveCaller` that always returns the same constant value"""

    # Note that method docstrings are inherited correctly (unlike dataclass fields)
    def drive_type(self):
        return 'const'

    const_value: Union[MBVar, float, int]
    """Value that will be output by the drive"""

    def __str__(self):
        return f'''{self.drive_header()}, {self.const_value}'''
    
class ClosestNextDriveCaller(DriveCaller):
    '''
    This drive returns a non-zero value when called for the ﬁrst time with an argument greater that or equal
    to the current threshold value, which is computed starting from initial_time and incrementing it each
    time by as many increment values as required to pass the current value of Time. As soon as a threshold
    value is exceeded, as many values of increment as required to pass the current value of Time are added,
    and the process repeats.
    This drive caller is useful within the output meter statement
    '''

    initial_time: Union[float, MBVar] = 0.
    final_time: Union[float, MBVar, Literal['forever']]
    increment: Union[DriveCaller]
    
    def drive_type(self) -> str:
        return 'closest next'

    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.initial_time}, {self.final_time}'
        s += ',\n\t# increment drive'
        if self.increment.idx is not None and self.increment.idx >= 0:
            s += f'\n\treference, {self.increment.idx}'
        else:
            s += f'\n\t{self.increment}'
        return s

class CosineDriveCaller(DriveCaller):    
    initial_time: Union[float, MBVar] = 0.0    
    angular_velocity: Union[float, MBVar]    
    amplitude: Union[float, MBVar]
    number_of_cycles: Union[float, MBVar, Literal['half', 'one', 'forever']]    
    initial_value: Union[float, MBVar] = 0.0

    def drive_type(self) -> str:
        return 'cosine'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.initial_time}, {self.angular_velocity}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value}'
        return s

class CubicDriveCaller(DriveCaller):
    const_coef: Union[MBVar, float]
    linear_coef: Union[MBVar, float]
    parabolic_coef: Union[MBVar, float]
    cubic_coef: Union[MBVar, float]
    
    def drive_type(self) -> str:
        return 'cubic'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.const_coef}, {self.linear_coef}, {self.parabolic_coef}, {self.cubic_coef}'
        return s

class DirectDriveCaller(DriveCaller):
    '''
    Transparently returns the input value; the arglist is empty. It is useful in conjunction with those drive
    callers that require their output to be fed into another drive caller, like the dof, node and element drive
    callers, when the output needs to be used as is.
    '''
    
    def drive_type(self) -> str:
        return 'direct'
    
    def __str__(self):
        return f'{self.drive_header()}'

class DiscreteFilterDriveCaller(DriveCaller):
    """
    Filters the output of the ancillary drive caller <input_drive> according to the discrete ﬁlter coeﬃcients
    """
    
    n_a: Union[int, MBVar]    
    """number of regression coeﬃcients"""

    a: List[Union[float, MBVar]]    
    """list of regression coeﬃcients"""

    b_0: Union[float, MBVar]    
    """direct transmission coeﬃcient, must always be present; set to zero if not needed"""

    n_b: Union[int, MBVar]    
    """number of input coeﬃcients"""

    b: List[Union[float, MBVar]]    
    """list of input coeﬃcients"""

    input_drive: DriveCaller
    """ancillary drive caller"""

    @field_validator('a')
    def validate_a_coefficients(cls, v, info: FieldValidationInfo):
        n_a = info.data.get('n_a', 0)
        if len(v) != n_a:
            raise ValueError(f"Length of 'a' list ({len(v)}) must match n_a ({n_a})")
        return v
    
    @field_validator('b')
    def validate_b_coefficients(cls, v, info: FieldValidationInfo):
        n_b = info.data.get('n_b', 0)
        if len(v) != n_b:
            raise ValueError(f"Length of 'b' list ({len(v)}) must match n_b ({n_b})")
        return v
    
    def drive_type(self) -> str:
        return 'discrete filter'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.n_a}'
        for a_i in self.a:
            s += f', {a_i}'
        s += f',\n\t{self.b_0}'
        s += f',\n\t{self.n_b}'
        for b_i in self.b:
            s += f', {b_i}'
        if self.input_drive.idx is not None and self.input_drive.idx >= 0:
            s += f',\n\treference, {self.input_drive.idx}'
        else:
            s += f',\n\t{self.input_drive}'
        return s

class DofDriveCaller(DriveCaller):
    '''
    a NodeDof, namely the reference to a degree of freedom of a node, is read. Then a recursive call to a
    drive data is read. The driver returns the value of the <func_drive> using the value of the NodeDof
    as input instead of the time. This can be used as a sort of explicit feedback, to implement fancy springs
    (where a force is driven through a function by the displacement of the node it is applied to) or an active
    control system
    '''
    
    driving_dof: NodeDof    
    func_drive: DriveCaller
    
    def drive_type(self) -> str:
        return 'dof'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.driving_dof}'
        if self.func_drive.idx is not None and self.func_drive.idx >= 0:
            s += f',\n\treference, {self.func_drive.idx}'
        else:
            s += f',\n\t{self.func_drive}'
        return s

class DoubleRampDriveCaller(DriveCaller):
    a_slope: Union[float, MBVar]    
    a_initial_time: Union[float, MBVar]
    a_final_time: Union[float, MBVar]    
    d_slope: Union[float, MBVar]    
    d_initial_time: Union[float, MBVar]    
    d_final_time: Union[float, MBVar, Literal['forever']]    
    initial_value: Union[float, MBVar]
    
    def drive_type(self) -> str:
        return 'double ramp'
    
    def __init__(self, **kwargs):
        # Check if a_initial_time wasn't explicitly provided
        if 'a_initial_time' not in kwargs:
            warnings.warn(
                "DoubleRampDriveCaller: <a_initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['a_initial_time'] = 0.0
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.a_slope}, {self.a_initial_time}, {self.a_final_time}'
        s += f',\n\t{self.d_slope}, {self.d_initial_time}, {self.d_final_time}'
        s += f',\n\t{self.initial_value}'
        return s

class DoubleStepDriveCaller(DriveCaller):
    initial_time: Union[float, MBVar]
    final_time: Union[float, MBVar]
    step_value: Union[float, MBVar]
    initial_value: Union[float, MBVar]

    def drive_type(self) -> str:
        return 'double step'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f',\n\t{self.initial_time}, {self.final_time}'
        s += f',\n\t{self.step_value}, {self.initial_value}'
        return s

class DriveDriveCaller(DriveCaller):    
    drive_caller1: DriveCaller
    drive_caller2: DriveCaller
    
    def drive_type(self) -> str:
        return 'drive'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        if self.drive_caller1.idx is not None and self.drive_caller1.idx >= 0:
            s += f',\n\treference, {self.drive_caller1.idx}'
        else:
            s += f',\n\t{self.drive_caller1}'
        if self.drive_caller2.idx is not None and self.drive_caller2.idx >= 0:
            s += f',\n\treference, {self.drive_caller2.idx}'
        else:
            s += f',\n\t{self.drive_caller2}'
        return s

class ElementDriveCaller(DriveCaller):    
    element: Element    
    private_data: str    
    func_drive: Union[DriveCaller, Literal['direct']]
    
    def drive_type(self) -> str:
        return 'element'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.element.idx}, {self.element.element_type()}'
        s += f', string, "{self.private_data}"'
        if isinstance(self.func_drive, str) and self.func_drive == 'direct':
            s += f', {self.func_drive}'
        elif hasattr(self.func_drive, 'idx') and self.func_drive.idx is not None and self.func_drive.idx >= 0:
            s += f', reference, {self.func_drive.idx}'
        else:
            s += f', {self.func_drive}'
        return s

class ExponentialDriveCaller(DriveCaller):
    """
    This drive yields a function that resembles the response of a ﬁrst-order system to a step input. Its
    value corresponds to initial_value for t < initial_time. For t ≥ initial_time, it grows to
    initial_value+amplitude_value exponentially. The growth rate is governed by time_constant_value
    """
    
    amplitude_value: Union[float, MBVar]    
    time_constant_value: Union[float, MBVar]    
    initial_time: Union[float, MBVar]    
    initial_value: Union[float, MBVar]

    def drive_type(self) -> str:
        return 'exponential'

    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
        
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.amplitude_value}, {self.time_constant_value}, {self.initial_time}, {self.initial_value}'
        return s

class FileDriveDrive(DriveCaller):
    # TODO: needs FileDrive before
    pass

class FourierSeriesDriveCaller(DriveCaller):
    """
    This drive corresponds to a Fourier series of fundamental angular velocity ω, truncated after n terms,
    over a given number of cycles P and starting at a given initial time

    f(t) = a_0/2 + ∑(a_k*cos(kω(t-t0)) + b_k*sin(kω(t-t0))) for k=1 to n
    """
    
    initial_time: Union[float, MBVar]    
    angular_velocity: Union[float, MBVar]    
    number_of_terms: Union[int, MBVar]    
    a_0: Union[float, MBVar]
    coefficients: List[Union[float, MBVar]]    
    number_of_cycles: Union[int, MBVar, Literal['one', 'forever']]    
    initial_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'angular_velocity', 'a_0', 'initial_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: Field must be an MBVar of type real or a float'
                f'\n-------------------\n'
            )
        return v
    
    @field_validator('number_of_terms')
    def validate_integer_mbvar(cls, v):
        if isinstance(v, MBVar) and 'integer' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: number_of_terms must be an MBVar of type integer or an int'
                f'\n-------------------\n'
            )
        return v
    
    @field_validator('coefficients')
    def validate_coefficients(cls, v, info: FieldValidationInfo):
        number_of_terms = info.data.get('number_of_terms')
        if isinstance(number_of_terms, MBVar):
            expected_length = 2 * number_of_terms.expression
        else:
            expected_length = 2 * number_of_terms
        if len(v) != expected_length:
            raise ValueError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: coefficients list should have {expected_length} elements '
                f'(a_1, b_1, a_2, b_2, ..., a_n, b_n) for {number_of_terms} terms'
                f'\n-------------------\n'
            )
        
        # Validate all coefficients are real numbers or MBVars of type real
        for i, coef in enumerate(v):
            if isinstance(coef, MBVar) and 'real' not in coef.var_type:
                raise TypeError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: Coefficient at index {i} must be an MBVar of type real or a float'
                    f'\n-------------------\n'
                )
        return v
    
    def drive_type(self) -> str:
        return 'fourier series'

    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
        
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}, {self.angular_velocity}, {self.number_of_terms}'
        # Add a_0 term
        s += f',\n\t{self.a_0}'     
        if isinstance(self.number_of_terms, MBVar):
            expected_length = 2 * self.number_of_terms.expression
        else:
            expected_length = 2 * self.number_of_terms   
        for i in range(0, expected_length, 2):
            if i+1 < expected_length:
                # Add a_k, b_k pair
                s += f',\n\t{self.coefficients[i]}, {self.coefficients[i+1]}'
        s += f',\n\t{self.number_of_cycles}, {self.initial_value}'
        return s
    
class FrequencySweepDriveCaller(DriveCaller):
    """
    this drive recursively calls two other drives that supply the angular velocity and the amplitude of the
    oscillation
    """
    
    initial_time: Union[float, MBVar]    
    angular_velocity_drive: DriveCaller    
    amplitude_drive: DriveCaller    
    initial_value: Union[float, MBVar]    
    final_time: Union[float, MBVar, Literal['forever']]    
    final_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'initial_value', 'final_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: Field must be an MBVar of type real or a float'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'frequency sweep'

    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
        
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}'
        if hasattr(self.angular_velocity_drive, 'idx') and self.angular_velocity_drive.idx is not None and self.angular_velocity_drive.idx >= 0:
            s += f',\n\treference, {self.angular_velocity_drive.idx}'
        else:
            s += f',\n\t{self.angular_velocity_drive}'
        if hasattr(self.amplitude_drive, 'idx') and self.amplitude_drive.idx is not None and self.amplitude_drive.idx >= 0:
            s += f',\n\treference, {self.amplitude_drive.idx}'
        else:
            s += f',\n\t{self.amplitude_drive}'
        s += f'\n\t{self.initial_value}, {self.final_time}, {self.final_value}'
        return s
    
class GiNaCDriveCaller(DriveCaller):
    """
    The GiNaC drive caller evaluates mathematical expressions.
    The function expression is evaluated and differentiated, if needed, as a function of the variable 
    passed as the optional symbol. If none is passed, 'Var' is used.
    """
    
    expression: Union[str, MBVar]    
    symbol: Optional[Union[str, MBVar]] = None
    
    @field_validator('expression')
    def validate_expression(cls, v):
        if isinstance(v, MBVar) and 'string' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'GiNaCDriveCaller: expression must be a string or an MBVar of type string'
                '\n-------------------\n'
            )
        return v
    
    @field_validator('symbol')
    def validate_symbol(cls, v):
        if v is not None and isinstance(v, MBVar) and 'string' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'GiNaCDriveCaller: symbol must be a string or an MBVar of type string'
                '\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'ginac'
    
    def __str__(self):
        s = self.drive_header()
        if self.symbol is not None:
            if isinstance(self.symbol, MBVar):
                s += f', symbol, {self.symbol}'
            else:
                s += f', symbol, "{self.symbol}"'
        if isinstance(self.expression, MBVar):
            s += f', {self.expression}'
        else:
            s += f', "{self.expression}"'            
        return s

class LinearDriveCaller(DriveCaller):
    """
    The Linear drive caller implements a linear function of time:
    f(t) = const_coef + slope_coef · t
    """
    
    const_coef: Union[float, MBVar]    
    slope_coef: Union[float, MBVar]
    
    @field_validator('const_coef', 'slope_coef')
    def validate_coefficients(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'LinearDriveCaller: coefficients must be real numbers or MBVars of type real'
                '\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'linear'
    
    def __str__(self):
        s = self.drive_header()
        s += f', {self.const_coef}, {self.slope_coef}'
        return s

class MeterDriveCaller(DriveCaller):
    """
    The Meter drive caller has value zero except for every 'steps_between_spikes' steps, 
    where it assumes unit value.
    """
    
    initial_time: Union[float, MBVar]    
    final_time: Union[float, MBVar, Literal['forever']]    
    steps_between_spikes: Optional[Union[int, MBVar]] = None
    
    @field_validator('initial_time', 'final_time')
    def validate_initial_time(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: Field must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    @field_validator('steps_between_spikes')
    def validate_steps(cls, v):
        if v is not None:
            if isinstance(v, MBVar) and 'integer' not in v.var_type:
                raise TypeError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: steps_between_spikes must be an integer or an MBVar of type integer'
                    f'\n-------------------\n'
                )
            elif isinstance(v, int) and v <= 0:
                raise ValueError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: steps_between_spikes must be a positive integer'
                    f'\n-------------------\n'
                )
        return v
    
    def drive_type(self) -> str:
        return 'meter'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        super().__init__(**kwargs)
    
    def __str__(self):
        s = self.drive_header()
        s += f', {self.initial_time}, {self.final_time}'
        if self.steps_between_spikes is not None:
            s += f', steps, {self.steps_between_spikes}'
        return s
            
class MultDriveCaller(DriveCaller):
    """
    The Mult drive caller multiplies the value of two subordinate drives.
    """
    
    drive_1: DriveCaller    
    drive_2: DriveCaller
        
    def drive_type(self) -> str:
        return 'mult'
    
    def __str__(self):
        s = self.drive_header()
        if hasattr(self.drive_1, 'idx') and self.drive_1.idx is not None and self.drive_1.idx >= 0:
            s += f',\n\treference, {self.drive_1.idx}'
        else:
            s += f',\n\t{self.drive_1}'
        if hasattr(self.drive_2, 'idx') and self.drive_2.idx is not None and self.drive_2.idx >= 0:
            s += f',\n\treference, {self.drive_2.idx}'
        else:
            s += f',\n\t{self.drive_2}'
        return s
    
class NodeDriveCaller(DriveCaller):
    """    
    The driver returns the value of the func_drive using the value of the node's private data
    as input instead of the time. This can be used as a sort of explicit feedback, to implement
    fancy springs (where a force is driven through a function by the rotation of a joint) or an
    active control system.
    """
    
    node: Node        
    private_data: str    
    func_drive: Union[DriveCaller, Literal['direct']]
    
    def drive_type(self) -> str:
        return 'node'
    
    def __str__(self):
        s = f'{self.drive_header()}'        
        s += f', {self.node.idx}, {self.node.node_type}'
        s += f', string, "{self.private_data}"'
        if hasattr(self.func_drive, 'idx') and self.func_drive.idx is not None and self.func_drive.idx >= 0:
            s += f', reference, {self.func_drive.idx}'
        else:
            s += f', {self.func_drive}'
        return s
    
class NullDriveCaller(DriveCaller):
    """Zero valued drive caller; the arglist is empty."""
    
    def drive_type(self) -> str:
        return 'null'
    
    def __str__(self):
        return f'{self.drive_header()}'
    
class ParabolicDriveCaller(DriveCaller):
    """
    The Parabolic drive caller implements a quadratic function of time:
    f(t) = const_coef + linear_coef · t + parabolic_coef · t²
    """
    
    const_coef: Union[float, MBVar]    
    linear_coef: Union[float, MBVar]    
    parabolic_coef: Union[float, MBVar]
    
    @field_validator('const_coef', 'linear_coef', 'parabolic_coef')
    def validate_coefficients(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'ParabolicDriveCaller: coefficients must be real numbers or MBVars of type real'
                '\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'parabolic'
    
    def __str__(self):
        s = self.drive_header()
        s += f', {self.const_coef}, {self.linear_coef}, {self.parabolic_coef}'
        return s
    
class PeriodicDriveCaller(DriveCaller):
    """
    Represents a periodic drive function that is zero before the initial time and follows 
    f(t) = func_drive(t - initial_time - period * floor((t - initial_time) / period)) for t ≥ initial_time.
    """

    initial_time: Union[float, MBVar]
    period: Union[float, MBVar]
    func_drive: DriveCaller
    
    @field_validator('initial_time', 'period')
    def validate_time_params(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: time parameters must be real numbers or MBVars of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'periodic'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0   
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}, {self.period}'
        if hasattr(self.func_drive, 'idx') and self.func_drive.idx is not None and self.func_drive.idx >= 0:
            s += f',\n\treference, {self.func_drive.idx}'
        else:
            s += f',\n\t{self.func_drive}'
        return s
    
class PiecewiseLinearDriveCaller(DriveCaller):
    """    
    The function performs linear interpolation between defined (point, value) pairs.
    The first and last point/value pairs are extrapolated if a value beyond the extremes is required.
    """
    
    num_points: Union[int, MBVar]    
    points_values: List[Tuple[Union[float, MBVar], Union[float, MBVar]]]
    """List of (point, value) coordinate pairs defining the piecewise linear function"""

    @field_validator('num_points')
    def validate_num_points(cls, v):
        if isinstance(v, MBVar) and 'integer' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'PiecewiseLinearDriveCaller: num_points must be an integer or an MBVar of type integer'
                '\n-------------------\n'
            )
        if v < 2:
            raise ValueError(
                '\n-------------------\nERROR: '
                'PiecewiseLinearDriveCaller: num_points must be at least 2'
                '\n-------------------\n'
            )
        return v

    @model_validator(mode='after')
    def validate_points_values_length(self):
        """Validate that the number of points matches the length of points_values"""
        if isinstance(self.num_points, int) and len(self.points_values) != self.num_points:
            raise ValueError(
                '\n-------------------\nERROR: '
                f'PiecewiseLinearDriveCaller: number of (point, value) pairs ({len(self.points_values)}) '
                f'does not match num_points ({self.num_points})'
                '\n-------------------\n'
            )
        return self
    
    def drive_type(self) -> str:
        return 'piecewise linear'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.num_points}'
        for point, value in self.points_values:
            s += f',\n\t{point}, {value}'
        return s
    
class PostponedDriveCaller(DriveCaller):
    """    
    This drive is a stub for a drive that cannot be defined early in the input file 
    because it occurs when the data manager is not yet available.
    A drive caller with the same label must be defined before this drive caller is first used.
    """
    
    label: Union[int, MBVar]
    """Label that identifies this drive for later definition"""
        
    def drive_type(self) -> str:
        return 'postponed'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.label}'
        return s
        
class RampDriveCaller(DriveCaller):
    """
    The Ramp drive caller implements a ramp function with specified slope.
    
    f(t) = initial_value                           if t < initial_time
           initial_value + slope·(t-initial_time)  if initial_time ≤ t ≤ final_time
           initial_value + slope·(final_time-initial_time)  if t > final_time
    """
    
    slope: Union[float, MBVar]    
    initial_time: Union[float, MBVar]  
    final_time: Union[float, MBVar, Literal['forever']]    
    initial_value: Union[float, MBVar]
    
    @field_validator('slope', 'initial_time', 'final_time', 'initial_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <{v}> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'ramp'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)

    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.slope}, {self.initial_time}, {self.final_time}, {self.initial_value}'
        return s

class RandomDriveCaller(DriveCaller):
    """
    The Random drive caller generates pseudo-random numbers.
    Numbers are uniformly distributed in the interval [mean_value - amplitude_value, mean_value + amplitude_value).
    The output can be held for a specified number of steps, and the random seed can be specified.
    """
    
    amplitude_value: Union[float, MBVar]    
    mean_value: Union[float, MBVar]    
    initial_time: Union[float, MBVar]  
    final_time: Union[float, MBVar, Literal['forever']]    
    steps_to_hold_value: Optional[Union[int, MBVar]] = None    
    seed_value: Optional[Union[int, MBVar, Literal['time']]] = None
    
    @field_validator('amplitude_value', 'mean_value', 'initial_time', 'final_time')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <{v}> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    @field_validator('steps_to_hold_value')
    def validate_steps(cls, v):
        if v is not None:
            if isinstance(v, MBVar) and 'integer' not in v.var_type:
                raise TypeError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: <steps_to_hold_value> must be an integer or an MBVar of type integer'
                    f'\n-------------------\n'
                )
            if v <= 0:
                raise ValueError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: <steps_to_hold_value> must be a positive integer'
                    f'\n-------------------\n'
                )
        return v
    
    @field_validator('seed_value')
    def validate_seed(cls, v):
        if v is not None and v != 'time':
            if isinstance(v, MBVar) and 'integer' not in v.var_type:
                raise TypeError(
                    f'\n-------------------\nERROR: '
                    f'{cls.__name__}: <seed_value> must be an integer, "time", or an MBVar of type integer'
                    f'\n-------------------\n'
                )
        return v
    
    def drive_type(self) -> str:
        return 'random'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.amplitude_value}, {self.mean_value}, {self.initial_time}, {self.final_time}'
        if self.steps_to_hold_value is not None:
            s += f', steps, {self.steps_to_hold_value}'
        if self.seed_value is not None:
            s += f', seed, {self.seed_value}'
        return s
        
class SampleAndHoldDriveCaller(DriveCaller):
    """    
    When trigger is non-zero, the value of function is recorded after convergence at the end of the time
    step, and returned whenever the drive is called afterwards. When trigger is zero, the last recorded value
    is returned. When initial_value is provided, if the trigger is initially zero, this value is returned 
    until trigger becomes non-zero.
    """
    
    function: DriveCaller    
    trigger: DriveCaller    
    initial_value: Optional[Union[float, MBVar]] = None
    
    @field_validator('initial_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <initial_value> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'sample and hold'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        if hasattr(self.function, 'idx') and self.function.idx is not None and self.function.idx >= 0:
            s += f',\n\treference, {self.function.idx}'
        else:
            s += f',\n\t{self.function}'
        if hasattr(self.trigger, 'idx') and self.trigger.idx is not None and self.trigger.idx >= 0:
            s += f',\n\treference, {self.trigger.idx}'
        else:
            s += f',\n\t{self.trigger}'
        if self.initial_value is not None:
            s += f', initial value, {self.initial_value}'
        return s

class ScalarFunctionDriveCaller(DriveCaller):
    pass
        
class SineDriveCaller(DriveCaller):
    """
    The Sine drive caller implements a sinusoidal function:
    f(t) = initial_value + amplitude · sin(angular_velocity · (t - initial_time))
    """
    
    initial_time: Union[float, MBVar]  
    angular_velocity: Union[float, MBVar]    
    amplitude: Union[float, MBVar]    
    number_of_cycles: Union[int, MBVar, Literal['half', 'one', 'forever']]    
    initial_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'initial_value', 'angular_velocity', 'amplitude')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: Field must be an MBVar of type real or a float'
                f'\n-------------------\n'
            )
        return v
    
    @field_validator('number_of_cycles')
    def validate_cycles(cls, v):
        if isinstance(v, str) and v not in ['half', 'one', 'forever']:
            raise ValueError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: string value for number_of_cycles must be one of: '
                f'"half", "one", "forever"'
                f'\n-------------------\n'
            )
        elif isinstance(v, MBVar) and 'integer' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: number_of_cycles must be an integer, '
                f'one of ("half", "one", "forever"), or an MBVar of type integer'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'sine'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)

    def __str__(self):
        s = self.drive_header()
        s += f', {self.initial_time}, {self.angular_velocity}, {self.amplitude}'
        s += f', {self.number_of_cycles}, {self.initial_value}'
        return s

class StepDriveCaller(DriveCaller):
    """    
    f(t) = 0              if t < initial_time
           step_value     if t >= initial_time
    """
    
    initial_time: Union[float, MBVar]
    step_value: Union[float, MBVar]  
    initial_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'step_value', 'initial_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <{v}> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'step'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}, {self.step_value}, {self.initial_value}'
        return s

class Step5DriveCaller(DriveCaller):
    initial_time: Union[float, MBVar]
    initial_value: Union[float, MBVar]
    final_time: Union[float, MBVar]
    final_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'initial_value', 'final_time', 'final_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <{v}> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'step5'

    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning 
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}, {self.initial_value}, {self.final_time}, {self.final_value}'
        return s

class StringDriveCaller(DriveCaller):    
    expression: Union[str, MBVar]
    
    @field_validator('expression')
    def validate_expression(cls, v):
        if isinstance(v, MBVar) and 'string' not in v.var_type:
            raise TypeError(
                '\n-------------------\nERROR: '
                'StringDriveCaller: expression must be a string or an MBVar of type string'
                '\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'string'
    
    def __str__(self):
        s = self.drive_header()
        if isinstance(self.expression, MBVar):
            s += f', "{self.expression.expression}"'
        else:
            s += f', "{self.expression}"'
        return s
    
class TanhDriveCaller(DriveCaller):
    """    
    f(t) = initial_value + amplitude · tanh(nd_slope · (t - initial_time))
    """
    
    initial_time: Union[float, MBVar]
    amplitude: Union[float, MBVar]  
    nd_slope: Union[float, MBVar]
    initial_value: Union[float, MBVar]
    
    @field_validator('initial_time', 'amplitude', 'nd_slope', 'initial_value')
    def validate_real_mbvar(cls, v):
        if isinstance(v, MBVar) and 'real' not in v.var_type:
            raise TypeError(
                f'\n-------------------\nERROR: '
                f'{cls.__name__}: <{v}> must be a real number or an MBVar of type real'
                f'\n-------------------\n'
            )
        return v
    
    def drive_type(self) -> str:
        return 'tanh'
    
    def __init__(self, **kwargs):
        # Check if initial_time wasn't explicitly provided
        if 'initial_time' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_time> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_time'] = 0.0
            
        # Check if initial_value wasn't explicitly provided
        if 'initial_value' not in kwargs:
            warnings.warn(
                f"{self.__class__.__name__}: <initial_value> is not set, assuming 0.0.",
                UserWarning
            )
            kwargs['initial_value'] = 0.0
            
        super().__init__(**kwargs)
    
    def __str__(self):
        s = f'{self.drive_header()}'
        s += f', {self.initial_time}, {self.amplitude}, {self.nd_slope}, {self.initial_value}'
        return s

class TimeDriveCaller(DriveCaller):
    """
    Yields the current time.
    """
    
    def drive_type(self) -> str:
        return 'time'
    
    def __str__(self):
        s = f'{self.drive_header()}'
        return s

class TimestepDriveCaller(DriveCaller):
    """
    Yields the current timestep.
    """
    
    def drive_type(self) -> str:
        return 'timestep'
    
    def __str__(self):
        return f'{self.drive_header()}'
    
class UnitDriveCaller(DriveCaller):
    """Always 1"""
    
    def drive_type(self) -> str:
        return 'unit'
    
    def __str__(self):
        return f'{self.drive_header()}'
                
class TplDriveCaller(DriveCaller):
    pass

if imported_pydantic:
    DriveDisplacement.model_rebuild()
    DriveDisplacementPin.model_rebuild()
    DriveHinge.model_rebuild()

    AngularAcceleration.model_rebuild()
    AngularVelocity.model_rebuild()
    AxialRotation.model_rebuild()
    Brake.model_rebuild()
    ImposedDisplacement.model_rebuild()
    ImposedDisplacement.model_rebuild()

class ConstitutiveLaw(MBEntity):
    """
    Abstract class for C++ type `ConstitutiveLaw`. Every time a "deformable"
    entity requires a constitutive law, a template constitutive law is read. This has been implemented by
    means of C++ templates in order to allow the definition of a general constitutive law when possible.

    Constitutive laws are also used in non-structural components, to allow some degree of generality in
    defining input/output relationships. Some constitutive laws are meaningful only when related to some
    precise dimensionality. In some special cases, general purpose elements use 1D constitutive laws
    to express an arbitrary dependence of some value on a scalar state of the system.

    The meaning of the input and output parameters of a constitutive law is dictated by the entity that
    uses it. In general, the user should refer to the element the constitutive law is being instantiated for in
    order to understand what the input and the output parameters are supposed to be.
    """

    class LawType(Enum):
        SCALAR_ISOTROPIC_LAW = "scalar isotropic law"
        D3_ISOTROPIC_LAW = "3D isotropic law"
        D6_ISOTROPIC_LAW = "6D isotropic law"

    idx: Optional[Union[MBVar, int]] = None
    """Index of this constitutive law to reuse with references"""

    law_type: LawType
    prestress: Optional[List] = None
    prestrain: Optional[List] = None

    @abstractmethod
    def const_law_name(self) -> str:
        """Name of the specific constitutive law"""
        pass

    @property
    def dim(self) -> int:
        """Determine the dimensionality based on the constitutive law name"""
        if self.law_type == ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW:
            return 1
        elif self.law_type == ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW:
            return 3
        elif self.law_type == ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW:
            return 6
        else: 
            raise ValueError(f"Unknown constitutive law name: {self.law_type}") 
    
    @staticmethod
    def validate_matrix(matrix: list, name: str, supported_dims: set = {3, 6}) -> None:
        """
        Validate a matrix to ensure it's a square matrix of a supported dimension.
        Default supported dimensions are 3x3 and 6x6.
        """
        if not isinstance(matrix, list):
            raise TypeError(f"{name.capitalize()} must be a list of lists.")    # runtime guard
        
        N = len(matrix)
        if N not in supported_dims:
            raise ValueError(f"Unsupported size for {name} matrix. Expected dimensions {supported_dims}, got {N}x?.")
        
        for i, row in enumerate(matrix):
            if not isinstance(row, list):
                raise TypeError(f"Row {i} of {name} matrix must be a list.")
            if len(row) != N:
                raise ValueError(f"{name.capitalize()} matrix must be square. Expected {N}x{N}, but row {i} has length {len(row)}.")

        
    def const_law_header(self) -> str:
        """Common syntax for start of any constitutive law"""
        if self.idx is not None and self.idx >= 0:
            return f'constitutive law: {self.idx}, name, "{self.law_type.value}",' \
                   f'\n\t{self.dim}, {self.const_law_name()}'
        else:
            return self.const_law_name()
        
    def const_law_footer(self) -> str:
        """Common syntax for end of any constitutive law"""
        s = ''
        if self.prestress is not None:
            s += f',\n\tprestress, {", ".join(str(i) for i in self.prestress)}'
        if self.prestrain is not None:
            s += f',\n\tprestrain, {", ".join(str(i) for i in self.prestrain)}'
        return s
    
class LinearElastic(ConstitutiveLaw):
    def const_law_name(self) -> str:
        if self.dim == 1:
            return 'linear elastic'
        else:
            return 'linear elastic isotropic'
    
    stiffness: Union[MBVar, float]
    """The isotropic stiffness coefficient"""
    
    def __str__(self):
        s = f'{self.const_law_header()}, {self.stiffness}'
        s += self.const_law_footer()
        return s

class LinearElasticGeneric(ConstitutiveLaw):
    stiffness: Union[float, MBVar, List[List[Union[float, MBVar]]]]
    
    def const_law_name(self) -> str:
        return 'linear elastic generic'

    @model_validator(mode='before')
    @classmethod
    def validate_stiffness_matrix(cls, values):
        stiffness = values.get('stiffness')
        # Only validate if it's a list (matrix)
        if isinstance(stiffness, list):
            cls.validate_matrix(stiffness, 'stiffness', supported_dims={1, 3, 6})
        return values

    def __str__(self):
        s = self.const_law_header()
        if isinstance(self.stiffness, (float, MBVar)):
            s += f', {self.stiffness}'
        elif isinstance(self.stiffness, list):
            matrix_str = ''
            for row in self.stiffness:
                row_str = ', '.join(map(str, row))
                matrix_str += f',\n\t{row_str}'
            s += matrix_str
        else:
            raise TypeError(f"{self.__class__.__name__}: Invalid type for stiffness matrix")
        s += self.const_law_footer()
        return s

class LinearElasticGenericAxialTorsionCoupling(ConstitutiveLaw):
    stiffness: List[List[Union[float, MBVar]]]
    coupling_coef: Union[float, MBVar]

    def const_law_name(self) -> str:
        return 'linear elastic generic axial torsion coupling'

    @model_validator(mode='before')
    @classmethod
    def validate_stiffness_matrix(cls, values):
        stiffness = values.get('stiffness')
        if isinstance(stiffness, list):
            cls.validate_matrix(stiffness, 'stiffness', supported_dims={6})
        return values

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        matrix_str = ''
        for row in self.stiffness:
            row_str = ', '.join(map(str, row))
            matrix_str += f',\n\t{row_str}'
        base_str += matrix_str
        base_str += f',\n\t{self.coupling_coef}'
        base_str += self.const_law_footer()
        return base_str
    
class CubicElasticGeneric(ConstitutiveLaw):
    stiffness_1: Union[float, MBVar, List[Union[float, MBVar]]]
    stiffness_2: Union[float, MBVar, List[Union[float, MBVar]]]
    stiffness_3: Union[float, MBVar, List[Union[float, MBVar]]]

    @model_validator(mode='before')
    @classmethod
    def validate_stiffness_forms(cls, values):
        s1, s2, s3 = values.get('stiffness_1'), values.get('stiffness_2'), values.get('stiffness_3')
        is_s1_list = isinstance(s1, list)
        is_s2_list = isinstance(s2, list)
        is_s3_list = isinstance(s3, list)
        # Ensure all three stiffnesses are of the same type (all scalars or all lists)
        if not (is_s1_list == is_s2_list == is_s3_list):
            raise TypeError("stiffness_1, stiffness_2, and stiffness_3 must all be scalars (for 1D) or all lists (for 3D).")
        # If they are lists (3D vector case), validate that they are all vectors of length 3
        if is_s1_list:
            if not (len(s1) == 3 and len(s2) == 3 and len(s3) == 3):
                raise ValueError("For 3D vector form, stiffness_1, stiffness_2, and stiffness_3 must each be a list of 3 elements.")
        return values

    def const_law_name(self) -> str:
        return 'cubic elastic generic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        if isinstance(self.stiffness_1, (float, MBVar)):
            # Scalar case
            base_str += f', {self.stiffness_1}, {self.stiffness_2}, {self.stiffness_3}'
        else:
            # 3x1 Vector case
            s1_str = ', '.join(map(str, self.stiffness_1))
            s2_str = ', '.join(map(str, self.stiffness_2))
            s3_str = ', '.join(map(str, self.stiffness_3))
            base_str += f',\n\t{s1_str},\n\t{s2_str},\n\t{s3_str}'
        base_str += self.const_law_footer()
        return base_str

class InverseSquareElastic(ConstitutiveLaw):
    """
    Inverse square elastic constitutive law
    """
    
    stiffness: Union[MBVar, float]
    ref_length: Union[MBVar, float]
    
    def const_law_name(self) -> str:
        return 'inverse square elastic'
    
    def __str__(self):
        s = f'{self.const_law_header()}, {self.stiffness}, {self.ref_length}'
        s += self.const_law_footer()
        return s
    
class LogElastic(ConstitutiveLaw):
    """
    Logarithmic elastic constitutive law
    """

    stiffness: Union[float, MBVar]

    def const_law_name(self) -> str:
        return 'log elastic'

    def __str__(self):
        s = f'{self.const_law_header()}, {self.stiffness}'
        s += self.const_law_footer()
        return s

class LinearElasticBistop(ConstitutiveLaw):
    """
    Linear elastic bistop constitutive law
    """

    stiffness: Union[float, MBVar]
    initial_status: Optional[Union[bool, str]]
    activating_condition: DriveCaller
    deactivating_condition: DriveCaller

    def const_law_name(self) -> str:
        return 'linear elastic bistop'

    # TODO: Ensure the string representation is correct
    def __str__(self):
        base_str = f'{self.const_law_header()}, {self.stiffness}'
        
        if self.initial_status is not None:
            base_str += f',\n\tinitial status, {self.initial_status},'

        base_str += '\n\t# activation condition drive'
        if self.activation_condition.idx is None:
            base_str += f'\n\t{self.activation_condition},'
        else:
            base_str += f'\n\treference, {self.activation_condition.idx},'
        base_str += '\n\t# activation condition drive'
        if self.deactivating_condition.idx is None:
            base_str += f'\n\t{self.deactivating_condition},'
        else:
            base_str += f'\n\treference, {self.deactivating_condition.idx}'
        base_str += self.const_law_footer()
        return base_str

class DoubleLinearElastic(ConstitutiveLaw):
    """
    Double linear elastic constitutive law
    """

    stiffness_1: Union[MBVar, float]
    upper_strain: Union[MBVar, float]
    lower_strain: Union[MBVar, float]
    stiffness_2: Union[MBVar, float]
    
    def const_law_name(self) -> str:
        return 'double linear elastic'
    
    def __str__(self):
        s = f'{self.const_law_header()}, {self.stiffness_1}, {self.upper_strain}, {self.lower_strain}, {self.stiffness_2}'
        s += self.const_law_footer()
        return s

class IsotropicHardeningElastic(ConstitutiveLaw):
    """
    Isotropic hardening elastic constitutive law
    """

    stiffness: Union[MBVar, float]
    reference_strain: Union[MBVar, float]
    linear_stiffness: Optional[Union[MBVar, float]]
    
    def const_law_name(self) -> str:
        return 'isotropic hardening elastic'
    
    def __str__(self):
        if self.linear_stiffness is not None:
            s = f'{self.const_law_header()}, {self.stiffness}, {self.reference_strain}, linear stiffness, {self.linear_stiffness}'
            s += self.const_law_footer()
            return s
        else:
            s = f'{self.const_law_header()}, {self.stiffness}, {self.reference_strain}'
            s += self.const_law_footer()
            return s
        
class LinearViscous(ConstitutiveLaw):
    """
    Linear viscous constitutive law
    """

    viscosity: Union[MBVar, float]    
    
    def const_law_name(self) -> str:
        if self.dim == 1:
            return 'linear viscous'
        else:
            return 'linear viscous isotropic'

    def __str__(self):
        s = f'{self.const_law_header()}, {self.viscosity}'
        s += self.const_law_footer()
        return s

class LinearViscousGeneric(ConstitutiveLaw):
    """
    Linear viscous generic constitutive law
    """

    viscosity: Union[float, MBVar, List[List[Union[float, MBVar]]]]

    def const_law_name(self) -> str:
        return 'linear viscous generic'

    def __str__(self):
        if isinstance(self.viscosity, (float, MBVar)):
            s = f'{self.const_law_header()}, {self.viscosity}'
            s += self.const_law_footer()
            return s
        elif isinstance(self.viscosity, list):
            N = len(self.viscosity)
            if N == 1:
                s = f'{self.const_law_header()}, {self.viscosity[0][0]}'
                s += self.const_law_footer()
                return s
            elif N == 3 or N == 6:
                matrix_str = ''
                for i in range(N):
                    row_str = ', '.join(str(self.viscosity[i][j]) for j in range(N))
                    matrix_str += f',\n\t{row_str}'
                s = f'{self.const_law_header()}{matrix_str}'
                s += self.const_law_footer()
                return s
            else:
                raise ValueError(f"{self.__class__.__name__}: Unsupported size of viscosity matrix")
        else:
            raise TypeError(f"{self.__class__.__name__}: Invalid type for viscosity matrix")
 
class LinearViscoelastic(ConstitutiveLaw):
    """
    Linear viscoelastic constitutive law
    """

    stiffness: Union[MBVar, float]

    viscosity: Union[MBVar, float]
    """The viscosity coefficient"""

    factor: Optional[Union[MBVar, float]]
    """Factor for proportional viscosity"""
    
    def const_law_name(self) -> str:
        if self.dim == 1:
            return 'linear viscoelastic'
        else:
            return 'linear viscoelastic isotropic'

    def __str__(self):
        if self.viscosity is not None:
            s = f'{self.const_law_header()}, {self.stiffness}, {self.viscosity}'
            s += self.const_law_footer()
            return s
        elif self.factor is not None:
            s = f'{self.const_law_header()}, {self.stiffness}, proportional, {self.factor}'
            s += self.const_law_footer()
            return s
        else:
            raise ValueError(f"{self.__class__.__name__}: Either viscosity or factor must be provided for Linear viscoelastic law")
        
class LinearViscoelasticGeneric(ConstitutiveLaw):
    """
    Linear viscoelastic generic constitutive law
    """

    stiffness: List[List[Union[float, MBVar]]]
    viscosity: Optional[List[List[Union[float, MBVar]]]] = None
    factor: Optional[Union[float, MBVar]] = None

    @model_validator(mode='before')
    def check_viscosity_factor(cls, values: 'FieldValidationInfo') -> Any:
        stiffness = values.get('stiffness')
        viscosity = values.get('viscosity')
        factor = values.get('factor')
        # Ensure either viscosity or factor is provided, but not both
        if viscosity is not None and factor is not None:
            raise ValueError("Either viscosity or factor must be provided, not both.")
        if viscosity is None and factor is None:
            raise ValueError("One of viscosity or factor must be provided.")

        def validate_matrix(matrix: List[List[Union[float, MBVar]]], name: str) -> None:
            """Validate the matrix to ensure it's a square matrix of size 3x3 or 6x6."""
            if matrix:
                N = len(matrix)
                if N not in {3, 6}:
                    raise ValueError(f"Unsupported size of {name} matrix. Expected 3x3 or 6x6, got {N}x{N}.")
                for row in matrix:
                    if len(row) != N:
                        raise ValueError(f"{name.capitalize()} matrix must be square. Expected {N}x{N}, but found a row with length {len(row)}.")
        # Validate stiffness matrix
        if stiffness:
            validate_matrix(stiffness, 'stiffness')
        # Validate viscosity matrix
        if viscosity:
            validate_matrix(viscosity, 'viscosity')
        return values
    
    def const_law_name(self) -> str:
        return 'linear viscoelastic generic'
    
    def __str__(self):
        base_str = f'{self.const_law_header()}'
        if isinstance(self.stiffness, (float, MBVar)):
            base_str += f', {self.stiffness}'
        elif isinstance(self.stiffness, list):
            N = len(self.stiffness)
            if N == 1:
                base_str += f', {self.stiffness[0][0]}'
            elif N == 3 or N == 6:
                matrix_str = ''
                for i in range(N):
                    row_str = ', '.join(str(self.stiffness[i][j]) for j in range(N))
                    matrix_str += f',\n\t{row_str}'
                base_str += f'{matrix_str}'
            else:
                raise ValueError("Unsupported size of stiffness matrix")
        else:
            raise TypeError("Invalid type for stiffness matrix")
        if self.viscosity is not None:
            if isinstance(self.viscosity, (float, MBVar)):
                base_str += f', {self.viscosity}'
            elif isinstance(self.viscosity, list):
                N = len(self.viscosity)
                if N == 1:
                    base_str += f', {self.viscosity[0][0]}'
                elif N == 3 or N == 6:
                    matrix_str = ''
                    for i in range(N):
                        row_str = ', '.join(str(self.viscosity[i][j]) for j in range(N))
                        matrix_str += f',\n\t{row_str}'
                    base_str += f'{matrix_str}'
                else:
                    raise ValueError("Unsupported size of viscosity matrix")
            else:
                raise TypeError("Invalid type for viscosity matrix")
        elif self.factor is not None:
            base_str += f', proportional, {self.factor}'
        base_str += self.const_law_footer()
        return base_str
   
class LinearTimeVariantViscoelasticGeneric(ConstitutiveLaw):
    """
    Linear time variant viscoelastic generic constitutive law
    """

    stiffness: Union[float, MBVar, List[List[Union[float, MBVar]]]]
    stiffness_scale: DriveCaller
    viscosity: Optional[Union[float, MBVar, List[List[Union[float, MBVar]]]]] = None
    factor: Optional[Union[float, MBVar]] = None
    viscosity_scale: DriveCaller

    def const_law_name(self) -> str:
        return 'linear time variant viscoelastic generic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        
        # String representation for stiffness
        if isinstance(self.stiffness, (float, MBVar)):
            base_str += f', {self.stiffness}'
        elif isinstance(self.stiffness, list):
            N = len(self.stiffness)
            if N == 1:
                base_str += f', {self.stiffness[0][0]}'
            elif N == 3 or N == 6:
                matrix_str = ''
                for i in range(N):
                    row_str = ', '.join(str(self.stiffness[i][j]) for j in range(N))
                    matrix_str += f',\n\t{row_str}'
                base_str += f'{matrix_str}'
            else:
                raise ValueError("Unsupported size of stiffness matrix")
        else:
            raise TypeError("Invalid type for stiffness matrix")

        # String representation for stiffness scale
        if self.stiffness_scale.idx is None:
            base_str += f',\n\t{self.stiffness_scale},'
        else:
            base_str += f',\n\treference, {self.stiffness_scale.idx},'
        
        # String representation for viscosity
        if self.viscosity is not None:
            if isinstance(self.viscosity, (float, MBVar)):
                base_str += f', {self.viscosity}'
            elif isinstance(self.viscosity, list):
                N = len(self.viscosity)
                if N == 1:
                    base_str += f', {self.viscosity[0][0]}'
                elif N == 3 or N == 6:
                    matrix_str = ''
                    for i in range(N):
                        row_str = ', '.join(str(self.viscosity[i][j]) for j in range(N))
                        matrix_str += f',\n\t{row_str}'
                    base_str += f',\n{matrix_str}'
                else:
                    raise ValueError("Unsupported size of viscosity matrix")
            else:
                raise TypeError("Invalid type for viscosity matrix")
        elif self.factor is not None:
            base_str += f', proportional, {self.factor}'
        else:
            raise ValueError("Either viscosity or factor must be provided")

        # String representation for viscosity scale
        if self.viscosity_scale.idx is None:
            base_str += f',\n\t{self.viscosity_scale}'
        else:
            base_str += f',\n\treference, {self.viscosity_scale.idx}'
        base_str += self.const_law_footer()
        return base_str

class LinearViscoelasticGenericAxialTorsionCoupling(ConstitutiveLaw):
    """
    Linear viscoelastic generic axial torsion coupling constitutive law
    """

    stiffness: List[List[Union[float, MBVar]]]
    viscosity: Optional[List[List[Union[float, MBVar]]]] = None
    factor: Optional[Union[float, MBVar]] = None
    coupling_coef: float

    def const_law_name(self) -> str:
        return 'linear viscoelastic generic axial torsion coupling'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        
        # String representation for stiffness
        if isinstance(self.stiffness, list):
            N = len(self.stiffness)
            if N != 6:
                raise ValueError("Stiffness matrix must be 6x1")
            matrix_str = ', '.join(str(self.stiffness[i][0]) for i in range(N))
            base_str += f',\n\t{matrix_str}'
        else:
            raise TypeError("Invalid type for stiffness matrix")

        # String representation for viscosity or factor
        if self.viscosity is not None:
            if isinstance(self.viscosity, list):
                N = len(self.viscosity)
                if N != 6:
                    raise ValueError("Viscosity matrix must be 6x1")
                matrix_str = ', '.join(str(self.viscosity[i][0]) for i in range(N))
                base_str += f',\n\t{matrix_str}'
            else:
                raise TypeError("Invalid type for viscosity matrix")
        elif self.factor is not None:
            base_str += f', proportional, {self.factor}'
        else:
            raise ValueError("Either viscosity or factor must be provided")

        # Adding the coupling coefficient
        base_str += f',\n\t{self.coupling_coef}'
        base_str += self.const_law_footer()
        return base_str

class CubicViscoelasticGeneric(ConstitutiveLaw):
    """
    Cubic viscoelastic generic constitutive law
    """

    stiffness_1: Union[float, MBVar, List[Union[float, MBVar]]]
    stiffness_2: Union[float, MBVar, List[Union[float, MBVar]]]
    stiffness_3: Union[float, MBVar, List[Union[float, MBVar]]]
    viscosity: Union[float, MBVar, List[Union[float, MBVar]]]

    def const_law_name(self) -> str:
        return 'cubic viscoelastic generic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        if isinstance(self.stiffness_1, (float, MBVar)):
            base_str += f', {self.stiffness_1}, {self.stiffness_2}, {self.stiffness_3}, {self.viscosity}'
        elif isinstance(self.stiffness_1, list):
            N = len(self.stiffness_1)
            if N == 3:
                stiffness_1_str = ', '.join(str(self.stiffness_1[i]) for i in range(N))
                stiffness_2_str = ', '.join(str(self.stiffness_2[i]) for i in range(N))
                stiffness_3_str = ', '.join(str(self.stiffness_3[i]) for i in range(N))
                viscosity_str = ', '.join(str(self.viscosity[i]) for i in range(N))
                base_str += f',\n\t{stiffness_1_str},\n\t{stiffness_2_str},\n\t{stiffness_3_str},\n\t{viscosity_str}'
            else:
                raise ValueError("Unsupported size of stiffness and viscosity vectors")
        else:
            raise TypeError("Invalid type for stiffness and viscosity values")
        base_str += self.const_law_footer()
        return base_str
    
class DoubleLinearViscoelastic(ConstitutiveLaw):
    """
    Double linear viscoelastic constitutive law
    """

    stiffness_1: Union[MBVar, float]
    upper_strain: Union[MBVar, float]
    lower_strain: Union[MBVar, float]
    stiffness_2: Union[MBVar, float]
    viscosity: Union[MBVar, float]
    viscosity_2: Optional[Union[MBVar, float]] = None

    def const_law_name(self) -> str:
        return 'double linear viscoelastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}, {self.stiffness_1}, {self.upper_strain}, {self.lower_strain}, {self.stiffness_2}, {self.viscosity}'
        if self.viscosity_2 is not None:
            base_str += f', second damping, {self.viscosity_2}'
        base_str += self.const_law_footer()
        return base_str
        
class TurbulentViscoelastic(ConstitutiveLaw):
    """
    Turbulent viscoelastic constitutive law
    """

    stiffness: Union[MBVar, float]
    parabolic_viscosity: Union[MBVar, float]
    threshold: Optional[Union[MBVar, float]] = None
    linear_viscosity: Optional[Union[MBVar, float]] = None

    def const_law_name(self) -> str:
        return 'turbulent viscoelastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}, {self.stiffness}, {self.parabolic_viscosity}'
        if self.threshold is not None:
            base_str += f', {self.threshold}'
            if self.linear_viscosity is not None:
                base_str += f', {self.linear_viscosity}'
        base_str += self.const_law_footer()
        return base_str

class LinearViscoelasticBistop(ConstitutiveLaw):
    """
    Linear viscoelastic bistop constitutive law
    """

    stiffness: Union[float, MBVar]
    viscosity: Union[float, MBVar]
    initial_status: Optional[Union[bool, str]] = None
    activating_condition: DriveCaller
    deactivating_condition: DriveCaller

    def const_law_name(self) -> str:
        return 'linear viscoelastic bistop'

    def __str__(self):
        base_str = f'{self.const_law_header()}, {self.stiffness}, {self.viscosity}'
        if self.initial_status is not None:
            base_str += f',\n\tinitial status, {self.initial_status},'
        base_str += '\n\t# activation condition drive'
        if self.activating_condition.idx is None:
            base_str += f'\n\t{self.activating_condition},'
        else:
            base_str += f'\n\treference, {self.activating_condition.idx},'
        base_str += '\n\t# deactivation condition drive'
        if self.deactivating_condition.idx is None:
            base_str += f'\n\t{self.deactivating_condition}'
        else:
            base_str += f'\n\treference, {self.deactivating_condition.idx}'
        base_str += self.const_law_footer()
        return base_str

class SymbolicElastic(ConstitutiveLaw):
    """
    Symbolic elastic constitutive law
    """

    epsilon: Union[str, List[str]]
    expression: Union[str, List[str]]

    def const_law_name(self) -> str:
        return 'symbolic elastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'

        if isinstance(self.epsilon, str):
            epsilon_list = [self.epsilon]
        else:
            epsilon_list = self.epsilons
        if isinstance(self.expressions, str):
            expression_list = [self.expression]
        else:
            expression_list = self.expression
        if len(epsilon_list) != len(expression_list):
            raise ValueError("The number of epsilons must match the number of expressions")
        
        epsilon_str = ', '.join(f'"{epsilon}"' for epsilon in epsilon_list)
        base_str += f',\n\tepsilon, {epsilon_str}'
        expression_str = ', '.join(f'"{expression}"' for expression in expression_list)
        base_str += f',\n\texpression, {expression_str}'
        base_str += self.const_law_footer()
        return base_str

class SymbolicViscous(ConstitutiveLaw):
    """
    Symbolic viscous constitutive law
    """

    epsilon_prime: Union[str, List[str]]
    expression: Union[str, List[str]]

    def const_law_name(self) -> str:
        return 'symbolic viscous'

    def __str__(self):
        base_str = f'{self.const_law_header()}'

        if isinstance(self.epsilon_prime, str):
            epsilon_prime_list = [self.epsilon_prime]
        else:
            epsilon_prime_list = self.epsilon_prime
        if isinstance(self.expression, str):
            expression_list = [self.expression]
        else:
            expression_list = self.expression
        if len(epsilon_prime_list) != len(expression_list):
            raise ValueError("The number of epsilon_primes must match the number of expressions")

        epsilon_prime_str = ', '.join(f'"{epsilon_prime}"' for epsilon_prime in epsilon_prime_list)
        base_str += f',\n\tepsilon prime, {epsilon_prime_str}'
        expression_str = ', '.join(f'"{expression}"' for expression in expression_list)
        base_str += f',\n\texpression, {expression_str}'
        base_str += self.const_law_footer()
        return base_str
    
class SymbolicViscoelastic(ConstitutiveLaw):
    """
    Symbolic viscoelastic constitutive law
    """

    epsilon: Union[str, List[str]]
    epsilon_prime: Union[str, List[str]]
    expression: Union[str, List[str]]

    def const_law_name(self) -> str:
        return 'symbolic viscoelastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'

        if isinstance(self.epsilon, str):
            epsilon_list = [self.epsilon]
        else:
            epsilon_list = self.epsilon
        if isinstance(self.epsilon_prime, str):
            epsilon_prime_list = [self.epsilon_prime]
        else:
            epsilon_prime_list = self.epsilon_prime
        if isinstance(self.expression, str):
            expression_list = [self.expression]
        else:
            expression_list = self.expression
        if len(epsilon_list) != len(epsilon_prime_list) or len(epsilon_list) != len(expression_list):
            raise ValueError("The number of epsilons, epsilon_primes, and expressions must match")

        epsilon_str = ', '.join(f'"{epsilon}"' for epsilon in epsilon_list)
        epsilon_prime_str = ', '.join(f'"{epsilon_prime}"' for epsilon_prime in epsilon_prime_list)
        expression_str = ', '.join(f'"{expression}"' for expression in expression_list)
        base_str += f',\n\tepsilon, {epsilon_str}'
        base_str += f',\n\tepsilon prime, {epsilon_prime_str}'
        base_str += f',\n\texpression, {expression_str}'
        base_str += self.const_law_footer()
        return base_str
    
class SymbolicViscoelastic(ConstitutiveLaw):
    """
    Symbolic viscoelastic constitutive law
    """

    epsilon: Union[str, List[str]]
    epsilon_prime: Union[str, List[str]]
    expression: Union[str, List[str]]

    def const_law_name(self) -> str:
        return 'symbolic viscoelastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'

        if isinstance(self.epsilon, str):
            epsilon_list = [self.epsilon]
        else:
            epsilon_list = self.epsilon
        if isinstance(self.epsilon_prime, str):
            epsilon_prime_list = [self.epsilon_prime]
        else:
            epsilon_prime_list = self.epsilon_prime
        if isinstance(self.expression, str):
            expression_list = [self.expression]
        else:
            expression_list = self.expression
        if len(epsilon_list) != len(epsilon_prime_list) or len(epsilon_list) != len(expression_list):
            raise ValueError("The number of epsilons, epsilon_primes, and expressions must match")

        epsilon_str = ', '.join(f'"{epsilon}"' for epsilon in epsilon_list)
        epsilon_prime_str = ', '.join(f'"{epsilon_prime}"' for epsilon_prime in epsilon_prime_list)
        expression_str = ', '.join(f'"{expression}"' for expression in expression_list)
        base_str += f',\n\tepsilon, {epsilon_str}'
        base_str += f',\n\tepsilon prime, {epsilon_prime_str}'
        base_str += f',\n\texpression, {expression_str}'
        base_str += self.const_law_footer()
        return base_str
    
class AnnElastic(ConstitutiveLaw):
    """
    Ann elastic constitutive law
    """

    file_name: str

    def const_law_name(self) -> str:
        return 'ann elastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        base_str += f',\n\t"{self.file_name}"'
        base_str += self.const_law_footer()
        return base_str

class AnnViscoelastic(ConstitutiveLaw):
    """
    Ann viscoelastic constitutive law
    """

    file_name: str

    def const_law_name(self) -> str:
        return 'ann viscoelastic'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        base_str += f',\n\t"{self.file_name}"'
        base_str += self.const_law_footer()
        return base_str

class ArrayConstitutiveLaw(ConstitutiveLaw):
    """
    Array constitutive law wrapper linearly combines the output of multiple constitutive laws.
    """

    number: int
    wrapped_const_laws: List[ConstitutiveLaw]

    def const_law_name(self) -> str:
        return 'array'

    def __str__(self):
        if self.number == 1:
            return str(self.wrapped_const_laws[0])

        base_str = f'{self.const_law_header()}, {self.number}'
        for law in self.wrapped_const_laws:
            base_str += f',\n\t{str(law)}'
        base_str += self.const_law_footer()
        return base_str

class BistopConstitutiveLaw(ConstitutiveLaw):
    """
    Bistop wrapper applies the logic of the bistop to a generic underlying constitutive law.
    """

    initial_status: Optional[Union[bool, str]]
    activating_condition: DriveCaller
    deactivating_condition: DriveCaller
    wrapped_const_law: ConstitutiveLaw

    def const_law_name(self) -> str:
        return 'bistop'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        if self.initial_status is not None:
            base_str += f',\n\tinitial status, {self.initial_status},'

        base_str += '\n\t# activation condition drive'
        if self.activating_condition.idx is None:
            base_str += f'\n\t{self.activating_condition},'
        else:
            base_str += f'\n\treference, {self.activating_condition.idx},'
        base_str += '\n\t# deactivation condition drive'
        if self.deactivating_condition.idx is None:
            base_str += f'\n\t{self.deactivating_condition},'
        else:
            base_str += f'\n\treference, {self.deactivating_condition.idx},'

        base_str += f'\n\t{str(self.wrapped_const_law)}'
        base_str += self.const_law_footer()
        return base_str

class InvariantAngularWrapper(ConstitutiveLaw):
    """
    Invariant angular wrapper for 3D constitutive laws used within the "attached" variant of the deformable hinge joint.
    """

    xi: Union[float, int, MBVar]
    wrapped_const_law: ConstitutiveLaw

    def const_law_name(self) -> str:
        return 'invariant angular'

    def __str__(self):
        base_str = f'{self.const_law_header()}'
        base_str += f',\n\t{self.xi}'
        base_str += f',\n\t{str(self.wrapped_const_law)}'
        base_str += self.const_law_footer()
        return base_str
    
class NamedConstitutiveLaw(MBEntity):
    """
    Adapter for using a constitutive law that is not yet implemented in the preprocessor
    as a regular `ConstitutiveLaw` subclass with argument checking.

    This should only be used temporarily, and may be removed in the future without prior notice.
    """
    
    content: str
    """Text that will be output to MBDyn for this law"""

    def __init__(self, law: Union[str, list]):
        warnings.warn(
            "Using a string for constitutive laws is not recommended " + \
            "and may be removed in the future. " + \
            "Consider using ConstitutiveLaw instances for better support.",
            UserWarning
        )
        if isinstance(law, list):
            law = ', '.join(str(l) for l in law)
        else:
            law = str(law)
        super().__init__(content=law)

    def __str__(self):
        return self.content

if imported_pydantic:
    DeformableAxial.model_rebuild()
    DeformableHinge.model_rebuild()
    Rod.model_rebuild()
    RodWithOffset.model_rebuild()
    RodBezier.model_rebuild()
    ViscousBody.model_rebuild()
    DeformableDisplacement.model_rebuild()
    DeformableJoint.model_rebuild()
    Beam.model_rebuild()
    AerodynamicBeam.model_rebuild()
    AerodynamicBody.model_rebuild()


class FileDriver(MBEntity):
    """
    Abstract class for file drivers. The file drivers are defined by the statement:
    file : <file_arglist> ;
    A comprehensive family of file drivers is available.
    """
    
    idx: Union[MBVar, int]
    """Index of this file driver"""

    @abstractmethod
    def driver_type(self) -> str:
        """Every file driver must have a type"""
        pass

    def file_header(self) -> str:
        """Common syntax for start of any file driver"""
        return f'file: {self.idx}, {self.driver_type()}'
    
class FixedStep(FileDriver):
    """
    Fixed Step file driver
    """
    
    class InterpolationType(Enum):
        LINEAR = "linear"
        CONST = "const"

    class BailoutType(Enum):
        NONE = "none"
        UPPER = "upper"
        LOWER = "lower"
        ANY = "any"

    class PadZeroesType(Enum):
        YES = 'yes'
        NO = 'no'
    
    steps_number: Union[int, MBVar, str]  # 'count' or specific number of steps
    columns_number: Union[int, MBVar]
    initial_time: Union[float, MBVar, str]  # 'from file' or specific initial time
    time_step: Union[float, MBVar, str]  # 'from file' or specific time step
    interpolation: Optional[InterpolationType]
    pad_zeroes: Optional[PadZeroesType]
    bailout: Optional[BailoutType]
    file_name: str

    def driver_type(self) -> str:
        return "fixed step"
    
    def __str__(self):
        base_str = f'{self.file_header()},\n\t'
        base_str += f'{self.steps_number},\n\t'
        base_str += f'{self.columns_number},\n\t'
        base_str += f'initial time, {self.initial_time},\n\t'
        base_str += f'time step, {self.time_step},\n\t'
        if self.interpolation:
            base_str += f'interpolation, {self.interpolation.value},\n\t'
        if self.pad_zeroes == FixedStep.PadZeroesType.NO and self.bailout:
            raise ValueError("Cannot set both 'pad zeroes' to 'no' and 'bailout' at the same time")
        if self.pad_zeroes:
            base_str += f'pad zeroes, {self.pad_zeroes},\n\t'
        elif self.bailout:
            base_str += f'bailout, {self.bailout},\n\t'
        base_str += f'"{self.file_name}"'
        return base_str
    
class VariableStep(FileDriver):
    """
    Variable Step file driver
    """
    
    channels_number: Union[int, MBVar]
    interpolation: Optional[FixedStep.InterpolationType]
    pad_zeroes: Optional[FixedStep.PadZeroesType]
    bailout: Optional[FixedStep.BailoutType]
    file_name: str

    def driver_type(self) -> str:
        return "variable step"
    
    def __str__(self):
        base_str = f'{self.file_header()},\n\t'
        base_str += f'{self.channels_number},\n\t'
        if self.interpolation:
            base_str += f'interpolation, {self.interpolation.value},\n\t'
        if self.pad_zeroes == FixedStep.PadZeroesType.NO and self.bailout:
            raise ValueError("Cannot set both 'pad zeroes' to 'no' and 'bailout' at the same time")
        if self.pad_zeroes:
            base_str += f'pad zeroes, {self.pad_zeroes},\n\t'
        elif self.bailout:
            base_str += f'bailout, {self.bailout},\n\t'
        base_str += f'"{self.file_name}"'
        return base_str


class Data(MBEntity):
    problem: Union[Literal["initial value"], Literal["inverse dynamics"]] = "initial value"

    def __str__(self):
        s = 'begin: data;\n'
        s += f'\tproblem: {self.problem};\n'
        s += 'end: data;'
        return s
    

# Initial Value 
class Strategy(MBEntity):
    """
    Abstract base class for all Strategies
    """

    @abstractmethod
    def strategy_type(self) -> str:
        """Every strategy class must define this to return its MBDyn syntax name"""
        raise NotImplementedError("called strategy_type of abstract Element")

    def strategy_header(self) -> str:
        """common syntax for start of any strategy"""
        return f'stratefy: {self.strategy_type()}'

class StrategyFactor(Strategy):
    reduction_factor: Union[float, MBVar]
    steps_before_reduction: Union[int, MBVar]
    raise_factor: Union[float, MBVar]
    steps_before_raise: Union[int, MBVar]
    min_iterations: Union[int, MBVar]
    max_iterations: Optional[Union[int, MBVar]] = None

    def __str__(self):
        s = f'{self.strategy_header()}, {self.steps_before_reduction}'
        s += f', {self.raise_factor}, {self.steps_before_raise}, {self.min_iterations}'
        if self.max_iterations is not None:
            s += f', {self.max_iterations}'
        return s

class StrategyChange(Strategy):
    time_step_pattern: DriveCaller

    def __str__(self):
        s = f'{self.strategy_header()}, {self.time_step_pattern}'
        return s

class StrategyNoChange(Strategy):
    def __str__(self):
        s = f'{self.strategy_header()}'
        return s
    
class Tolerance(MBEntity):
    residual_tolerance: Union[Literal['null'], float, MBVar]
    residual_test: Optional[Literal['none', 'norm', 'minmax']] = None
    scaling: Optional[str] = None
    solution_tolerance: Optional[Union[Literal['null'], float, MBVar]] = None
    solution_test: Optional[Literal['none', 'norm', 'minmax']] = None

    @field_validator('scaling')
    def validate_scaling(cls, v, info):
        residual_test = info.data.get('residual_test')
        if residual_test is None and v is not None:
            raise ValueError("scaling should be None if residual_test is not provided")
        if v is not None and v.lower() != 'scale':
            raise ValueError("scaling must be a string literal: 'scale'")
        return v.lower() if v is not None else v

    @field_validator('solution_test', mode='before')
    def validate_solution_test(cls, v, info):
        solution_tolerance = info.data.get('solution_tolerance')
        if solution_tolerance is None and v is not None:
            raise ValueError("solution_test should be None if solution_tolerance is not provided")
        return v

    def __str__(self):
        s = f'tolerance: {self.residual_tolerance}'
        if self.residual_test is not None:
            s += f', test, {self.residual_test}'
            if self.scaling is not None:
                s += f', {self.scaling}'
        if self.solution_tolerance is not None:
            s += f', {self.solution_tolerance}'
            if self.solution_test is not None:
                s += f', test, {self.solution_test}'
        return s

class MaxIterations(MBEntity):
    '''Error out after max_iterations without passing the convergence test. The default value is zero.'''

    max_iterations: int = 0
    optional_keywords: Optional[Literal['at most']] = None

    def __str__(self):
        s = f'max iterations: {self.max_iterations}'
        if self.optional_keywords is not None:
            s += f', {self.optional_keywords}'
        return s

class Method(MBEntity):
    """Base class for every method"""

    @abstractmethod
    def __str__(self) -> str:
        """Has to be overridden to output the MBDyn syntax"""
        pass

class CrankNicolson(Method):
    def __str__(self):
        return 'method: crank nicolson'

class MethodWithRadius(Method):
    differential_radius: DriveCaller
    algebraic_radius: Optional[DriveCaller] = None
    
    def __str__(self):
        s = f'method: {self.__class__.__name__.lower()}, {self.differential_radius}'
        if self.algebraic_radius is not None:
            s += f', {self.algebraic_radius}'
        return s

class MS(MethodWithRadius):
    """The 'ms' method (also referred to as 'ms2') allows for tuning algorithmic dissipation."""
    pass  # Inherits the behavior from MethodWithRadius

class MS2(MethodWithRadius):
    pass  # Inherits the behavior from MethodWithRadius

class MS3(MethodWithRadius):
    """The 'ms3' method is a three-step method allowing for algorithmic dissipation tuning."""
    pass  # Inherits the behavior from MethodWithRadius

class MS4(MethodWithRadius):
    """The 'ms4' method is a four-step method allowing for algorithmic dissipation tuning."""
    pass  # Inherits the behavior from MethodWithRadius

class Hope(MethodWithRadius):
    """The 'hope' method is a multi-stage method combining Crank-Nicolson and 'ms' methods."""
    pass  # Inherits the behavior from MethodWithRadius

class SS2(MethodWithRadius):
    pass

class SS3(MethodWithRadius):
    pass

class SS4(MethodWithRadius):
    pass

class Bathe(MethodWithRadius):
    pass

class MSSTC3(MethodWithRadius):
    pass

class MSSTC4(MethodWithRadius):
    pass

class MSSTC5(MethodWithRadius):
    pass

class MSSTH3(MethodWithRadius):
    pass

class MSSTH4(MethodWithRadius):
    pass

class MSSTH5(MethodWithRadius):
    pass

class Hybrid(MethodWithRadius):
    default_hybrid_method: Literal['implicit euler', 'crank nicolson', 'ms2', 'hope']
    
    def __str__(self):
        s = f'method: hybrid, {self.default_hybrid_method}, {self.differential_radius}'
        if self.algebraic_radius is not None:
            s += f', {self.algebraic_radius}'
        return s

class DIRK33(Method):
    def __str__(self):
        return 'method: DIRK33'

class DIRK43(Method):
    def __str__(self):
        return 'method: DIRK43'

class DIRK54(Method):
    def __str__(self):
        return 'method: DIRK54'

class BDF(Method):
    order: Optional[Union[int, MBVar]] = None  # 1 or 2

    def __str__(self):
        s = 'method: bdf'
        if self.order is not None:
            s += f', order, {self.order}'
        return s

class ImplicitEuler(Method):
    def __str__(self):
        return 'method: implicit euler'
    
class NonlinearSolver(MBEntity):
    """The nonlinear solver solves a nonlinear problem F (x) = 0."""

    @abstractmethod
    def nonlinear_solver_name(self) -> str:
        """Name of the specific nonlinear solver"""
        pass

    def nonlinear_solver_header(self) -> str:
        """Common syntax for start of any nonlinear solver"""
        return f'nonlinear solver: {self.nonlinear_solver_name()}'
    
class NewtonRaphson(NonlinearSolver):
    pass
    

class MethodforEigenanalysis(MBEntity):
    '''Base class for the methods used in Eigenanalysis'''

    @abstractmethod
    def __str__(self) -> str:
        """Has to be overridden to output the MBDyn syntax"""
        pass

class UseLapack(MethodforEigenanalysis):
    balance: Optional[Literal['no', 'scale', 'permute', 'all']] = None

    def __str__(self):
        s = 'use lapack'
        if self.balance is not None:
            s += f', balance, {self.balance}'
        return s
    
class UseArpack(MethodforEigenanalysis):
    nev: Union[int, MBVar]
    ncv: Union[int, MBVar]
    tol: Union[float, MBVar]
    max_iter: Optional[Union[int, MBVar]] = 300

    @field_validator('tol')
    def check_tolerance(cls, v):
        if v < 0:
            raise ValueError("Tolerance (tol) must be positive. Use zero for machine precision.")
        return v

    def __str__(self):
        s = f'use arpack, {self.nev}, {self.ncv}, {self.tol}'
        if self.max_iter != 300:
            s += f', max iterations, {self.max_iter}'
        return s

class UseJdqz(MethodforEigenanalysis):
    nev: Union[int, MBVar]
    ncv: Union[int, MBVar]
    tol: Union[float, MBVar]

    @field_validator('tol')
    def check_tolerance(cls, v):
        if v < 0:
            raise ValueError("Tolerance (tol) must be positive. Use zero for machine precision.")
        return v

    def __str__(self):
        return f'use arpack, {self.nev}, {self.ncv}, {self.tol}'

class UseExternal(MethodforEigenanalysis):
    def __str__(self):
        return 'use external'

class Eigenanalysis(MBEntity):
    '''
    Performs the direct eigenanalysis of the problem. This functionality is experimental. Direct 
    eigenanalysis based on the matrices of the system only makes sense when the system is in a 
    steady conﬁguration, so the user needs to ensure this conﬁguration has been reached.
    Moreover, not all elements currently contribute to the Jacobian matrix of the system, so YMMV. In case
    of rotating systems, a steady conﬁguration could be reached when the model is expressed in a relative
    reference frame, using the rigid body kinematics card.
    '''

    # Mode Options Enum for eigenvalue sorting criteria
    class ModeOptions(Enum):
        SMALLEST_MAGNITUDE = "smallest magnitude"
        LARGEST_MAGNITUDE = "largest magnitude"
        LARGEST_REAL_PART = "largest real part"
        SMALLEST_REAL_PART = "smallest real part"
        LARGEST_IMAGINARY_PART = "largest imaginary part"
        SMALLEST_IMAGINARY_PART = "smallest imaginary part"

    num_times: Optional[Union[int, MBVar]] = None
    when: Union[float, MBVar, List[Union[float, MBVar]]]
    suffix_width: Optional[Union[float, MBVar, Literal['compute']]] = None
    suffix_format: Optional[str] = None
    output_full_matrices: Optional[bool] = None
    output_sparse_matrices: Optional[bool] = None
    output_eigenvectors: Optional[bool] = None
    output_geometry: Optional[bool] = None
    matrix_precision: Optional[Union[float, MBVar]] = None
    results_precision: Optional[Union[float, MBVar]] = None
    parameter: Optional[Union[float, MBVar]] = None
    mode_options: Optional[ModeOptions] = None
    lower_frequency_limit: Optional[Union[float, MBVar]] = None
    upper_frequency_limit: Optional[Union[float, MBVar]] = None
    method: Optional[MethodforEigenanalysis] = None

    @model_validator(mode='after')
    def check_when_and_num_times(cls, self):
        when = self.when
        num_times = self.num_times
        if isinstance(when, list) and num_times is None:
            raise ValueError("If 'when' is given as a list, 'num_times' must also be provided.")
        return self

    def add_optional_field(self, s, field_name, field_value):
        if field_value is True:
            return s + f',\n\t{field_name}'
        elif field_value is not None and field_value is not False:
            return s + f',\n\t{field_name}, {field_value}'
        return s

    def __str__(self):
        s = 'eigenanalysis: '
        if isinstance(self.when, List):
            s += f'\n\tlist, {self.num_times}, '
            s += ', '.join(str(i) for i in self.when)
        else:
            s += f'\n\t{self.when}'
        
        s = self.add_optional_field(s, 'suffix width', self.suffix_width)
        s = self.add_optional_field(s, 'suffix format', self.suffix_format)
        s = self.add_optional_field(s, 'output full matrices', self.output_full_matrices)
        s = self.add_optional_field(s, 'output sparse matrices', self.output_sparse_matrices)
        s = self.add_optional_field(s, 'output eigenvectors', self.output_eigenvectors)
        s = self.add_optional_field(s, 'output geometry', self.output_geometry)
        s = self.add_optional_field(s, 'matrix output precision', self.matrix_precision)
        s = self.add_optional_field(s, 'results output precision', self.results_precision)
        s = self.add_optional_field(s, 'parameter', self.parameter)
        s = self.add_optional_field(s, 'mode', self.mode_options)
        s = self.add_optional_field(s, 'lower frequency limit', self.lower_frequency_limit)
        s = self.add_optional_field(s, 'upper frequency limit', self.upper_frequency_limit)
        if self.method is not None:
            s += f',\n\t{self.method}'
        
        return s

        # TODO: Delete these lines of code after testing of the new approach
        # if self.suffix_width is not None:
        #     s += f',\n\tsuffix width, {self.suffix_width}'
        # if self.suffix_format is not None:
        #     s += f',\n\tsuffix format, {self.suffix_format}'
        # if self.output_full_matrices is True:
        #     s += f',\n\toutput full matrices'
        # if self.output_sparse_matrices is True:
        #     s += f',\n\toutput sparse matrices'
        # if self.output_eigenvectors is True:
        #     s += f',\n\toutput eigenvectors'
        # if self.output_geometry is True:
        #     s += f',\n\toutput geometry'
        # if self.matrix_precision is not None:
        #     s += f',\n\tmatrix output precision, {self.matrix_precision}'
        # if self.results_precision is not None:
        #     s += f',\n\tresults output precision, {self.results_precision}'
        # if self.parameter is not None:
        #     s += f',\n\tparameter, {self.parameter}'
        # if self.mode_options is not None:
        #     s += f',\n\tmode, {self.mode_options}'
        # if self.lower_frequency_limit is not None:
        #     s += f',\n\tlower frequency limit, {self.lower_frequency_limit}'
        # if self.upper_frequency_limit is not None:
        #     s += f',\n\tupper frequency limit, {self.upper_frequency_limit}'
        # if self.method is not None:
        #     s += f',\n\t{self.method}'
        # return s


class LinearSolver(MBEntity):
    solver_name: Literal[
        'naive', 'umfpack', 'klu', 'y12', 'lapack', 'superlu', 'taucs', 
        'pardiso', 'pardiso_64', 'watson', 'pastix', 'qr', 'spqr', 
        'aztecoo', 'amesos', 'siconos dense', 'siconos sparse'
    ]

    # General solver settings
    storage_mode: Optional[Literal['map', 'cc', 'dir', 'grad']] = None
    ordering: Optional[Literal['colamd', 'mmdata', 'amd', 'given', 'metis']] = None
    
    # Threading configuration
    multithread: Optional[Literal['mt', 'multithread']] = None
    threads: Optional[Union[int, MBVar]] = None
    
    # Solver-specific parameters
    workspace_size: Optional[Union[int, MBVar]] = None
    pivot_factor: Optional[Union[float, MBVar]] = None
    drop_tolerance: Optional[Union[float, MBVar]] = None
    block_size: Optional[Union[int, MBVar]] = None
    
    # Scaling options
    scale: Optional[Literal[
        'no', 'always', 'once', 'row max', 'row sum', 'column max', 
        'column sum', 'lapack', 'iterative', 'row max column max'
    ]] = None
    scale_tolerance: Optional[Union[float, MBVar]] = None
    scale_max_iter: Optional[Union[int, MBVar]] = None
    
    # Refinement and tolerance settings
    refine_tolerance: Optional[Union[float, MBVar]] = None
    refine_max_iter: Optional[Union[int, MBVar]] = None

    # Preconditioner options
    preconditioner: Optional[Literal[
        'umfpack', 'klu', 'lapack', 'ilut', 'superlu', 'mumps',
        'scalapack', 'dscpack', 'pardiso', 'paraklete', 'taucs', 'csparse'
    ]] = None

    @model_validator(mode="before")
    def check_solver_specific_parameters(cls, values):
        # Check if multithread is set but threads is None
        if values.get('multithread') is not None and values.get('threads') is None:
            raise ValueError("If multithread is set, threads must also be specified.")
        
        solver = values.get('solver_name')
        
        # Check parameters specific to 'umfpack'
        if solver == 'umfpack':
            values.setdefault('block_size', 32)
            if values.get('workspace_size') is not None:
                raise ValueError("workspace_size is ignored for umfpack solver.")
            if values.get('drop_tolerance') is None:
                values['drop_tolerance'] = 0.0  # Default drop tolerance for umfpack

        # Enforce ordering for naive solver
        if solver == 'naive' and values.get('ordering') is None:
            raise ValueError("The naive solver requires an ordering option for robustness, e.g., 'colamd'.")

        # Enforce refine_max_iter for certain solvers
        if solver in ['pardiso', 'pardiso_64', 'pastix'] and values.get('refine_max_iter') is None:
            raise ValueError(f"{solver} requires refine_max_iter for stability.")

        # Enforce ignore of certain parameters for specific solvers
        if solver in ['naive', 'y12', 'lapack'] and values.get('workspace_size') is not None:
            raise ValueError(f"workspace_size is ignored for {solver} solver.")
        
        if solver == 'klu' and values.get('scale') not in [None, 'always', 'once']:
            raise ValueError("KLU solver supports only 'always' or 'once' scale options.")

        # Ensure that iterative refinement settings are consistent
        if solver in ['umfpack', 'pastix'] and values.get('refine_max_iter') and values.get('refine_tolerance') is None:
            raise ValueError("Refinement tolerance is required if refine_max_iter is set.")

        # TODO: Have to check thoroughly 
        # # Check for supported keywords by solvers
        # supported_keywords = {
        #     'umfpack': ['map', 'cc', 'dir', 'drop_tolerance', 'block_size', 'scale', 'refine_max_iter'],
        #     'klu': ['map', 'cc', 'dir', 'scale', 'refine_max_iter'],
        #     'y12': ['map', 'dir'],
        #     'superlu': ['map', 'cc', 'scale'],
        #     'pastix': ['map', 'cc', 'scale', 'refine_max_iter'],
        #     'naive': ['cc', 'scale', 'colamd', 'mmdata'],
        #     'spqr': ['colamd', 'amd', 'metis', 'given'],
        # }

        # if solver in supported_keywords and 'keywords' in values:
        #     invalid_keywords = [kw for kw in values['keywords'] if kw not in supported_keywords[solver]]
        #     if invalid_keywords:
        #         raise ValueError(f"The following keywords are not supported by the {solver} solver: {', '.join(invalid_keywords)}.")

        # Check for pivot factor validity
        if values.get('pivot_factor') is not None:
            if not (0.0 <= values['pivot_factor'] <= 1.0):
                raise ValueError("pivot_factor must be between 0.0 and 1.0.")

        # Check for drop tolerance with unsupported solvers
        if solver != 'umfpack' and values.get('drop_tolerance') is not None:
            raise ValueError("drop_tolerance can only be used with the umfpack solver.")

        # Check for inconsistent scaling options
        if solver not in ['naive', 'klu', 'umfpack', 'pastix'] and values.get('scale') is not None:
            raise ValueError(f"scale option is not supported by the {solver} solver.")

        # Check for valid block size for umfpack
        if solver != 'umfpack' and values.get('block_size') is not None:
            raise ValueError("block_size can only be used with the umfpack solver.")

        return values

    def __str__(self):
        s = f'linear solver: {self.solver_name}'
        if self.storage_mode is not None:
            s += f', {self.storage_mode}'
        if self.ordering is not None:
            s += f', {self.ordering}'
        if self.multithread is not None:
            s += f',\n\t\t{self.multithread}, {self.threads}'
        if self.workspace_size is not None:
            s += f',\n\t\tworkspace size, {self.workspace_size}'
        if self.pivot_factor is not None:
            s += f',\n\t\tpivot factor, {self.pivot_factor}'
        if self.drop_tolerance is not None:
            s += f',\n\t\tdrop tolerance, {self.drop_tolerance}'
        if self.block_size is not None:
            s += f',\n\t\tblock size, {self.block_size}'
        if self.scale is not None:
            s += f',\n\t\tscale, {self.scale}'
            if self.scale_tolerance is not None:
                s += f',\n\t\t\tscale tolerance, {self.scale_tolerance}'
            if self.scale_max_iter is not None:
                s += f',\n\t\t\tscale iterations, {self.scale_max_iter}'
        if self.refine_tolerance is not None:
            s += f',\n\t\ttolerance, {self.refine_tolerance}'
        if self.refine_max_iter is not None:
            s += f',\n\t\tmax iterations, {self.refine_max_iter}'
        if self.preconditioner is not None:
            s += f',\n\t\tpreconditioner, {self.preconditioner}'
        return s
    
class Threads(MBEntity):
    mode: Literal['auto', 'disable', 'assembly', 'solver'] = 'auto'
    threads: Optional[Union[int, MBVar]] = None

    @model_validator(mode='after')
    def check_threads_provided(cls, model):
        mode = model.mode
        threads = model.threads

        if mode in ['assembly', 'solver'] and threads is None:
            raise ValueError("threads must be provided if mode is 'assembly' or 'solver'.")
        elif mode in ['auto', 'disable'] and threads is not None: 
            raise ValueError("threads must be None if mode is 'auto' or 'disable'.")
        return model
    
    def __str__(self):
        s = f'threads: {self.mode}'
        if self.threads:
            s += f', {self.threads}'
        return s

class DerivativesCoefficient(MBEntity):
    coefficient: Optional[Union[float, MBVar]] = None
    is_auto: Optional[bool] = False
    max_iterations: Optional[Union[int, MBVar]] = None
    factor: Optional[Union[float, MBVar]] = None

    @model_validator(mode='before')
    @classmethod
    def validate_auto_case(cls, data: dict):
        is_auto = data.get('is_auto', False)  # Default to False if not provided
        coefficient = data.get('coefficient')
        factor = data.get('factor')
        max_iterations = data.get('max_iterations')

        if not is_auto:
            if coefficient is None:
                raise ValueError("When 'auto' is not selected, a numeric value for 'coefficient' must be specified.")
            if factor is not None or max_iterations is not None:
                raise ValueError("`factor` and `max_iterations` can only be specified when 'auto' is selected.")
        
        return data

    def __str__(self):
        s = 'derivatives coefficient: '
        if self.is_auto:
            s += f"{f'{self.coefficient}, ' if self.coefficient is not None else ''}auto"
        else:
            s += f'{self.coefficient}'
        if self.max_iterations is not None:
            s += f",\n\tmax iterations, {self.max_iterations}"
        if self.factor is not None:
            s += f",\n\tfactor, {self.factor}"
        return s


class OutputSettings(MBEntity):
    items: List[Literal[
        "iterations", "residual", "solution", "jacobian matrix", 
        "messages", "counter", "bailout", "matrix condition number", 
        "solver condition number", "cpu time", "none"
    ]]

    # TODO: Have to get a review
    @model_validator(mode="before")
    @classmethod
    def validate_items(cls, data: dict):
        items = data.get("items", [])
        # Ensure the 'none' keyword is used alone or as the first item
        if "none" in items and items[0] != "none":
            raise ValueError("If 'none' is specified, it must be the first item.")
        return data

    def __str__(self):
        return f"output: {', '.join(self.items)}"

class InitialValue(MBEntity):
    '''
    At present, the main problem is initial value, which solves initial value problems by means of generic
    integration schemes that can be cast in a broad family of multistep and, experimentally, Implicit Runge-
    Kutta-like schemes
    '''

    initial_time: Union[float, MBVar]
    final_time: Union[float, MBVar, Literal["forever"]]
    strategy: Optional[Union[StrategyChange, StrategyFactor, StrategyNoChange]] = None
    min_time_step: Optional[Union[float, MBVar]] = None
    max_time_step: Optional[Union[float, MBVar, Literal["unlimited"]]] = None
    time_step: Union[float, MBVar]
    tolerance: Tolerance
    max_iterations: MaxIterations
    modify_residual_test: Optional[Union[bool, int]] = False    # 0 / 1 / True / False
    method: Optional[Method] = None
    eigenanalysis: Optional[Eigenanalysis] = None
    linear_solver: Optional[LinearSolver] = None
    threads: Optional[Threads] = None
    derivatives_tolerance: Optional[Union[float, MBVar]] = None
    derivatives_max_iterations: Optional[Union[int, MBVar]] = None
    derivatives_coefficient: Optional[DerivativesCoefficient] = None
    output_settings: Optional[OutputSettings] = None
    output_meter: Optional[DriveCaller] = None

    @field_validator('modify_residual_test')
    def set_modify_residual_test(cls, v):
        if isinstance(v, (int, bool)):
            if v in [0, False]:
                return None
            elif v in [1, True]:
                return "modify residual test"
            else:
                raise ValueError("modify_residual_test must be 0, 1, True, or False.")
        else:
            raise TypeError("modify_residual_test must be of type int or bool.")
        
    def __str__(self):
        s = "begin: initial value;\n"
        s += f"\tinitial time: {self.initial_time};\n"
        s += f"\tfinal time: {self.final_time};\n"
        if self.strategy:
            s += f"\t{self.strategy};\n"
        if self.min_time_step:
            s += f"\tmin time step: {self.min_time_step};\n"
        if self.max_time_step:
            s += f"\tmax time step: {self.max_time_step};\n"
        s += f"\ttime step: {self.time_step};\n"
        s += f"\t{self.max_iterations};\n"
        s += f"\t{self.tolerance};\n"
        if self.modify_residual_test:
            s += f"\tmodify residual test;\n"
        if self.method:
            s += f"\t{self.method};\n"
        if self.eigenanalysis:
            s += f"\t{self.eigenanalysis};\n"
        if self.linear_solver:
            s += f"\t{self.linear_solver};\n"
        if self.threads:
            s += f"\t{self.threads};\n"
        if self.derivatives_tolerance:
            s += f"\tderivatives tolerance: {self.derivatives_tolerance};\n"
        if self.derivatives_max_iterations:
            s += f"\tderivatives max iterations: {self.derivatives_max_iterations};\n"
        if self.derivatives_coefficient:
            s += f"\t{self.derivatives_coefficient};\n"
        if self.output_settings:
            s += f"\t{self.output_settings};\n"
        if self.output_meter:
            s += f"\toutput meter: {self.output_meter};\n"
        s += "end: initial value;\n\n"
        return s

# Control Data
class Print(MBEntity):
    items: List[Literal[
        "dof stats", "dof description", "equation description", "description", 
        "element connection", "node connection", "connection", "all", "none"
    ]] = []
    item_to_file: Optional[List[bool]] = None  

    def __str__(self):
        s = "print: "

        # Check if item_to_file has the same length as items
        item_to_file_filled = self.item_to_file or [False] * len(self.items)
        if len(item_to_file_filled) < len(self.items):
            item_to_file_filled.extend([False] * (len(self.items) - len(item_to_file_filled)))
            
        # Loop through each item and build the string with appropriate formatting
        for idx, item in enumerate(self.items):
            s += f"{item}"
            # Check if 'to file' is specified for this item using item_to_file_filled
            if item_to_file_filled[idx]:  
                s += ", to file"
            # Add a comma between items except for the last item
            if idx < len(self.items) - 1:  
                s += ", "
        return s
    
class OutputResults(MBEntity):
    '''
    This deprecated statement was intended for producing output in formats compatible with other software.
    Most of them are produced in form of post-processing, based on the default raw output.
    '''

    file_format: Literal["classic", "classic64", "nc4", "nc4classic"] = "nc4"
    sync: bool = False
    text: bool = False

    def __str__(self):
        s = f'output results: netcdf, {self.file_format}'
        if self.sync:
            s += ',\n\tsync'
        else: 
            s += ',\n\tno sync'
        if self.text:
            s += ', text'
        else: 
            s += ', no text'
        return s
    
class ConstRBK(MBEntity):
    position: Optional[Position] = None
    orientation: Optional[Position] = None
    velocity: Optional[Position] = None
    angular_velocity: Optional[Position] = None
    acceleration: Optional[Position] = None
    angular_acceleration: Optional[Position] = None

    def __str__(self):
        s = 'const'
        if self.position:
            s += f',\n\tposition, {self.position}'
        if self.orientation:
            s += f',\n\torientation, {self.orientation}'
        if self.velocity:
            s += f',\n\tvelocity, {self.velocity}'
        if self.angular_velocity:
            s += f',\n\tangular velocity, {self.angular_velocity}'
        if self.acceleration:
            s += f',\n\tacceleration, {self.acceleration}'
        if self.angular_acceleration:
            s += f',\n\tangular acceleration, {self.angular_acceleration}'
        return s

class DriveRBK(MBEntity):
    position: Optional[List] = None
    orientation: Optional[List] = None
    velocity: Optional[List] = None
    angular_velocity: Optional[List] = None
    acceleration: Optional[List] = None
    angular_acceleration: Optional[List] = None

    def __str__(self):
        s = 'drive'
        if self.position:
            s += f',\n\tposition, {", ".join(str(i) for i in self.position)}'
        if self.orientation:
            s += f',\n\torientation, {", ".join(str(i) for i in self.orientation)}'
        if self.velocity:
            s += f',\n\tvelocity, {", ".join(str(i) for i in self.velocity)}'
        if self.angular_velocity:
            s += f',\n\tangular velocity, {", ".join(str(i) for i in self.angular_velocity)}'
        if self.acceleration:
            s += f',\n\tacceleration, {", ".join(str(i) for i in self.acceleration)}'
        if self.angular_acceleration:
            s += f',\n\tangular acceleration, {", ".join(str(i) for i in self.angular_acceleration)}'
        return s


class ControlData(MBEntity):
    '''
    This section is read by the manager of all the bulk simulation data, namely the nodes, the drivers and
    the elements. It is used to set some global parameters closely related to the behavior of these entities,
    to tailor the initial assembly of the joints in case of structural simulations, and to tell the data manager
    how many entities of every type it should expect from the following sections. Historically this is due to
    the fact that the data structure for nodes and elements is allocated at the beginning with ﬁxed size. This
    is going to change, giving raise to a "free" and resizeable structure. But this practice is to be considered
    reliable since it allows a sort of double-check on the entities that are inserted.
    '''

    use_auto_differentiation: Optional[Union[bool, int]] = False    # 0 / 1 / True / False
    skip_initial_joint_assembly: Optional[Union[bool, int]] = False    # 0 / 1 / True / False
    simulation_title: Optional[str] = None
    print: Optional[Print] = None
    output_frequency: Optional[Union[int, MBVar]] = None
    output_meter: Optional[DriveCaller] = None
    output_results: Optional[OutputResults] = None
    default_orientation: Union[Literal["euler123", "euler313", "euler321", "orientation vector", "orientation matrix"]] = "euler123"
    model: Literal["static"] = "static"
    rbk_data: Optional[Union[ConstRBK, DriveRBK]] = None

    ## Model Counter Cards
    # Nodes
    abstract_nodes: Optional[Union[int, str]] = None
    electric_nodes: Optional[Union[int, str]] = None
    hydraulic_nodes: Optional[Union[int, str]] = None
    parameter_nodes: Optional[Union[int, str]] = None
    structural_nodes: Optional[Union[int, str]] = None
    thermal_nodes: Optional[Union[int, str]] = None

    # Drivers
    file_drivers: Optional[Union[int, str]] = None

    # Elements
    aerodynamic_elements: Optional[Union[int, str]] = None
    aeromodals: Optional[Union[int, str]] = None
    air_properties: Optional[Union[int, str]] = None
    automatic_structural_elements: Optional[Union[int, str]] = None
    beams: Optional[Union[int, str]] = None
    bulk_elements: Optional[Union[int, str]] = None
    electric_bulk_elements: Optional[Union[int, str]] = None
    electric_elements: Optional[Union[int, str]] = None
    external_elements: Optional[Union[int, str]] = None
    forces: Optional[Union[int, str]] = None
    genels: Optional[Union[int, str]] = None
    gravity: Optional[Union[int, str]] = None
    hydraulic_elements: Optional[Union[int, str]] = None
    induced_velocity_elements: Optional[Union[int, str]] = None
    joints: Optional[Union[int, str]] = None
    joint_regularizations: Optional[Union[int, str]] = None
    loadable_elements: Optional[Union[int, str]] = None
    output_elements: Optional[Union[int, str]] = None
    plates: Optional[Union[int, str]] = None # Not present in the manual, but present in CrankPanel_v2.mbd
    solids: Optional[Union[int, str]] = None
    surface_loads: Optional[Union[int, str]] = None
    rigid_bodies: Optional[Union[int, str]] = None

    @field_validator('use_auto_differentiation', mode='after')
    def set_use_auto_differentiation(cls, v):
        if isinstance(v, (int, bool)):  # Check if v is an int or a bool
            if v in [0, False]:
                return None  # Return None for 0 or False
            elif v in [1, True]:
                return "use auto differentiation"  # Return specific string for 1 or True
            else:
                raise ValueError("use_auto_differentiation must be 0, 1, True, or False.")
        else:
            raise TypeError("use_auto_differentiation must be of type int or bool.")

    @field_validator('skip_initial_joint_assembly', mode='after')
    def set_skip_initial_joint_assembly(cls, v):
        if isinstance(v, (int, bool)):  # Check if v is an int or a bool
            if v in [0, False]:
                return None  # Return None for 0 or False
            elif v in [1, True]:
                return "skip initial joint assembly"  # Return specific string for 1 or True
            else:
                raise ValueError("skip_initial_joint_assembly must be 0, 1, True, or False.")
        else:
            raise TypeError("skip_initial_joint_assembly must be of type int or bool.")

    def __str__(self):
        s = 'begin: control data;\n'
        if self.use_auto_differentiation:
            s += f'\tuse automatic differentiation;\n'
        if self.skip_initial_joint_assembly:
            s += f'\tskip initial joint assembly;\n'
        if self.simulation_title:
            s += f'\ttitle: {self.simulation_title};\n'
        if self.print:
            s += f'\t{self.print};\n'
        if self.output_frequency:
            s += f'\toutput frequency: {self.output_frequency};\n'
        if self.output_meter:
            s += f'\toutput meter: {self.output_meter};\n'
        if self.output_results:
            s += f'\t{self.output_results};\n'

        s += f'\tdefault orientation: {self.default_orientation};\n'
        s += f'\tmodel: {self.model};\n'

        if self.rbk_data:
            s += f'\trigid body kinematics: {self.rbk_data};\n'

        # Model Counter Cards - Nodes
        if self.abstract_nodes:
            s += f'\tabstract nodes: {self.abstract_nodes};\n'
        if self.electric_nodes:
            s += f'\telectric nodes: {self.electric_nodes};\n'
        if self.hydraulic_nodes:
            s += f'\thydraulic nodes: {self.hydraulic_nodes};\n'
        if self.parameter_nodes:
            s += f'\tparameter nodes: {self.parameter_nodes};\n'
        if self.structural_nodes:
            s += f'\tstructural nodes: {self.structural_nodes};\n'
        if self.thermal_nodes:
            s += f'\tthermal nodes: {self.thermal_nodes};\n'
        
        # Model Counter Cards - Drivers
        if self.file_drivers:
            s += f'\tfile drivers: {self.file_drivers};\n'
        
        # Model Counter Cards - Elements
        if self.aerodynamic_elements:
            s += f'\taerodynamic elements: {self.aerodynamic_elements};\n'
        if self.aeromodals:
            s += f'\taeromodals: {self.aeromodals};\n'
        if self.air_properties:
            s += f'\tair properties: {self.air_properties};\n'
        if self.automatic_structural_elements:
            s += f'\tautomatic structural elements: {self.automatic_structural_elements};\n'
        if self.beams:
            s += f'\tbeams: {self.beams};\n'
        if self.bulk_elements:
            s += f'\tbulk elements: {self.bulk_elements};\n'
        if self.electric_bulk_elements:
            s += f'\telectric bulk elements: {self.electric_bulk_elements};\n'
        if self.electric_elements:
            s += f'\telectric elements: {self.electric_elements};\n'
        if self.external_elements:
            s += f'\texternal elements: {self.external_elements};\n'
        if self.forces:
            s += f'\tforces: {self.forces};\n'
        if self.genels:
            s += f'\tgenels: {self.genels};\n'
        if self.gravity:
            s += f'\tgravity: {self.gravity};\n'
        if self.hydraulic_elements:
            s += f'\thydraulic elements: {self.hydraulic_elements};\n'
        if self.induced_velocity_elements:
            s += f'\tinduced velocity elements: {self.induced_velocity_elements};\n'
        if self.joints:
            s += f'\tjoints: {self.joints};\n'
        if self.joint_regularizations:
            s += f'\tjoint regularizations: {self.joint_regularizations};\n'
        if self.loadable_elements:
            s += f'\tloadable elements: {self.loadable_elements};\n'
        if self.output_elements:
            s += f'\toutput elements: {self.output_elements};\n'
        if self.plates:   # Not present in the manual, but present in CrankPanel_v2.mbd
            s += f'\tplates: {self.plates};\n'
        if self.solids:
            s += f'\tsolids: {self.solids};\n'
        if self.surface_loads:
            s += f'\tsurface loads: {self.surface_loads};\n'
        if self.rigid_bodies:
            s += f'\trigid bodies: {self.rigid_bodies};\n'

        s += 'end: control data;\n\n'
        return s
