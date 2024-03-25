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

from __future__ import print_function, division
import sys
from numbers import Number, Integral
import pdb

if sys.version_info[0] < 3:
        import __builtin__ as builtins
else:
        import builtins

declared_ConstMBVars = {}
declared_IfndefMBVars = {}
declared_MBVars = {}

MBDynLib_simplify = True

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

class MBVar(terminal_expression):
    base_types = ('bool', 'integer', 'real', 'string')
    var_types = base_types + tuple(['const' + ' ' + ti for ti in base_types]) +\
                tuple(['ifndef' + ' ' + 'const' + ti for ti in base_types])
    def __init__(self, name, var_type, expression):
        assert(name)
        self.name = name
        self.var_type = var_type
        self.expression = expression
        assert (var_type in self.var_types), (
            '\n-------------------\nERROR:' + 
            ' MBVar: unknown variable type {}\n\t'.format(var_type) + 
            '\n-------------------\n'
        )
       #self.do_declare = do_declare
        #if self.do_declare:
        assert (name in declared_ConstMBVars) == False, (
            '\n-------------------\nERROR:' + 
            ' re-defining an already declared const variable:\n\t' + 
            var_type + ' ' + name + 
            '\n-------------------\n')
        assert (name in declared_IfndefMBVars) == False, (
            '\n-------------------\nERROR:' + 
            ' re-defining an already declared ifndef variable:\n\t' + 
            var_type + ' ' + name + 
            '\n-------------------\n')
        self.declare()
    def __get__(self):
        return get_value(self.expression)
    def __trunc__(self):
        y = self.__get__()
        assert isinstance(y, int), (
                'Error, __trunc__  required for expression \n\'' + 
                str(self) + 
                '\'\nof type ' + str(type(y)) +
                ' \n')
        return self.expression.__trunc__()
    def __index__(self):
        return self.expression.__trunc__()
    def __str__(self):
        return str(self.name)
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
                '\n-------------------\nERROR:' + 
                ' re-defining an already declared variable of type ' + str(declared_MBVars[self.name].var_type) + '\n' + 
                'with different type ' + str(self.var_type) +
                '\n-------------------\n')
            if ('string' in self.var_type):
                print('set: ' + self.name + ' = \"' + str(self.expression) + '\";')
            else:
                print('set: ' + self.name + ' = ' + str(self.expression) + ';')
        else:
            declared_MBVars[self.name] = self
            if ('string' in self.var_type):
                print('set: ' + self.var_type + ' ' + self.name + ' = \"' + str(self.expression) + '\";')
            else:
                print('set: ' + self.name + ' = ' + str(self.expression) + ';')
        #globals()[self.name] = self    
        #__builtins__[self.name] = self    
        setattr(builtins, self.name, self)


class ConstMBVar(MBVar):
    def __init__(self, name, var_type, value):
        MBVar.__init__(self, name, 'const ' + var_type, value)
    def declare(self):
        #assert self.do_declare == True, (
        #    '\n-------------------\nERROR:' +
        #    ' declaring either temporary '
        #    'or already declared variable:\n\t' + 
        #    self.var_type + ' ' + self.name + 
        #    '\n-------------------\n')
        MBVar.declare(self)
        #self.do_declare = False
        declared_ConstMBVars[self.name] = self


class IfndefMBVar(MBVar):
    def __init__(self, name, var_type, value):
        if name in declared_MBVars:
            pass
        else:
            MBVar.__init__(self, name, 'ifndef ' + var_type, value)

class null:
    def __str__(self):
        s = 'null'
        return s

class eye:
    def __str__(self):
        s = 'eye'
        return s

class Reference:
    def __init__(self, idx, pos, orient, vel, angvel):
        assert isinstance(pos, Position), (
            '\n-------------------\nERROR:'+
            ' the position of a reference must be ' +
            ' an instance of the Position class;' +
            '\n-------------------\n')
        assert isinstance(orient, Position), (
            '\n-------------------\nERROR:' +
            ' the orientation of a reference must be ' +
            ' an instance of the Position class;' +
            '\n-------------------\n')
        assert isinstance(vel, Position), (
            '\n-------------------\nERROR:' +
            ' the velocity of a reference must be ' +
            ' an instance of the Position class;' +
            '\n-------------------\n')
        assert isinstance(angvel, Position), (
            '\n-------------------\nERROR:' +
            ' the angulare velocity of a reference must be ' +
            ' an instance of the Position class;' +
            '\n-------------------\n')
        self.idx = idx
        self.position = pos
        self.orientation = orient
        self.velocity = vel
        self.angular_velocity = angvel
    def __str__(self):
        s = 'reference: '
        s = s + str(self.idx) + ', \n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity) + ';\n'
        return s

class Position:
    def __init__(self, ref, rel_pos):
        self.reference = ref
        if isinstance(rel_pos, list):
            self.relative_position = rel_pos
        else:
            self.relative_position = [rel_pos]
    def __str__(self):
        s = ''
        if self.reference != '':
            s = 'reference, ' + str(self.reference) + ', '
        s = s + ', '.join(str(i) for i in self.relative_position)
        return s
    def isnull(self):
        return (self.reference == '') and isinstance(self.relative_position[0], null)
    def iseye(self):
        return (self.reference == '') and isinstance(self.relative_position[0], eye)

class Node:
    def __init__(self, idx, pos, orient, vel, angular_vel, node_type = 'dynamic',
            scale = 'default', output = 'yes'):
        assert isinstance(pos, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial position of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(orient, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial orientation of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(vel, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial velocity of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(angular_vel, Position), (
            '\n-------------------\nERROR:' + 
            ' the initial angular velocity of a node must be ' +  
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert node_type in ('dynamic', 'static',), (
            '\n-------------------\nERROR:' + 
            ' unrecognised or unsupported node type;' + 
            '\n-------------------\n')
        self.idx = idx
        self.position = pos
        self.orientation = orient
        self.velocity = vel
        self.angular_velocity = angular_vel
        self.node_type = node_type
        self.scale = scale
        self.output = output
    def __str__(self):
        s = 'structural: ' + str(self.idx) + ', ' + str(self.node_type) + ',\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.orientation) + ',\n'
        s = s + '\t' + str(self.velocity) + ',\n'
        s = s + '\t' + str(self.angular_velocity)
        if self.scale != 'default':
            s = s + ',\n\tscale, ' + str(self.scale)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DynamicNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'dynamic')

class StaticNode(Node):
    def __init__(self, idx, pos, orient, vel, angular_vel):
        Node.__init__(self, idx, pos, orient, vel, angular_vel, 'static')

class DisplacementNode():
    def __init__(self, idx, pos, vel, node_type = 'dynamic',
            scale = 'default', output = 'yes'):
        self.idx = idx
        self.position = pos
        self.velocity = vel
        self.node_type = node_type
        self.scale = scale
        self.output = output
    def __str__(self):
        s = 'structural: ' + str(self.idx) + ', ' + str(self.node_type) + ' displacement,\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + str(self.velocity)
        if self.scale != 'default':
            s = s + ',\n\t scale, ' + str(self.scale)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DynamicDisplacementNode(DisplacementNode):
    def __init__(self, idx, pos, vel):
        DisplacementNode.__init__(self, idx, pos, vel, 'dynamic')

class StaticDisplacementNode(DisplacementNode):
    def __init__(self, idx, pos, vel):
        DisplacementNode.__init__(self, idx, pos, vel, 'static')

class PointMass:
    def __init__(self, idx, node, mass, output = 'yes'):
        self.idx = idx
        self.node = node
        self.mass = mass
        self.output = output
    def __str__(self):
        s = 'body: ' + str(self.idx) + ', ' + str(self.node) + ', ' + str(self.mass)
        if self.output != 'yes':
            s = s + ', output, ' + str(self.output)
        s = s + ';\n'
        return s

class Element:
    idx = -1

class Body(Element):
    def __init__(self, idx, node, mass, position, inertial_matrix, inertial = null,
            output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' +
            ' in defining a body, the center of mass relative position ' + 
            ' mass must be an instance of the Position class;' + 
            '\n-------------------\n')
        assert isinstance(inertial_matrix, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a body, the inertial matrix' + 
            ' must be a list;' + 
            '\n-------------------\n')
        self.idx = idx
        self.type = 'body'
        self.node = node
        self.mass = mass
        self.position = position
        self.inertial_matrix = inertial_matrix
        self.inertial = inertial
        self.output = output
    def __str__(self):
        s = 'body: ' + str(self.idx) + ', ' + str(self.node) + ',\n'
        s = s + '\t' + str(self.mass) + ',\n'
        s = s + '\t' + str(self.position) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.inertial_matrix) 
        if self.inertial != null:
            s = s + ',\n'
            if isinstance(self.inertial, list):
                s = s + ', '.join(str(i) for i in self.inertial_matrix)
            else:
                s = s + ', ' + self.inertial_matrix
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralForce(Element):
    def __init__(self, idx, node, ftype, position, force_drive, 
            force_orientation = [], moment_orientation = [],
            moment_drive = [], output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural force, the relative arm must be' +
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert ftype in {'absolute', 'follower', 'total'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural force: ' + str(ftype) + 
            '\n-------------------\n')
        if ftype == 'total':
            assert isinstance(force_orientation, Position), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total force, the force orientation ' +
                ' must be an instance of the Position class;' + 
                '\n-------------------\n')
            assert isinstance(moment_orientation, Position), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total force, the moment orientation ' +
                ' must be an instance of the Position class;' + 
                '\n-------------------\n')
        self.idx = idx
        self.type = 'force'
        self.node = node
        self.ftype = ftype
        self.position = position
        self.force_drive = force_drive
        self.force_orientation = force_orientation
        self.moment_orientation = moment_orientation
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'force: ' + str(self.idx) + ', ' + self.ftype
        s = s + ',\n\t' + str(self.node)
        s = s + ',\n\t\tposition, ' + str(self.position)
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation)
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation)
            s = s + ',\n\t\tforce, ' + ', '.join(str(i) for i in self.force_drive)
            s = s + ',\n\t\tmoment, ' + ', '.join(str(i) for i in self.moment_drive)
        else: # ftype = { absolute|follower }
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.force_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralInternalForce(Element):
    def __init__(self, idx, nodes, ftype, positions, force_drive, 
            force_orientation = [], moment_orientation = [],
            moment_drive = [], output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal force with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions) , (
            '\n-------------------\nERROR:' + 
            ' in defining a structural internal force all relative arms ' +
            ' must be instances of the Position class;' + 
            '\n-------------------\n')
        assert ftype in {'absolute', 'follower', 'total'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural internal force: ' + str(ftype) + 
            '\n-------------------\n')
        if ftype == 'total':
            assert all(isinstance(pos, Position) for pos in force_orientation), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total internal force all the ' +
                ' force orientations must be instances of the Position class;' + 
                '\n-------------------\n')
            assert all(isinstance(pos, Position) for pos in moment_orientation), (
                '\n-------------------\nERROR:' + 
                ' in defining a structural total internal force all the ' +
                ' moment orientations must be instances of the Position class;' + 
                '\n-------------------\n')
        self.idx = idx
        self.type = 'force'
        self.nodes = nodes
        self.ftype = ftype
        self.positions = positions
        self.force_drive = force_drive
        self.force_orientation = force_orientation
        self.moment_orientation = moment_orientation
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'force: ' + str(self.idx) + ', ' + self.ftype + ' internal'
        s = s + ',\n\t' + str(self.nodes[0])
        s = s + ',\n\t\tposition, ' + str(self.positions[0])
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation[0])
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation[0])
        s = s + ',\n\t' + str(self.nodes[1])
        s = s + ',\n\t\tposition, ' + str(self.positions[1])
        if self.ftype == 'total':
            s = s + ',\n\t\tforce orientation, ' + str(self.force_orientation[1])
            s = s + ',\n\t\tmoment orientation, ' + str(self.moment_orientation[1])
            s = s + ',\n\t\tforce, ' + ', '.join(str(i) for i in self.force_drive)
            s = s + ',\n\t\tmoment, ' + ', '.join(str(i) for i in self.moment_drive)
        else: # ftype = { absolute|follower }
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.force_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralCouple(Element):
    def __init__(self, idx, node, ctype, position, moment_drive, output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural couple, the relative arm must be' +
            ' an instance of the Position class;' + 
            '\n-------------------\n')
        assert ctype in {'absolute', 'follower'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural couple: ' + str(ctype) + 
            ';\n-------------------\n')
        self.idx = idx
        self.type = 'couple'
        self.node = node
        self.ctype = ctype
        self.position = position
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'couple: ' + str(self.idx) + ', ' + self.ctype
        s = s + ',\n\t' + str(self.node)
        if len(self.position):
            s = s + ',\n\t\tposition, ' + str(self.position)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.moment_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class StructuralInternalCouple(Element):
    def __init__(self, idx, nodes, ctype, positions, moment_drive, output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal couple with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(positions) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a structural internal couple with ' + str(len(positions)) + 
            ' relative positions (!= 2);' +
            '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions), (
            '\n-------------------\nERROR:' + 
            ' in defining a structural internal couple all the relative positions ' +
            ' must be instances of the Position class;' + 
            '\n-------------------\n')
        assert ctype in {'absolute', 'follower'}, (
            '\n-------------------\nERROR:' + 
            ' unrecognised type of structural internal couple: ' + str(ctype) + 
            '\n-------------------\n')
        self.idx = idx
        self.type = 'couple'
        self.nodes = nodes
        self.ctype = ctype
        self.positions = positions
        self.moment_drive = moment_drive
        self.output = output
    def __str__(self):
        s = 'couple: ' + str(self.idx) + ', ' + self.ctype + ' inernal'
        s = s + ',\n\t' + str(self.nodes[0])
        s = s + ',\n\t\tposition, ' + str(self.position[0])
        s = s + ',\n\t' + str(self.nodes[1])
        s = s + ',\n\t\tposition, ' + str(self.position[1])
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.moment_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class Clamp(Element):
    def __init__(self, idx, node, pos = Position('', 'node'), 
            orient = Position('', 'node'), output = 'yes'):
        self.idx = idx
        self.type = 'joint'
        self.node = node
        self.position = pos
        self.orientation = orient
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', clamp, ' + str(self.node) + ',\n'
        s = s + '\tposition, ' + str(self.position) + ',\n'
        s = s + '\torientation, ' + str(self.orientation)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class TotalJoint(Element):
    def __init__(self, idx, nodes, positions, \
            position_orientations, rotation_orientations, \
            position_constraints, orientation_constraints, \
            position_drive, orientation_drive,
            output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' nodes must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative positions must be given in a list' + 
            '\n-------------------\n')    
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(position_orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative position orientations must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == len(position_orientations), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(position_orientations)) + ' position orientations;\n' +
            '\n-------------------\n')
        assert isinstance(rotation_orientations, list), (
            '\n-------------------\nERROR:' + 
            ' in defining a total joint, the' +
            ' relative rotation orientations must be given in a list' + 
            '\n-------------------\n')
        assert len(nodes) == len(rotation_orientations), (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(rotation_orientations)) + ' rotation orientations;\n' +
            '\n-------------------\n')
        assert isinstance(position_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' position constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(position_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' +
            str(len(position_constraints)) + ' position constraints;\n' +
            '\n-------------------\n')
        assert isinstance(orientation_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' orientation constraints must be given as a list;' + 
            '\n-------------------\n')    
        assert len(orientation_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' +
            str(len(orientation_constraints)) + ' orientation constraints;\n' +
            '\n-------------------\n')
        assert all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint all offsets must be instances of ' + 
            ' the class Position;\n' +
            '\n-------------------\n')
        self.idx = idx
        self.type = 'joint'
        self.nodes = nodes
        self.positions = positions
        self.position_orientations = position_orientations
        self.rotation_orientations = rotation_orientations
        self.position_constraints = position_constraints
        self.orientation_constraints = orientation_constraints
        self.position_drive = position_drive
        self.orientation_drive = orientation_drive
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', total joint'
        for (node, pos, pos_or, rot_or) in zip(self.nodes, self.positions,
                self.position_orientations, self.rotation_orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s = s + ',\n\t\tposition, ' + str(pos)
            if not(pos_or.iseye()):
                s = s + ',\n\t\tposition orientation, ' + str(pos_or)
            if not(rot_or.iseye()):
                s = s + ',\n\t\trotation orientation, ' + str(rot_or)
        if sum(self.position_constraints):
            s = s + ',\n\tposition constraint, '\
                    + ', '.join(str(pc) for pc in self.position_constraints)
            if isinstance(self.position_drive, list):
                s = s + ',\n\t\t' + ', '.join(str(i) for i in self.position_drive)
            else:
                s = s + ',\n\t\t' + str(self.position_drive)
        if sum(self.orientation_constraints):
            s = s + ',\n\torientation constraint, '\
                    + ', '.join(str(oc) for oc in self.orientation_constraints)
            if isinstance(self.orientation_drive, list):
                s = s + ',\n\t\t' + ', '.join(str(i) for i in self.orientation_drive)
            else:
                s = s + ',\n\t\t', + str(self.orientation_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s


class TotalPinJoint(Element):
    def __init__(self, idx, node, 
            positions, position_orientations, rotation_orientations, 
            position_constraints, orientation_constraints, 
            position_drive, orientation_drive,
            output = 'yes'):
        if not isinstance(positions, list):
            positions = [positions]
        assert (len(positions) in [1, 2]) and all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative positions must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        if not isinstance(position_orientations, list):
            position_orientations = [position_orientations]
        assert ((len(position_orientations) in [1, 2]) and all([isinstance(pos, Position) for pos in position_orientations])), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative position orientations must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        if not isinstance(rotation_orientations, list):
            rotation_orientations = [rotation_orientations]
        assert ((len(rotation_orientations) in [1, 2]) and all([isinstance(pos, Position) for pos in rotation_orientations])), (
            '\n-------------------\nERROR:' +
            ' in defining a total pin joint, ' + 
            ' relative rotation orientations must be given as a single instance' + 
            ' of the Position class or as a list of Position instances' + 
            '\n-------------------\n')
        assert isinstance(position_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' position constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(position_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(position_constraints)) + 
            ' position constraints;' + '\n-------------------\n')
        assert isinstance(orientation_constraints, list), (
            '\n-------------------\nERROR:' +
            ' in defining a total joint, ' 
            ' orientation constraints must be given as a list;' + 
            '\n-------------------\n')
        assert len(orientation_constraints) == 3, (
            '\n-------------------\nERROR:' +
            ' defining a total joint with ' + str(len(orientation_constraints)) + 
            ' orientation constraints;' + '\n-------------------\n')
        self.idx = idx
        self.type = 'joint'
        self.node = node
        self.positions = positions
        self.position_orientations = position_orientations
        self.rotation_orientations = rotation_orientations
        self.position_constraints = position_constraints
        self.orientation_constraints = orientation_constraints
        self.position_drive = position_drive
        self.orientation_drive = orientation_drive
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', total pin joint'
        s = s + ',\n\t' + str(self.node)
        if not(self.positions[0].isnull()):
            s = s + ',\n\t\tposition, ' + str(self.positions[0])
        if not(self.position_orientations[0].iseye()):
            s = s + ',\n\t\tposition orientation, ' + str(self.position_orientations[0])
        if not(self.rotation_orientations[0].iseye()):
            s = s + ',\n\t\trotation orientation, ' + str(self.rotation_orientations[0])
        if len(self.positions) == 2 and not(self.positions[1].isnull()):
            s = s + ',\n\t# GROUND'
            s = s + '\n\t\tposition, ' + str(self.positions[1])
        if len(self.position_orientations) == 2 and not(self.position_orientations[1].iseye()):
            s = s + ',\n\t\tposition orientation, ' + str(self.position_orientations[1])
        if len(self.rotation_orientations) == 2 and not(self.rotation_orientations[1].iseye()):
            s = s + ',\n\t\trotation orientation, ' + str(self.rotation_orientations[1])
        if sum(self.position_constraints):
            s = s + ',\n\tposition constraint, '\
                    + ', '.join(str(pc) for pc in self.position_constraints)
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.position_drive)
        if sum(self.orientation_constraints):
            s = s + ',\n\torientation constraint, '\
                    + ', '.join(str(oc) for oc in self.orientation_constraints)
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.orientation_drive)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class JointRegularization(Element):
    def __init__(self, idx, coefficients):
        assert (isinstance(coefficients, list) and len(coefficients) >= 1) or (isinstance(coefficients, Number)), (
            '\n-------------------\nERROR:' + 
            ' joint regularization needs at least one' +
            ' coefficient ' + '\n-------------------\n')
        self.idx = idx
        self.type = 'joint regularization'
        self.coefficients = coefficients
    def __str__(self):
        s = 'joint regularization: ' + str(self.idx) + ", tikhonov"
        if isinstance(self.coefficients, list):
            s = s + 'list, ' + ', '.join(str(co) for co in self.coefficients)
        else:
            s = s + ', ' + str(self.coefficients)
        s = s + ';\n'
        return s


class Rod(Element):
    def __init__(self, idx, nodes, positions, const_law, length = 'from nodes', 
            output = 'yes'):
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a rod with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a rod with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        if not isinstance(positions, list):
            positions = [positions]
        assert all([isinstance(pos, Position) for pos in positions]), (
            '\n-------------------\nERROR:' +
            ' in defining a rod all offsets must be instances of ' + 
            ' the class Position;\n' +
            '\n-------------------\n')
        self.idx = idx
        self.type = 'joint'
        self.nodes = nodes
        self.positions = positions
        self.const_law = const_law
        self.length = length
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', rod'
        for (node, position) in zip(self.nodes, self.positions):
            s = s + ',\n\t' + str(node)
            if not(position.isnull()):
                s = s + ',\n\t\tposition, ' + str(position)
        s = s + ',\n\t' + str(self.length) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.const_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DeformableDiaplacement(Element):
    def __init__(self, idx, nodes, positions, orientations, const_law, output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' +
            ' in defining a deformable displacement joint, the' +
            ' nodes must be given in a list' +
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' +
            ' defining a deformable displacement joint with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' +
            ' in defining a deformable displacement joint, the' +
            ' relative positions must be given in a list' +
            '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a deformable displacement joint with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(orientations, list), (
            '\n-------------------\nERROR:' +
            ' in defining a deformable displacement joint, the' +
            ' relative position orientations must be given in a list' +
            '\n-------------------\n')
        self.idx = idx
        self.type = 'joint'
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.constitutive_law = const_law
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', deformable displacement'
        for (node, pos, orient) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s + s + ',\n\t\tposition, ' + str(pos)
            if not(self.pos_or.iseye()):
                s + s + ',\n\t\torientation, ' + str(orient)
        s = s + '\n\t'
        if isinstance(self.constitutive_law, str):
            s = s + self.constitutive_law
        else:
            s = s + ', '.join(str(i) for i in self.constitutive_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class DeformableHinge(Element):
    def __init__(self, idx, nodes, positions, orientations, const_law, output = 'yes'):
        assert isinstance(nodes, list), (
            '\n-------------------\nERROR:' +
            ' in defining a deformable hinge, the' +
            ' nodes must be given in a list' +
            '\n-------------------\n')
        assert len(nodes) == 2, (
            '\n-------------------\nERROR:' +
            ' defining a deformable hinge with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert isinstance(positions, list), (
            '\n-------------------\nERROR:' +
            ' in defining a displacement hinge, the' +
            ' relative positions must be given in a list' +
            '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a deformable hinge with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert isinstance(orientations, list), (
            '\n-------------------\nERROR:' +
            ' in defining a deformable hinge, the' +
            ' relative position orientations must be given in a list' +
            '\n-------------------\n')
        self.idx = idx
        self.type = 'joint'
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.constitutive_law = const_law
        self.output = output
    def __str__(self):
        s = 'joint: ' + str(self.idx) + ', deformable hinge'
        for (node, pos, orient) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node)
            if not(pos.isnull()):
                s + s + ',\n\t\tposition, ' + str(pos)
            if not(orient.iseye()):
                s + s + ',\n\t\torientation, ' + str(orient)
        s = s + ',\n\t'
        if isinstance(self.constitutive_law, str):
            s = s + self.constitutive_law
        else:
            s = s + ', '.join(str(i) for i in self.constitutive_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class Shell(Element):
    def __init__(self, shell_type, idx, nodes, const_law, output = 'yes'):
        self.idx = idx
        self.type = shell_type
        self.nodes = nodes
        if isinstance(const_law, list):
            self.const_law = const_law
        else:
            self.const_law = [const_law]
        self.output = output
    def __str__(self):
        s = str(self.type) + ': ' + str(self.idx) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.nodes) + ',\n'
        s = s + '\t' + ', '.join(str(i) for i in self.const_law)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s
        
class Beam(Element):
    def __init__(self, idx, nodes, positions, orientations, const_laws_orientations,
            const_laws, output = 'yes'):
        assert len(nodes) == 3 or len(nodes) == 2, (
            '\n-------------------\nERROR:' + 
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes' + '\n-------------------\n')
        assert len(nodes) == len(positions), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(positions)) + ' relative positions;\n' +
            '\n-------------------\n')
        assert len(nodes) == len(orientations), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(nodes)) +
            ' nodes and ' + str(len(orientations)) + ' relative orientations;\n' +
            '\n-------------------\n')
        assert len(const_laws_orientations) == len(const_laws), (
            '\n-------------------\nERROR:' +
            ' defining a beam with ' + str(len(const_laws)) +
            ' coonstitutive laws and ' + str(len(const_laws_orientations)) + ' constitutive law orientations;' +
            '\n-------------------\n')
        if len(nodes) == 2:
            self.type = 'beam2'
        else:
            self.type = 'beam3'
        self.idx = idx
        self.nodes = nodes
        self.positions = positions
        self.orientations = orientations
        self.const_laws_orientations = const_laws_orientations
        self.const_laws = const_laws
        self.output = output
    def __str__(self):
        s = str(self.type) + ': ' + str(self.idx)
        for (node, position, orientation) in zip(self.nodes, self.positions, self.orientations):
            s = s + ',\n\t' + str(node) + ',\n\t\tposition, ' + str(position) + ',\n\t\torientation, ' + str(orientation)
        for (cl_or, cl) in zip(self.const_laws_orientations, self.const_laws):
            s = s + ',\n\t' + str(cl_or) + ',\n\t' 
            if isinstance(cl, str):
                s = s + cl 
            else: 
                s  = s + ', '.join(str(i) for i in cl)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        s = s + ';\n'
        return s

class AerodynamicBody(Element):
    def __init__(self, idx, node, 
            position, orientation, span,
            chord, aero_center, b_c_point, twist, integration_points,
            induced_velocity = [], tip_loss = [], control = [], 
            airfoil_data = [], unsteady = [], 
            jacobian = 'no', custom_output = [], output = 'yes'):
        assert isinstance(position, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the' + 
            ' relative surface offset must be an instance of the' + 
            ' Position class;' + '\n-------------------\n')
        assert isinstance(orientation, Position), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' relative surface orientation must be an instance of the' 
            ' Position class;' + '\n-------------------\n')
        assert isinstance(span, (Number)) or (isinstance(span, MBVar) and (span.var_type in ('real', 'const real'))), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body, the' + 
            ' surface span must be numeric' + 
            '\n-------------------\n')
        assert (isinstance(integration_points, Integral) and (integration_points > 0)) \
                or (isinstance(integration_points, MBVar) and integration_points.var_type in ('integer', 'const integer')), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with ' + str(integration_points) +
            ' integration_points' + '\n-------------------\n')
        assert (induced_velocity == []) or isinstance(induced_velocity, (Integral, MBVar)), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' induced velocity elment tag must be an integer or MBVar;' 
            '\n-------------------\n')
        assert not(len(unsteady)) or ((len(unsteady) > 0)*'bielawa' == 'bielawa'), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with unrecognised unsteady flag'
            '\n-------------------\n')
        assert (jacobian in {'yes', 'no'}) or isinstance(jacobian, bool), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic body with unrecognised jacobian flag'
            '\n-------------------\n')
        self.idx = idx
        self.type = 'aerodynamic body'
        self.node = node
        self.position = position
        self.orientation = orientation
        self.span = span
        self.chord = chord
        self.aero_center = aero_center
        self.b_c_point = b_c_point
        self.twist = twist
        self.integration_points = integration_points
        self.induced_velocity = induced_velocity
        self.tip_loss = tip_loss
        self.control = control
        self.airfoil_data = airfoil_data
        self.unsteady = unsteady
        self.jacobian = jacobian
        self.custom_output = custom_output
        self.output = output
    def __str__(self):
        s = 'aerodynamic body: ' + str(self.idx)
        s = s + ',\n\t ' + str(self.node)
        if self.induced_velocity:
            s = s + ',\n\t\tinduced velocity ' + str(self.induced_velocity)
        s = s + ',\n\t\t' + str(self.position)
        s = s + ',\n\t\t' + str(self.orientation)
        s = s + ',\n\t\t' + str(self.span)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.chord)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.aero_center)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.b_c_point)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.twist)
        if len(self.tip_loss):
            s = s + ',\n\t\ttip loss, ' + ', '.join(str(i) for i in self.tip_loss)
        s = s + '\n\t\t' + str(self.integration_points)
        if len(self.control):
            s = s + ',\n\t\tcontrol, ' + ', '.join(str(i) for i in self.control)
        if len(self.airfoil_data):
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.airfoil_data)
        if len(self.unsteady):
            s = s + ',\n\t\tunsteady, ' + str(self.unsteady)
        if len(self.jacobian):
            s = s + ',\n\t\tjacobian, ' + str(self.jacobian)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        if len(self.custom_output):
            s = s + ',\n\tcustom output, ' + ', '.join(str(i) for i in self.custom_output)
        s = s + ';\n'
        return s


class AerodynamicBeam(Element):
    def __init__(self, idx, beam, 
            positions, orientations,
            chord, aero_center, b_c_point, twist, integration_points, 
            induced_velocity = [], tip_loss = [], control = [], 
            airfoil_data = [], unsteady = [], 
            jacobian = 'no', custom_output = [], output = 'yes'):
        assert len(positions) in {2,3}, (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(len(positions)) +
            ' relative surface offsets (not in [2,3])' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in positions), (
            ' in defining an aerodynamic beam the' + 
            ' relative surface offsets must be instances of the' + 
            ' Position class;' + '\n-------------------\n')
        assert len(orientations) in {2,3}, (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(len(orientations)) +
            ' relative surface orientations (not in [2,3])' + '\n-------------------\n')
        assert all(isinstance(pos, Position) for pos in orientations), (
            ' in defining an aerodynamic beam the' + 
            ' relative surface orientations must be instances of the' + 
            ' Position class;' + '\n-------------------\n')
        assert len(positions) == len(orientations), (
            '\n-------------------\nERROR:' + 
            ' definining an aerodynamic beam with ' + str(len(positions)) + 
            ' relative surface offsets and ' + str(len(orientations)) + 
            ' relative surface orientations' + '\n-------------------\n')
        assert (isinstance(integration_points, Integral) and (integration_points > 0))\
                or (isinstance(integration_points, MBVar) and integration_points.var_type in ('integer', 'const integer')), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with ' + str(integration_points) +
            ' integration_points' + '\n-------------------\n')
        assert (induced_velocity == []) or isinstance(induced_velocity, (Integral, MBVar)), (
            '\n-------------------\nERROR:' + 
            ' in defining an aerodynamic body the '
            ' induced velocity elment tag must be an integer or an MBVar;' 
            '\n-------------------\n')
        assert not(len(unsteady)) or ((len(unsteady) > 0)*'bielawa' == 'bielawa'), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with unrecognised unsteady flag'
            '\n-------------------\n')
        assert (jacobian in {'yes', 'no'}) or isinstance(jacobian, bool), (
            '\n-------------------\nERROR:' + 
            ' defining an aerodynamic beam with unrecognised jacobian flag'
            '\n-------------------\n')
        self.idx = idx
        self.type = 'aerodynamic beam' + str(len(self.positions))
        self.beam = beam
        self.positions = positions
        self.orientations = orientations
        self.chord = chord
        self.aero_center = aero_center
        self.b_c_point = b_c_point
        self.twist = twist
        self.integration_points = integration_points
        self.induced_velocity = induced_velocity
        self.tip_loss = tip_loss
        self.control = control
        self.airfoil_data = airfoil_data
        self.unsteady = unsteady
        self.jacobian = jacobian
        self.custom_output = custom_output
        self.output = output
    def __str__(self):
        s = 'aerodynamic beam' + str(len(self.positions)) + ': ' + str(self.idx)
        s = s + ',\n\t ' + str(self.beam)
        if self.induced_velocity:
            s = s + ',\n\t\tinduced velocity ' + str(self.induced_velocity)
        for (pos, ori) in zip(self.positions, self.orientations):
            s = s + ',\n\t\t' + str(pos)
            s = s + ',\n\t\t' + str(ori)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.chord)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.aero_center)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.b_c_point)
        s = s + ',\n\t\t' + ', '.join(str(i) for i in self.twist)
        if len(self.tip_loss):
            s = s + ',\n\t\ttip loss, ' + ', '.join(str(i) for i in self.tip_loss)
        s = s + ',\n\t\t' + str(self.integration_points)
        if len(self.control):
            s = s + ',\n\t\tcontrol, ' + ', '.join(str(i) for i in self.control)
        if len(self.airfoil_data):
            s = s + ',\n\t\t' + ', '.join(str(i) for i in self.airfoil_data)
        if len(self.unsteady):
            s = s + ',\n\t\tunsteady, ' + str(self.unsteady)
        if self.jacobian == 'yes':
            s = s + ',\n\t\tjacobian, ' + str(self.jacobian)
        if self.output != 'yes':
            s = s + ',\n\toutput, ' + str(self.output)
        if len(self.custom_output):
            s = s + ',\n\tcustom output, ' + ', '.join(str(i) for i in self.custom_output)
        s = s + ';\n'
        return s

# General stuff
class NodeDof:
    idx = -1
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['node_label'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <node_label> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.node_label = kwargs['node_label']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <node_label> must be provided' + 
                    '\n-------------------\n'
            )
        try:
            if kwargs['node_type'] not in ('abstract', 'electric', 'hydraulic', 'parameter', 'structural', 'thermal'):
                raise ValueError(
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <node_type> must be either \'abstract\', \'electric\'' +
                    ' \'hydraulic\', \'parameter\', \'structural\', \'thermal\'' +
                    '\n-------------------\n'
                    )
            else:
                self.node_type = kwargs['node_type']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <node_type> must be provided' +
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['dof_number'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <dof_number> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.dof_number = kwargs['dof_number']
        except KeyError:
            pass
        try:
            if kwargs['dof_order'] not in ('algebraic', 'differential'):
                raise ValueError(
                    '\n-------------------\nERROR:' +
                    ' NodeDof: <dof_order> must either be an integer value or an MBVar' + 
                    '\n-------------------\n'
                    )
            else:
                self.dof_order = kwargs['dof_order']
        except KeyError:
            pass
    def __str__(self):
        s = '{}, {}'.format(self.node_label, self.node_type)
        if hasattr(self, 'dof_number'):
            s = s + ', {}'.format(self.dof_number)
        if hasattr(self, 'dof_order'):
            s = s + ', {}'.format(self.dof_order)
        return s

# Drives
class DriveCaller:
    idx = -1
    pass

class ArrayDriveCaller(DriveCaller):
    type = 'array'
    def __init__(self, *args, **kwargs):
        for arg in args:
            assert isinstance(arg, DriveCaller), (
                    '\n-------------------\nERROR:' +
                    ' ArrayDriveCaller: each argument of constructor must be' + 
                    ' a DriveCaller instance' + 
                    '\n-------------------\n')
        self.drives = args
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ArrayDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{},'.format(self.type)
        for drive in self.drives:
            if drive.idx < 0:
                s = s + '\n\t{}'.format(drive)
            else:
                s = s + '\n\treference, {}'.format(drive.idx)
        return s

class BistopDriveCaller(DriveCaller):
    type = 'bistop'
    def __init__(self, **kwargs):
        try:
            if kwargs['initial_status'] in ['active', 'inactive']:
                self.initial_status = kwargs['initial status']
            else:
                raise ValueError(
                        '\n------------------\nERROR:' + 
                        ' BistopDriveCaller: <initial_status> must be' + 
                        ' either \'active\' or \'inactive\'' + 
                        '\n------------------\n')
        except KeyError:
            self.initial_status = 'active'
            pass
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' BistopDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        assert isinstance(kwargs['activation_condition'], DriveCaller), (
                '\n-------------------\nERROR:' +
                ' BistopDriveCaller: <activation_condition> must be' + 
                ' a DriveCaller instance' + 
                '\n-------------------\n')

        assert isinstance(kwargs['deactivation_condition'], DriveCaller), (
                '\n-------------------\nERROR:' +
                ' BistopDriveCaller: <deactivation_condition> must be' + 
                ' a DriveCaller instance' + 
                '\n-------------------\n')
        self.activation_condition = kwargs['activation_condition']
        self.deactivation_condition = kwargs['deactivation_condition']
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{},\n\tinitial status, {},'.format(self.type, self.initial_status)
        s = s + '\n\t# activation condition drive'
        if self.activation_condition.idx < 0:
            s = s + '\n\t{}'.format(self.activation_condition)
        else:
            s = s + '\n\treference, {}'.format(self.activation_condition.idx)
        s = s + '\n\t# deactivation condition drive'
        if self.deactivation_condition.idx < 0:
            s = s + '\n\t{}'.format(self.deactivation_condition)
        else:
            s = s + '\n\treference, {}'.format(self.deactivation_condition.idx)
        return s


class ConstDriveCaller(DriveCaller):
    type = 'const'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ConstDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        
        try:
            assert isinstance(kwargs['const_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ConstDriveCaller: <const_value> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.const_value = kwargs['const_value']
        except KeyError:
            errprint(
                '\n-------------------\nERROR:' +
                ' ConstDriveCaller: <const_value> must be provided' + 
                '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx) 
        s = s + '{}, {}'.format(self.type, self.const_value)
        return s


class ClosestNextDriveCaller(DriveCaller):
    type = 'closest next'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_time = kwargs['initial_time']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' ClosestNextDriveCaller: <initial_time> is not set, assuming 0.' + 
                '\n-------------------\n')
            self.initial_time = 0.
        
        try:
            assert isinstance(kwargs['final_time'], (Number, MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <final_time> must either be a number, '
                    '\'forever\', or an MBVar' + 
                    '\n-------------------\n')
            self.final_time = kwargs['final_time']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <final_time> is not set' + 
                    '\n-------------------\n')
        try:
            assert isinstance(kwargs['increment'], DriveCaller), (
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <increment> should be a' + 
                    ' DriveCaller instance' + 
                    '\n-------------------\n')
            self.increment = kwargs['increment']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' ClosestNextDriveCaller: <increment> is not set' + 
                    '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ',\n\t{}, {},'.format(self.initial_time, self.final_time)
        s = s + '\n\t# increment drive'
        if self.increment.idx < 0:
            s = s + '\n\t{}'.format(self.increment)
        else:
            s = s + '\n\treference, {}'.format(self.increment.idx)
        return s


class CosineDriveCaller(DriveCaller):
    type = 'cosine'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_time = kwargs['initial_time']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CosineDriveCaller: <initial_time> not set, assuming 0.' + 
                    '\n-------------------\n'
                    )
            self.initial_time = 0.
            pass
        try:
            assert isinstance(kwargs['angular_velocity'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <angular_velocity> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.angular_velocity = kwargs['angular_velocity']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <angular_velocity> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['amplitude'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <amplitude> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.amplitude = kwargs['amplitude']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <amplitude> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['number_of_cycles'], (Number, MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <number_of_cycles> must either be a number,'
                    ' one in (\'half\', \'one\', \'forever\'), or an MBVar' + 
                    '\n-------------------\n')
            self.number_of_cycles = kwargs['number_of_cycles']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <number_of_cycles> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['initial_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CosineDriveCaller: <initial_value> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_value = kwargs['initial_value']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CosineDriveCaller: <initial_value> not provided, assuming 0.' + 
                    '\n-------------------\n'
            )
            self.initial_value = 0.
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}, {}, '.format(self.type, self.initial_time)
        s = s + '{}, {}, '.format(self.angular_velocity, self.amplitude)
        s = s + '{}, {}'.format(self.number_of_cycles, self.initial_value)
        return s

class CubicDriveCaller(DriveCaller):
    type = 'cubic'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CubicDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['const_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CubicDriveCaller: <const_coef> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.const_coef = kwargs['const_coef']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CubicDriveCaller: <const_coef> is required' + 
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['linear_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CubicDriveCaller: <linear_coef> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.linear_coef = kwargs['linear_coef']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CubicDriveCaller: <linear_coef> is required' + 
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['parabolic_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CubicDriveCaller: <parabolic_coef> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.parabolic_coef = kwargs['parabolic_coef']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CubicDriveCaller: <parabolic_coef> is required' + 
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['cubic_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' CubicDriveCaller: <cubic_coef> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.cubic_coef = kwargs['cubic_coef']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' CubicDriveCaller: <cubic_coef> is required' + 
                    '\n-------------------\n'
                    )
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}, {}, '.format(self.type, self.const_coef)
        s = s + '{}, {}, '.format(self.linear_coef, self.parabolic_coef)
        s = s + '{}'.format(self.cubic_coef)
        return s

class DirectDriveCaller(DriveCaller):
    type = 'direct'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DirectDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        return s

class DiscreteFilterDriveCaller(DriveCaller):
    type = 'discrete filter'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['n_a'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <n_a> must either be an integer or an MBVar' + 
                    '\n-------------------\n')
            self.n_a = kwargs['n_a']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' DiscreteFilterDriveCaller: <n_a> is required' + 
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['a'], list), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <a> must be a list of' + 
                    ' numbers or MBVars' + 
                    '\n-------------------\n')
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' DiscreteFilterDriveCaller: <n_a> is required' + 
                    '\n-------------------\n'
                    )
        for a_i in kwargs['a']:
            assert isinstance(a_i, (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: each component of <a> must be a' + 
                    ' number or an MBVar' + 
                    '\n-------------------\n'
                )
        self.a = kwargs['a']
        try:
            assert isinstance(kwargs['b_0'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <b_0> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.b_0 = kwargs['b_0']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' DiscreteFilterDriveCaller: <b_0> is required' +
                    ' set to 0 if not needed' +
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['n_b'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <n_b> must either be an integer or an MBVar' + 
                    '\n-------------------\n')
            self.n_b = kwargs['n_b']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' DiscreteFilterDriveCaller: <n_b> is required' + 
                    '\n-------------------\n'
                    )
        try:
            assert isinstance(kwargs['b'], list), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <b> must be a list of' + 
                    ' numbers or MBVars' + 
                    '\n-------------------\n')
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' DiscreteFilterDriveCaller: <n_b> is required' + 
                    '\n-------------------\n'
                    )
        for b_i in kwargs['b']:
            assert isinstance(b_i, (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: each component of <b> must be a' + 
                    ' number or an MBVar' + 
                    '\n-------------------\n'
                )
        self.b = kwargs['b']
        try:
            assert isinstance(kwargs['input_drive'], DriveCaller), (
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <input_drive> should be a' + 
                    ' DriveCaller instance' + 
                    '\n-------------------\n')
            self.input_drive = kwargs['input_drive']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' DiscreteFilterDriveCaller: <input_drive> is not set' + 
                    '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}, '.format(self.type)
        s = s + '\n\t{}'.format(self.n_a)
        for a_i in self.a:
            s = s + ',\n\t\t{}'.format(a_i)
        s = s + ',\n\t{}'.format(self.b_0)
        s = s + ',\n\t{}'.format(self.n_b)
        for b_i in self.b:
            s = s + ',\n\t\t{}'.format(b_i)
        if self.input_drive.idx >= 0:
            s = s + ',\n\treference, {}'.format(self.input_drive.idx)
        else:
            s = s + ',\n\t{}'.format(self.input_drive)
        return s


class DofDriveCaller(DriveCaller):
    type = 'dof'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DofDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n'
            )
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert(isinstance(kwargs['driving_dof'], NodeDof)), (
                    '\n-------------------\nERROR:' +
                    ' DofDriveCaller: <driving_dof> must be a NodeDof' + 
                    '\n-------------------\n'
            )
            self.driving_dof = kwargs['driving_dof']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' DofDriveCaller: <driving_dof> must be provided' + 
                    '\n-------------------\n'
            )
        try:
            assert(isinstance(kwargs['func_drive'], DriveCaller)), (
                    '\n-------------------\nERROR:' +
                    ' DofDriveCaller: <func_drive> must be a DriveCaller' + 
                    '\n-------------------\n'
            )
            self.func_drive = kwargs['func_drive']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' DofDriveCaller: <func_drive> must be provided' + 
                    '\n-------------------\n'
            )
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}'.format(self.idx)
        s = s + ', {}'.format(self.type)
        s = s + ',\n\t{}'.format(self.driving_dof)
        s = s + ',\n\t{}'.format(self.func_drive)
        return s

class DoubleRampDriveCaller(DriveCaller):
    type = 'double ramp'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['a_slope'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <a_slope> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.a_slope = kwargs['a_slope']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleRampDriveCaller: <a_slope> must be provided' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['a_initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <a_initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.a_initial_time = kwargs['a_initial_time']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleRampDriveCaller: <a_initial_time> is not set, assuming 0.' + 
                '\n-------------------\n')
            self.a_initial_time = 0.
            pass
        try:
            assert isinstance(kwargs['a_final_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <a_final_time> must either be a number'
                    ' or an MBVar' + 
                    '\n-------------------\n')
            self.a_final_time = kwargs['a_final_time']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <a_final_time> must be provided' + 
                    '\n-------------------\n')
        try:
            assert isinstance(kwargs['d_slope'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <d_slope> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.d_slope = kwargs['d_slope']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleRampDriveCaller: <d_slope> must be provided' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['d_initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <d_initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.d_initial_time = kwargs['d_initial_time']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleRampDriveCaller: <d_initial_time> must be provided' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['d_final_time'], (Number, MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <a_final_time> must either be a number,'
                    ' an MBVar, or \'forever\'' + 
                    '\n-------------------\n')
            self.d_final_time = kwargs['d_final_time']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <d_final_time> is not set' + 
                    '\n-------------------\n')
        try:
            assert isinstance(kwargs['initial_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleRampDriveCaller: <initial_value> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_value = kwargs['initial_value']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleRampDriveCaller: <initial_value> must be provided' + 
                '\n-------------------\n') # Why is it not assumed to be zero?
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ',\n\t{}, {}, {}'.format(self.a_slope, self.a_initial_time, self.a_final_time)
        s = s + ',\n\t{}, {}, {}'.format(self.d_slope, self.d_initial_time, self.d_final_time)
        s = s + ',\n\t{}'.format(self.initial_value)
        return s

class DoubleStepDriveCaller(DriveCaller):
    type = 'double step'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleStepDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleStepDriveCaller: <initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_time = kwargs['initial_time']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleStepDriveCaller: <initial_time> not set, assuming 0.' + 
                '\n-------------------\n')
            self.initial_time = 0.
            pass
        try:
            assert isinstance(kwargs['final_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleStepDriveCaller: <final_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.final_time = kwargs['final_time']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleStepDriveCaller: <final_time> is not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['step_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleStepDriveCaller: <step_value> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.step_value = kwargs['step_value']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleStepDriveCaller: <step_value> is not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['initial_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DoubleStepDriveCaller: <initial_value> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_value = kwargs['initial_value']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DoubleStepDriveCaller: <initial_value> is not set, assuming 0.' + 
                '\n-------------------\n')
            self.initial_value = 0.
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ',\n\t{}, {}'.format(self.initial_time, self.final_time)
        s = s + ',\n\t{}, {}'.format(self.step_value, self.initial_value)
        return s


class DriveDriveCaller(DriveCaller):
    type = 'drive'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' DriveDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['drive_caller1'], DriveCaller), (
                    '\n-------------------\nERROR:' +
                    ' DriveDriveCaller: <drive_caller1> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.drive_caller1 = kwargs['drive_caller1']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DriveDriveCaller: <drive_caller1> not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['drive_caller2'], DriveCaller), (
                    '\n-------------------\nERROR:' +
                    ' DriveDriveCaller: <drive_caller2> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.drive_caller2 = kwargs['drive_caller2']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' DriveDriveCaller: <drive_caller1> not set' + 
                '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        if self.drive_caller1.idx < 0:
            s = s + ',\n\t{}'.format(self.drive_caller1)
        else:
            s = s + ',\n\treference, {}'.format(self.drive_caller1.idx)
        if self.drive_caller2.idx < 0:
            s = s + ',\n\t{}'.format(self.drive_caller2)
        else:
            s = s + ',\n\treference, {}'.format(self.drive_caller2.idx)
        return s


class ElementDriveCaller(DriveCaller):
    type = 'element'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' ElementDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['element'], Element), (
                    '\n-------------------\nERROR:' +
                    ' ElementDriveCaller: <element> must either be an instance of Element' + 
                    '\n-------------------\n')
            self.element = kwargs['element']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' ElementDriveCaller: <element> not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['private_data'], str), (
                    '\n-------------------\nERROR:' +
                    ' ElementDriveCaller: <private_data> must either be a string' + 
                    '\n-------------------\n')
            self.private_data = kwargs['private_data']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' ElementDriveCaller: <final_time> is not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['func_drive'], (DriveCaller, str)), (
                    '\n-------------------\nERROR:' +
                    ' ElementDriveCaller: <func_drive> must either be a' +
                    ' DriveCaller or \'direct\'' + 
                    '\n-------------------\n')
            if isinstance(kwargs['func_drive'], str) and kwargs['func_drive'] != 'direct':
                raise ValueError(
                    '\n-------------------\nERROR:' +
                    ' ElementDriveCaller: <func_drive> must either be a' +
                    ' DriveCaller or \'direct\'' + 
                    '\n-------------------\n'
                    )
            self.func_drive = kwargs['func_drive']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' ElementDriveCaller: <func_drive> is not set' + 
                '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ', {}, {}'.format(self.element.idx, self.element.type)
        s = s + ', string, \"{}\"'.format(self.private_data)
        s = s + ', {}'.format(self.func_drive)
        return s


class ExponentialDriveCaller(DriveCaller):
    type = 'exponential'
    def __init__(self, **kwargs):
        try:
            arg = 'idx'
            assert isinstance(kwargs[arg], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) + 
                    ' an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs[arg]
        except KeyError:
            pass
        try:
            arg = 'amplitude_value'
            assert isinstance(kwargs[arg], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) +
                    ' a number of an MBVar' + 
                    '\n-------------------\n')
            self.amplitude_value = kwargs[arg]
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' {}: <{}> not set'.format(self.__class__.__name__, arg) + 
                '\n-------------------\n')
        try:
            arg = 'time_constant_value'
            assert isinstance(kwargs[arg], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) +
                    ' a number of an MBVar' + 
                    '\n-------------------\n')
            self.time_constant_value = kwargs[arg]
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' {}: <{}> not set'.format(self.__class__.__name__, arg) + 
                '\n-------------------\n')
        try:
            arg = 'initial_time'
            assert isinstance(kwargs[arg], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) +
                    ' a number of an MBVar' + 
                    '\n-------------------\n')
            self.initial_time = kwargs[arg]
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' {}: <{}> not set, assuming 0.'.format(self.__class__.__name__, arg) + 
                '\n-------------------\n')
            self.initial_time = 0.
            pass
        try:
            arg = 'initial_value'
            assert isinstance(kwargs[arg], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) +
                    ' a number of an MBVar' + 
                    '\n-------------------\n')
            self.initial_value = kwargs[arg]
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' {}: <{}> not set, assuming 0.'.format(self.__class__.__name__, arg) + 
                '\n-------------------\n')
            self.initial_value = 0.
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ', {}, {}, {}, {}'.format(
                    self.amplitude_value, 
                    self.time_constant_value,
                    self.initial_time,
                    self.initial_value
                    )
        return s


class FileDriveDrive(DriveCaller):
    # TODO: needs FileDrive before
    pass


class FourierSeriesDrive(DriveCaller):
    type = 'fourier series'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' FourierSeriesDrive: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            arg = 'initial_time'
            assert isinstance(kwargs[arg], (MBVar, Number)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a number or an MBVar' + 
                    '\n-------------------\n')
            if isinstance(kwargs[arg], MBVar) and ('real' not in kwargs[arg].var_type):
                raise TypeError(
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' an MBVar of type real' + 
                    '\n-------------------\n'
                )
            self.initial_time= kwargs[arg]
        except KeyError:
            pass
        try:
            arg = 'angular_velocity'
            assert isinstance(kwargs[arg], (MBVar, Number)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a number or an MBVar' + 
                    '\n-------------------\n')
            if isinstance(kwargs[arg], MBVar) and ('real' not in kwargs[arg].var_type):
                raise TypeError(
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' an MBVar of type real' + 
                    '\n-------------------\n'
                )
            self.angular_velocity = kwargs[arg]
        except KeyError:
            pass
        try:
            arg = 'number_of_terms'
            assert isinstance(kwargs[arg], (MBVar, Integral)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a number or an MBVar' + 
                    '\n-------------------\n')
            if isinstance(kwargs[arg], MBVar) and ('integer' not in kwargs[arg].var_type):
                raise TypeError(
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' an MBVar of type integer' + 
                    '\n-------------------\n'
                )
            self.number_of_terms = kwargs[arg]
        except KeyError:
            pass
        try:
            arg = 'number_of_cycles'
            assert isinstance(kwargs[arg], (MBVar, Integral, str)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a number or an MBVar' + 
                    '\n-------------------\n')
            if isinstance(kwargs[arg], MBVar) and kwargs[arg].var_type == 'string':
                if kwargs[arg] not in ('one', 'forever'):
                    raise ValueError(
                        '\n-------------------\nERROR:' +
                        ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                        ' either \'one\' or \'forever\', if of type string' + 
                        '\n-------------------\n'
                        )
            self.number_of_cycles = kwargs[arg]
        except KeyError:
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ', {}, {}, {}'.format(
                    self.initial_time, 
                    self.angular_velocity,
                    self.number_of_terms
                    )
        s = s + ',\n\t {}'.format(self.coefs)
        return s

class FrequencySweepDriveCaller(DriveCaller):
    type = 'frequency sweep'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n'
            )
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n')
            self.initial_time = kwargs['initial_time']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' FrequencySweepDriveCaller: <initial_time> not set, assuming 0.' + 
                    '\n-------------------\n'
                    )
            self.initial_time = 0.
            pass
        try:
            assert(isinstance(kwargs['angular_velocity_drive'], DriveCaller)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <angular_velocity_drive> must be a DriveCaller' + 
                    '\n-------------------\n'
            )
            self.angular_velocity_drive = kwargs['angular_velocity_drive']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <angular_velocity_drive> must be provided' + 
                    '\n-------------------\n'
            )
        try:
            assert(isinstance(kwargs['amplitude_drive'], DriveCaller)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <amplitude_drive> must be a DriveCaller' + 
                    '\n-------------------\n'
            )
            self.amplitude_drive = kwargs['amplitude_drive']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <amplitude_drive> must be provided' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['initial_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <initial_value> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.initial_value = kwargs['initial_value']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' FrequencySweepDriveCaller: <initial_value> not provided, assuming 0.' + 
                    '\n-------------------\n'
            )
            self.initial_value = 0.
            pass
        try:
            assert isinstance(kwargs['final_time'], (Number, MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <final_time> must either be a number, '
                    '\'forever\', or an MBVar' + 
                    '\n-------------------\n'
            )
            self.final_time = kwargs['final_time']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <final_time> is not set' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['final_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <final_value> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.final_value = kwargs['final_value']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' FrequencySweepDriveCaller: <final_value> is not set.' + 
                    '\n-------------------\n'
            )
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ',\n\t{}'.format(self.initial_time)
        s = s + ',\n\t{}'.format(self.angular_velocity_drive)
        s = s + ',\n\t{}'.format(self.amplitude_drive)
        s = s + ',\n\t{}, {}'.format(self.initial_value, self.final_time)
        s = s + ',\n\t{}'.format(self.final_value)
        return s

class GiNaCDriveCaller(DriveCaller):
    type = 'ginac'
    def __init__(self, **kwargs):
        try:
            arg = 'idx'
            assert isinstance(kwargs[arg], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must either be'.format(self.__class__.__name__, arg) + 
                    ' an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs[arg]
        except KeyError:
            pass
        try:
            arg = 'expression'
            assert isinstance(kwargs[arg], (MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a string or an MBVar' + 
                    '\n-------------------\n')
            if isinstance(kwargs[arg], MBVar) and ('string' not in kwargs[arg].var_type):
                raise TypeError(
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' an MBVar of type string' + 
                    '\n-------------------\n'
                )
            self.expression = kwargs[arg]
        except KeyError:
                errprint(
                '\n-------------------\nERROR' +
                ' {}: <{}> not set'.format(self.__class__.__name__, arg) + 
                '\n-------------------\n'
                )
        try:
            arg = 'symbol'
            assert isinstance(kwargs[arg], (MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' a string or an MBVar' + 
                    '\n-------------------\n')
            if isinstance (kwargs[arg], MBVar) and ('string' not in kwargs[arg].var_type):
                raise TypeError(
                    '\n-------------------\nERROR:' +
                    ' {}: <{}> must be'.format(self.__class__.__name__, arg) +
                    ' an MBVar of type string' + 
                    '\n-------------------\n'
                )
            self.symbol = kwargs[arg]
        except KeyError:
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        try:
            if isinstance(self.symbol, MBVar):
                s = s + ', symbol, {}'.format(self.symbol)
            else:
                s = s + ', symbol, \"{}\"'.format(self.symbol)
        except AttributeError:
            pass
        if isinstance(self.expression, MBVar):
            s = s + ', {}'.format(self.expression)
        else:
            s = s + ', \"{}\"'.format(self.expression)
        return s
    pass


class LinearDriveCaller(DriveCaller):
    type = 'linear'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' LinearDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n')
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['const_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' LinearDriveCaller: <const_coef> must either be a number of an MBVar' + 
                    '\n-------------------\n')
            self.const_coef = kwargs['const_coef']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' LinearDriveCaller: <const_coef> not set' + 
                '\n-------------------\n')
        try:
            assert isinstance(kwargs['slope_coef'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' LinearDriveCaller: <slope_coef> must either be a number of an MBVar' + 
                    '\n-------------------\n')
            self.slope_coef = kwargs['slope_coef']
        except KeyError:
            (
                '\n-------------------\nWARNING:' +
                ' LinearDriveCaller: <slope_coef> not set' + 
                '\n-------------------\n')
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}'.format(self.type)
        s = s + ', {}, {}'.format(self.const_coef, self.slope_coef)
        return s
    
class SineDriveCaller(DriveCaller):
    type = 'sine'
    def __init__(self, **kwargs):
        try:
            assert isinstance(kwargs['idx'], (Integral, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <idx> must either be an integer value or an MBVar' + 
                    '\n-------------------\n'
            )
            self.idx = kwargs['idx']
        except KeyError:
            pass
        try:
            assert isinstance(kwargs['initial_time'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <initial_time> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.initial_time = kwargs['initial_time']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' SineDriveCaller: <initial_time> not set, assuming 0.' + 
                    '\n-------------------\n'
            )
            self.initial_time = 0.
            pass
        try:
            assert isinstance(kwargs['angular_velocity'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <angular_velocity> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.angular_velocity = kwargs['angular_velocity']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <angular_velocity> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['amplitude'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <amplitude> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.amplitude = kwargs['amplitude']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <amplitude> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['number_of_cycles'], (Number, MBVar, str)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <number_of_cycles> must either be a number,'
                    ' one in (\'half\', \'one\', \'forever\'), or an MBVar' + 
                    '\n-------------------\n')
            self.number_of_cycles = kwargs['number_of_cycles']
        except KeyError:
            errprint(
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <amplitude> is required' + 
                    '\n-------------------\n'
            )
        try:
            assert isinstance(kwargs['initial_value'], (Number, MBVar)), (
                    '\n-------------------\nERROR:' +
                    ' SineDriveCaller: <initial_value> must either be a number or an MBVar' + 
                    '\n-------------------\n'
            )
            self.initial_value = kwargs['initial_value']
        except KeyError:
            errprint(
                    '\n-------------------\nWARNING:' +
                    ' SineDriveCaller: <initial_value> not provided, assuming 0.' + 
                    '\n-------------------\n'
            )
            self.initial_value = 0.
            pass
    def __str__(self):
        s = ''
        if self.idx >= 0:
            s = s + 'drive caller: {}, '.format(self.idx)
        s = s + '{}, {}, '.format(self.type, self.initial_time)
        s = s + '{}, {}, '.format(self.angular_velocity, self.amplitude)
        s = s + '{}, {}'.format(self.number_of_cycles, self.initial_value)
        return s

class Data:
    problem_type = ('INITIAL VALUE', 'INVERSE DYNAMICS')
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            if key == 'problem type':
                if value in self.problem_type:
                    self.type = value
                else:
                    raise ValueError('Unrecognised problem type')            

class InitialValueStrategy:
    strategy_type = ('NO CHANGE', 'FACTOR', 'CHANGE')
    def __init__(self, stype, **kwargs):
        if stype in self.strategy_type:
            self.type = value
        else:
            raise ValueError('Unrecognised strategy')
       
        if self.type == 'FACTOR':
            for key, value in kwargs.items():
                if key == 'reduction_factor':
                    self.reduction_factor = value
                if key == 'steps_before_reduction':
                    self.steps_before_reduction = value
                if key == 'raise_factor':
                    self.raise_factor = value
                if key == 'steps_before_raise':
                    self.steps_before_raise = value
                if key == 'minimum_iterations':
                    self.minimum_iterations = value
                if key == 'maximum_iterations':
                    self.maximum_iterations = value
        if self.self_type == 'CHANGE':
            self.time_step_pattern = DriveCaller('const', 1e-3);

class InitialValue:
    def __init__(self):
        self.initial_time = 0.
        self.final_time = 10.
        self.strategy = InitialValueStrategy()
        self.min_time_step = 1e-6
        self.max_time_step = 1.
        self.time_step = 1e-3