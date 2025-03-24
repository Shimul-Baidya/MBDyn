from io import StringIO
import unittest, warnings

import MBDynLib as l

try:
    import pydantic
except ImportError:
    pydantic = None


class ErrprintCalled(Exception):
    """Exception raised instead of `errprint` function outputting to stderr"""
    pass


def patched_errprint(*args, **kwargs):
    output = StringIO()
    print(*args, file=output, **kwargs)
    contents = output.getvalue()
    output.close()
    raise ErrprintCalled(contents)


# put our function instead of the original one (monkeypatching)
l.errprint = patched_errprint


class TestNodeClasses(unittest.TestCase):
    def setUp(self):
        # Create Position2 instances for testing correctly
        self.pos = l.Position2(relative_position=[1.0, 2.0, 3.0], reference='global')
        self.orient = l.Position2(relative_position=[l.eye()], reference='')
        self.vel = l.Position2(relative_position=[0.1, 0.2, 0.3], reference='global')
        self.ang_vel = l.Position2(relative_position=[l.null()], reference='')
    
    def test_node2_initialization(self):
        """Test that Node2 initializes correctly with default values"""
        node = l.Node2(idx=1, position=self.pos, orientation=self.orient, 
                    velocity=self.vel, angular_velocity=self.ang_vel)
        
        self.assertEqual(node.idx, 1)
        self.assertEqual(node.position, self.pos)
        self.assertEqual(node.orientation, self.orient)
        self.assertEqual(node.velocity, self.vel)
        self.assertEqual(node.angular_velocity, self.ang_vel)
        self.assertEqual(node.node_type, 'dynamic')  # Default value
        self.assertEqual(node.scale, 'default')
        self.assertEqual(node.output, 'yes')
    
    def test_dynamic_node2(self):
        """Test DynamicNode2 initialization and string representation"""
        node = l.DynamicNode2(idx=2, pos=self.pos, orient=self.orient, 
                          vel=self.vel, angular_vel=self.ang_vel, 
                          accelerations='yes')
        
        expected_str = (f"structural: 2, dynamic,\n"
                       f"\treference, global, 1.0, 2.0, 3.0,\n"
                       f"\t{self.orient},\n"
                       f"\treference, global, 0.1, 0.2, 0.3,\n"
                       f"\t{self.ang_vel},\n"
                       f"\taccelerations, yes;\n")
        
        self.assertEqual(str(node), expected_str)
    
    def test_static_node2(self):
        """Test StaticNode2 initialization and string representation"""
        node = l.StaticNode2(idx=3, pos=self.pos, orient=self.orient, 
                         vel=self.vel, angular_vel=self.ang_vel)
        
        expected_str = (f"structural: 3, static,\n"
                       f"\treference, global, 1.0, 2.0, 3.0,\n"
                       f"\t{self.orient},\n"
                       f"\treference, global, 0.1, 0.2, 0.3,\n"
                       f"\t{self.ang_vel};\n")
        
        self.assertEqual(str(node), expected_str)
    
    def test_modal_node(self):
        """Test ModalNode initialization and string representation"""
        node = l.ModalNode(idx=4, pos=self.pos, orient=self.orient, 
                       vel=self.vel, angular_vel=self.ang_vel)
        
        expected_str = (f"structural: 4, modal,\n"
                       f"\treference, global, 1.0, 2.0, 3.0,\n"
                       f"\t{self.orient},\n"
                       f"\treference, global, 0.1, 0.2, 0.3,\n"
                       f"\t{self.ang_vel};\n")
        
        self.assertEqual(str(node), expected_str)
    
    def test_displacement_node2(self):
        """Test DisplacementNode2 initialization"""
        node = l.DisplacementNode2(idx=5, position=self.pos, velocity=self.vel)
        
        self.assertEqual(node.idx, 5)
        self.assertEqual(node.position, self.pos)
        self.assertEqual(node.velocity, self.vel)
        self.assertEqual(node.node_type, 'dynamic')  # Default value
    
    def test_dynamic_displacement_node2(self):
        """Test DynamicDisplacementNode2 initialization and string representation"""
        node = l.DynamicDisplacementNode2(idx=6, pos=self.pos, vel=self.vel, 
                                      accelerations='yes')
        
        expected_str = (f"structural: 6, dynamic displacement,\n"
                       f"\treference, global, 1.0, 2.0, 3.0,\n"
                       f"\treference, global, 0.1, 0.2, 0.3,\n"
                       f"\taccelerations, yes;\n")
        
        self.assertEqual(str(node), expected_str)
    
    def test_static_displacement_node2(self):
        """Test StaticDisplacementNode2 initialization and string representation"""
        node = l.StaticDisplacementNode2(idx=7, pos=self.pos, vel=self.vel)
        
        expected_str = (f"structural: 7, static displacement,\n"
                       f"\treference, global, 1.0, 2.0, 3.0,\n"
                       f"\treference, global, 0.1, 0.2, 0.3;\n")
        
        self.assertEqual(str(node), expected_str)

class TestNodeDof(unittest.TestCase):
    def test_node_dof_creation_valid(self):
        """Test creating a NodeDof instance with valid data"""
        # Test with all required fields only
        node_dof = l.NodeDof(
            node_label=1,
            node_type='structural'
        )
        self.assertIsInstance(node_dof, l.NodeDof)
        self.assertEqual(node_dof.node_label, 1)
        self.assertEqual(node_dof.node_type, 'structural')
        self.assertIsNone(node_dof.dof_number)
        self.assertIsNone(node_dof.dof_order)

        # Test with all fields
        node_dof = l.NodeDof(
            node_label=2,
            node_type='electric',
            dof_number=3,
            dof_order='algebraic'
        )
        self.assertIsInstance(node_dof, l.NodeDof)
        self.assertEqual(node_dof.node_label, 2)
        self.assertEqual(node_dof.node_type, 'electric')
        self.assertEqual(node_dof.dof_number, 3)
        self.assertEqual(node_dof.dof_order, 'algebraic')

    def test_node_dof_str_method(self):
        """Test the __str__ method of NodeDof"""
        # Test with required fields only
        node_dof = l.NodeDof(
            node_label=1,
            node_type='structural'
        )
        expected_str = '1, structural'
        self.assertEqual(str(node_dof), expected_str)

        # Test with dof_number
        node_dof = l.NodeDof(
            node_label=2,
            node_type='electric',
            dof_number=3
        )
        expected_str = '2, electric, 3'
        self.assertEqual(str(node_dof), expected_str)

        # Test with dof_order
        node_dof = l.NodeDof(
            node_label=3,
            node_type='parameter',
            dof_order='differential'
        )
        expected_str = '3, parameter, differential'
        self.assertEqual(str(node_dof), expected_str)

        # Test with all fields
        node_dof = l.NodeDof(
            node_label=4,
            node_type='thermal',
            dof_number=5,
            dof_order='algebraic'
        )
        expected_str = '4, thermal, 5, algebraic'
        self.assertEqual(str(node_dof), expected_str)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_node_dof_missing_required_fields(self):
        """Test that NodeDof raises appropriate errors when required fields are missing"""
        # Missing node_label
        with self.assertRaises(Exception):
            l.NodeDof(
                node_type='structural'
            )
        
        # Missing node_type
        with self.assertRaises(Exception):
            l.NodeDof(
                node_label=1
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_node_dof_invalid_types(self):
        """Test that NodeDof raises appropriate errors when fields have invalid types"""
        # Invalid node_type
        with self.assertRaises(Exception):
            l.NodeDof(
                node_label=1,
                node_type='invalid_type'
            )
        
        # Invalid node_label type
        with self.assertRaises(Exception):
            l.NodeDof(
                node_label='invalid',  # Should be int or MBVar
                node_type='structural'
            )
        
        # Invalid dof_number type
        with self.assertRaises(Exception):
            l.NodeDof(
                node_label=1,
                node_type='structural',
                dof_number='invalid'  # Should be int or MBVar
            )
        
        # Invalid dof_order
        with self.assertRaises(Exception):
            l.NodeDof(
                node_label=1,
                node_type='structural',
                dof_order='invalid_order'  # Should be 'algebraic' or 'differential'
            )

    def test_node_dof_with_mbvar(self):
        """Test creating a NodeDof instance with MBVar"""
        if 'node_var' not in l.declared_MBVars:
            node_label_var = l.MBVar(name='node_var', var_type='integer', expression=100)
        else:
            node_label_var = l.declared_MBVars['node_var']

        if 'dof_var' not in l.declared_MBVars:
            dof_number_var = l.MBVar(name='dof_var', var_type='integer', expression=5)
        else:
            dof_number_var = l.declared_MBVars['dof_var']
        
        node_dof = l.NodeDof(
            node_label=node_label_var,
            node_type='structural',
            dof_number=dof_number_var
        )
        
        self.assertEqual(node_dof.node_label, node_label_var)
        self.assertEqual(node_dof.dof_number, dof_number_var)
        self.assertIn(str(node_label_var), str(node_dof))

class TestArrayDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive1 = l.ConstDriveCaller(const_value=42)
        self.const_drive2 = l.ConstDriveCaller(const_value=10)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=20)

    def test_array_drive_caller_creation_valid(self):
        """Test that ArrayDriveCaller works with valid input"""
        # Create with two drives without idx
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive2])
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(len(array_drive.drives), 2)
        self.assertEqual(array_drive.drives[0], self.const_drive1)
        self.assertEqual(array_drive.drives[1], self.const_drive2)
        
        # Create with a mix of drives with and without idx
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive_with_idx])
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(len(array_drive.drives), 2)
        
        # Create with specific idx for the array
        array_drive = l.ArrayDriveCaller(idx=10, drives=[self.const_drive1, self.const_drive2])
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(array_drive.idx, 10)
        
        # Create with single drive (minimum required)
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1])
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(len(array_drive.drives), 1)

    def test_array_drive_caller_str_representation(self):
        """Test the string representation of ArrayDriveCaller"""
        # Test without idx
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive2])
        expected_str = "array, 2,\n\tconst, 42,\n\tconst, 10"
        self.assertEqual(str(array_drive), expected_str)
        
        # Test with idx
        array_drive = l.ArrayDriveCaller(idx=10, drives=[self.const_drive1, self.const_drive2])
        expected_str = "drive caller: 10, array, 2,\n\tconst, 42,\n\tconst, 10"
        self.assertEqual(str(array_drive), expected_str)
        
        # Test with a mix of drives with and without idx
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive_with_idx])
        expected_str = "array, 2,\n\tconst, 42,\n\treference, 5"
        self.assertEqual(str(array_drive), expected_str)
        
        # Test with single drive
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1])
        expected_str = "array, 1,\n\tconst, 42"
        self.assertEqual(str(array_drive), expected_str)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_array_drive_caller_missing_required_field(self):
        """Test creating an ArrayDriveCaller instance missing a required field (drives)"""
        with self.assertRaises(Exception):
            l.ArrayDriveCaller(idx=10)  # Missing drives field

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_array_drive_caller_empty_drives(self):
        """Test creating an ArrayDriveCaller with an empty drives list, which should fail"""
        with self.assertRaises(ValueError) as context:
            l.ArrayDriveCaller(drives=[])
        self.assertIn("array drive must contain at least one drive caller", str(context.exception))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_array_drive_caller_invalid_drives_type(self):
        """Test invalid type for drives field"""
        with self.assertRaises(Exception):
            l.ArrayDriveCaller(drives="not a list")  # Should be a list

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_array_drive_caller_invalid_drive_items(self):
        """Test invalid items in drives list"""
        with self.assertRaises(Exception):
            l.ArrayDriveCaller(drives=[self.const_drive1, "not a drive"])  # One item is not a DriveCaller

    def test_array_drive_caller_nested(self):
        """Test creating a nested ArrayDriveCaller"""
        inner_array = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive2])
        outer_array = l.ArrayDriveCaller(drives=[inner_array, self.const_drive_with_idx])
        self.assertIsInstance(outer_array, l.ArrayDriveCaller)
        self.assertEqual(len(outer_array.drives), 2)
        self.assertIsInstance(outer_array.drives[0], l.ArrayDriveCaller)
        expected_str = "array, 2,\n\tarray, 2,\n\t\tconst, 42,\n\t\tconst, 10,\n\treference, 5"
        self.assertEqual(str(outer_array), expected_str)


class TestBistopDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive1 = l.ConstDriveCaller(const_value=1.0)
        self.const_drive2 = l.ConstDriveCaller(const_value=0.0)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=0.5)

    def test_bistop_drive_caller_creation_valid(self):
        """Test that BistopDriveCaller works with valid inputs"""
        # Default initial_status
        bistop_drive = l.BistopDriveCaller(
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        self.assertIsInstance(bistop_drive, l.BistopDriveCaller)
        self.assertEqual(bistop_drive.initial_status, 'active')
        self.assertEqual(bistop_drive.activation_condition, self.const_drive1)
        self.assertEqual(bistop_drive.deactivation_condition, self.const_drive2)
        
        # With explicit initial_status
        bistop_drive = l.BistopDriveCaller(
            initial_status='inactive',
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        self.assertEqual(bistop_drive.initial_status, 'inactive')
        
        # With idx
        bistop_drive = l.BistopDriveCaller(
            idx=10,
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        self.assertEqual(bistop_drive.idx, 10)

    def test_bistop_drive_caller_str_representation(self):
        """Test the string representation of BistopDriveCaller"""
        # Without idx, with default initial_status
        bistop_drive = l.BistopDriveCaller(
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        expected_str = "bistop,\n\tinitial status, active,\n\t# activation condition drive\n\tconst, 1.0\n\t# deactivation condition drive\n\tconst, 0.0"
        self.assertEqual(str(bistop_drive), expected_str)
        
        # With idx
        bistop_drive = l.BistopDriveCaller(
            idx=10,
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        expected_str = "drive caller: 10, bistop,\n\tinitial status, active,\n\t# activation condition drive\n\tconst, 1.0\n\t# deactivation condition drive\n\tconst, 0.0"
        self.assertEqual(str(bistop_drive), expected_str)
        
        # With reference drive callers
        bistop_drive = l.BistopDriveCaller(
            activation_condition=self.const_drive_with_idx,
            deactivation_condition=self.const_drive_with_idx
        )
        expected_str = "bistop,\n\tinitial status, active,\n\t# activation condition drive\n\treference, 5\n\t# deactivation condition drive\n\treference, 5"
        self.assertEqual(str(bistop_drive), expected_str)
        
        # With custom initial_status
        bistop_drive = l.BistopDriveCaller(
            initial_status='inactive',
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        expected_str = "bistop,\n\tinitial status, inactive,\n\t# activation condition drive\n\tconst, 1.0\n\t# deactivation condition drive\n\tconst, 0.0"
        self.assertEqual(str(bistop_drive), expected_str)

    def test_bistop_drive_caller_drive_type(self):
        """Test the drive_type method of BistopDriveCaller"""
        bistop_drive = l.BistopDriveCaller(
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        self.assertEqual(bistop_drive.drive_type(), 'bistop')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_bistop_drive_caller_missing_required_field(self):
        """Test creating a BistopDriveCaller instance missing a required field"""
        # Missing activation_condition
        with self.assertRaises(Exception):
            l.BistopDriveCaller(
                deactivation_condition=self.const_drive2
            )
        
        # Missing deactivation_condition
        with self.assertRaises(Exception):
            l.BistopDriveCaller(
                activation_condition=self.const_drive1
            )

    unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_bistop_drive_caller_invalid_initial_status(self):
        """Test creating a BistopDriveCaller with invalid initial_status"""
        with self.assertRaises(Exception):
            l.BistopDriveCaller(
                initial_status='invalid_status',  # Invalid value
                activation_condition=self.const_drive1,
                deactivation_condition=self.const_drive2
            )

    unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_bistop_drive_caller_invalid_drive_types(self):
        """Test creating a BistopDriveCaller with invalid drive types"""
        with self.assertRaises(Exception):
            l.BistopDriveCaller(
                activation_condition="not a drive",  # Invalid type
                deactivation_condition=self.const_drive2
            )
        
        with self.assertRaises(Exception):
            l.BistopDriveCaller(
                activation_condition=self.const_drive1,
                deactivation_condition="not a drive"  # Invalid type
            )

    def test_bistop_drive_caller_nested(self):
        """Test creating a BistopDriveCaller with nested BistopDriveCallers"""
        inner_bistop = l.BistopDriveCaller(
            activation_condition=self.const_drive1,
            deactivation_condition=self.const_drive2
        )
        
        outer_bistop = l.BistopDriveCaller(
            activation_condition=inner_bistop,
            deactivation_condition=self.const_drive_with_idx
        )
        
        self.assertIsInstance(outer_bistop, l.BistopDriveCaller)
        self.assertEqual(outer_bistop.activation_condition, inner_bistop)
        self.assertEqual(outer_bistop.deactivation_condition, self.const_drive_with_idx)
        
        # The string representation should properly nest the inner bistop drive
        expected_inner_bistop_str = str(inner_bistop)
        self.assertIn(expected_inner_bistop_str, str(outer_bistop))


class TestConstDrive(unittest.TestCase):
    def test_const_drive_caller(self):
        """Check that the new module can be used the same way as the current one"""

        # just value
        cdc = l.ConstDriveCaller(const_value=42)
        self.assertEqual(cdc.idx, None)
        self.assertEqual(cdc.const_value, 42)
        self.assertEqual(str(cdc), 'const, 42')

        # index and value
        cdc = l.ConstDriveCaller(idx=1, const_value=42)
        self.assertEqual(cdc.idx, 1)
        self.assertEqual(cdc.const_value, 42)
        self.assertEqual(str(cdc), 'drive caller: 1, const, 42')

        # can't use positional arguments, arguably better
        with self.assertRaises(TypeError):
            cdc = l.ConstDriveCaller(42, 1)
            self.assertEqual(cdc.idx, 1)
            self.assertEqual(cdc.const_value, 42)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_missing_arguments(self):
        with self.assertRaises(Exception):
            l.ConstDriveCaller()
        with self.assertRaises(Exception):
            l.ConstDriveCaller(idx=1)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_wrong_argument_types(self):
        with self.assertRaises(Exception):
            l.ConstDriveCaller(idx=1.0)

        with self.assertRaises(Exception):
            l.ConstDriveCaller(const_value='a')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_extra_arguments(self):
        with self.assertRaises(Exception):
            l.ConstDriveCaller(const_value=42, foo=1.0)
        # allows to catch typos in optional arguments
        with self.assertRaises(Exception):
            l.ConstDriveCaller(const_value=42, ibx=1)
        # that also prevents wrongly assigning other members in constructor
        with self.assertRaises(Exception):
            l.ConstDriveCaller(const_value=42, drive_type='bar')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_schema(self):
        pydantic.TypeAdapter(l.ConstDriveCaller).json_schema()

    def test_abstract_class(self):
        """Check that user can't create abstract classes, which are used only to share functionality (can't be part of MBDyn output)"""
        with self.assertRaises(TypeError):
            e = l.MBEntity()
        with self.assertRaises(TypeError):
            dc = l.DriveCaller2()


class TestClosestNextDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive = l.ConstDriveCaller(const_value=5.0)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=10.0)

    def test_closest_next_drive_caller_creation_valid(self):
        """Test that ClosestNextDriveCaller works with valid input"""
        # Create with default initial_time and other required fields
        closest_next_drive = l.ClosestNextDriveCaller(
            final_time='forever',
            increment=self.const_drive
        )
        self.assertIsInstance(closest_next_drive, l.ClosestNextDriveCaller)
        self.assertEqual(closest_next_drive.initial_time, 0.0)
        self.assertEqual(closest_next_drive.final_time, 'forever')
        self.assertEqual(closest_next_drive.increment, self.const_drive)
        
        # Create with custom initial_time and numeric final_time
        closest_next_drive = l.ClosestNextDriveCaller(
            initial_time=1.5,
            final_time=100.0,
            increment=self.const_drive
        )
        self.assertIsInstance(closest_next_drive, l.ClosestNextDriveCaller)
        self.assertEqual(closest_next_drive.initial_time, 1.5)
        self.assertEqual(closest_next_drive.final_time, 100.0)
        
        # Create with specific idx
        closest_next_drive = l.ClosestNextDriveCaller(
            idx=10,
            final_time='forever',
            increment=self.const_drive
        )
        self.assertIsInstance(closest_next_drive, l.ClosestNextDriveCaller)
        self.assertEqual(closest_next_drive.idx, 10)

    def test_closest_next_drive_caller_str_representation(self):
        """Test the string representation of ClosestNextDriveCaller"""
        # Test without idx
        closest_next_drive = l.ClosestNextDriveCaller(
            final_time='forever',
            increment=self.const_drive
        )
        expected_str = "closest next,\n\t0.0, forever,\n\t# increment drive\n\tconst, 5.0"
        self.assertEqual(str(closest_next_drive), expected_str)
        
        # Test with idx
        closest_next_drive = l.ClosestNextDriveCaller(
            idx=10,
            initial_time=2.0,
            final_time=50.0,
            increment=self.const_drive
        )
        expected_str = "drive caller: 10, closest next,\n\t2.0, 50.0,\n\t# increment drive\n\tconst, 5.0"
        self.assertEqual(str(closest_next_drive), expected_str)
        
        # Test with reference drive
        closest_next_drive = l.ClosestNextDriveCaller(
            final_time='forever',
            increment=self.const_drive_with_idx
        )
        expected_str = "closest next,\n\t0.0, forever,\n\t# increment drive\n\treference, 5"
        self.assertEqual(str(closest_next_drive), expected_str)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_closest_next_drive_caller_missing_required_field(self):
        """Test creating a ClosestNextDriveCaller instance missing a required field"""
        # Missing final_time
        with self.assertRaises(Exception):
            l.ClosestNextDriveCaller(
                increment=self.const_drive
            )
        
        # Missing increment
        with self.assertRaises(Exception):
            l.ClosestNextDriveCaller(
                final_time='forever'
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_closest_next_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for initial_time
        with self.assertRaises(Exception):
            l.ClosestNextDriveCaller(
                initial_time="invalid",
                final_time='forever',
                increment=self.const_drive
            )
        
        # Invalid type for final_time
        with self.assertRaises(Exception):
            l.ClosestNextDriveCaller(
                final_time='time',  # Should be float, MBVar, or 'forever'
                increment=self.const_drive
            )
        
        # Invalid type for increment
        with self.assertRaises(Exception):
            l.ClosestNextDriveCaller(
                final_time='forever',
                increment="not a drive"  # Should be a DriveCaller
            )

    def test_closest_next_drive_caller_nested(self):
        """Test creating nested ClosestNextDriveCaller"""
        # Create an ArrayDriveCaller as the increment
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive, self.const_drive_with_idx])
        
        closest_next_drive = l.ClosestNextDriveCaller(
            final_time='forever',
            increment=array_drive
        )
        
        self.assertIsInstance(closest_next_drive, l.ClosestNextDriveCaller)
        self.assertEqual(closest_next_drive.increment, array_drive)
        
        # Verify the string representation with nested drive
        expected_str = "closest next,\n\t0.0, forever,\n\t# increment drive\n\tarray, 2,\n\tconst, 5.0,\n\treference, 5"
        self.assertEqual(str(closest_next_drive), expected_str)


class TestCosineDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.angular_velocity = 2.0
        self.amplitude = 1.5
        self.number_of_cycles = 'forever'
        self.initial_value = 0.5
        
        if 'test_angular_velocity' not in l.declared_MBVars:
            self.angular_velocity_var = l.MBVar(name='test_angular_velocity', var_type='real', expression=3.14)
        else:
            self.angular_velocity_var = l.declared_MBVars['test_angular_velocity']
            
        if 'test_amplitude' not in l.declared_MBVars:
            self.amplitude_var = l.MBVar(name='test_amplitude', var_type='real', expression=2.5)
        else:
            self.amplitude_var = l.declared_MBVars['test_amplitude']
            
        if 'test_initial_time' not in l.declared_MBVars:
            self.initial_time_var = l.MBVar(name='test_initial_time', var_type='real', expression=1.0)
        else:
            self.initial_time_var = l.declared_MBVars['test_initial_time']

    def test_cosine_drive_caller_creation_valid(self):
        """Test that CosineDriveCaller works with valid input"""
        # Create with minimum required parameters (defaults for initial_time and initial_value)
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        self.assertIsInstance(cosine_drive, l.CosineDriveCaller)
        self.assertEqual(cosine_drive.initial_time, 0.0)
        self.assertEqual(cosine_drive.angular_velocity, self.angular_velocity)
        self.assertEqual(cosine_drive.amplitude, self.amplitude)
        self.assertEqual(cosine_drive.number_of_cycles, self.number_of_cycles)
        self.assertEqual(cosine_drive.initial_value, 0.0)
        
        # Create with all parameters specified
        cosine_drive = l.CosineDriveCaller(
            initial_time=1.0,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(cosine_drive, l.CosineDriveCaller)
        self.assertEqual(cosine_drive.initial_time, 1.0)
        self.assertEqual(cosine_drive.initial_value, self.initial_value)
        
        # Create with specific idx
        cosine_drive = l.CosineDriveCaller(
            idx=10,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        self.assertIsInstance(cosine_drive, l.CosineDriveCaller)
        self.assertEqual(cosine_drive.idx, 10)
        
        # Test with string literals for number_of_cycles
        for cycles in ['half', 'one', 'forever']:
            cosine_drive = l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles=cycles
            )
            self.assertEqual(cosine_drive.number_of_cycles, cycles)
            
        # Test with numeric value for number_of_cycles
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=2.5
        )
        self.assertEqual(cosine_drive.number_of_cycles, 2.5)

    def test_cosine_drive_caller_with_mbvars(self):
        """Test CosineDriveCaller with MBVar objects for parameters"""
        # Create with MBVar for angular_velocity
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity_var,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        self.assertIsInstance(cosine_drive, l.CosineDriveCaller)
        self.assertEqual(cosine_drive.angular_velocity, self.angular_velocity_var)
        
        # Create with MBVar for amplitude
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude_var,
            number_of_cycles=self.number_of_cycles
        )
        self.assertEqual(cosine_drive.amplitude, self.amplitude_var)
        
        # Create with MBVar for initial_time
        cosine_drive = l.CosineDriveCaller(
            initial_time=self.initial_time_var,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        self.assertEqual(cosine_drive.initial_time, self.initial_time_var)
        
        # Create with multiple MBVar parameters
        cosine_drive = l.CosineDriveCaller(
            initial_time=self.initial_time_var,
            angular_velocity=self.angular_velocity_var,
            amplitude=self.amplitude_var,
            number_of_cycles=self.number_of_cycles
        )
        self.assertEqual(cosine_drive.initial_time, self.initial_time_var)
        self.assertEqual(cosine_drive.angular_velocity, self.angular_velocity_var)
        self.assertEqual(cosine_drive.amplitude, self.amplitude_var)

    def test_cosine_drive_caller_str_representation(self):
        """Test the string representation of CosineDriveCaller"""
        # Test without idx, with default values for initial_time and initial_value
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        expected_str = "cosine,\n\t0.0, 2.0, 1.5, forever, 0.0"
        self.assertEqual(str(cosine_drive), expected_str)
        
        # Test with idx and custom values
        cosine_drive = l.CosineDriveCaller(
            idx=10,
            initial_time=1.0,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles='one',
            initial_value=0.5
        )
        expected_str = "drive caller: 10, cosine,\n\t1.0, 2.0, 1.5, one, 0.5"
        self.assertEqual(str(cosine_drive), expected_str)
        
        # Test with MBVar parameters
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity_var,
            amplitude=self.amplitude_var,
            number_of_cycles=self.number_of_cycles
        )
        expected_str = f"cosine,\n\t0.0, {self.angular_velocity_var}, {self.amplitude_var}, forever, 0.0"
        self.assertEqual(str(cosine_drive), expected_str)

    def test_cosine_drive_caller_drive_type(self):
        """Test the drive_type method of CosineDriveCaller"""
        cosine_drive = l.CosineDriveCaller(
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles
        )
        self.assertEqual(cosine_drive.drive_type(), 'cosine')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_cosine_drive_caller_missing_required_field(self):
        """Test creating a CosineDriveCaller instance missing a required field"""
        # Missing angular_velocity
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        
        # Missing amplitude
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                number_of_cycles=self.number_of_cycles
            )
        
        # Missing number_of_cycles
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_cosine_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for angular_velocity
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity="invalid",
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        
        # Invalid type for amplitude
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude="invalid",
                number_of_cycles=self.number_of_cycles
            )
        
        # Invalid type for number_of_cycles
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles='Invalid'  # Should be float, MBVar, or specific strings
            )
        
        # Invalid value for number_of_cycles
        with self.assertRaises(Exception):
            l.CosineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles="invalid_value"  # Should be 'half', 'one', or 'forever' if a string
            )

class TestCubicDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create variables for testing
        if 'test_const_coef' not in l.declared_MBVars:
            self.const_coef_var = l.MBVar(name='test_const_coef', var_type='real', expression=1.0)
        else:
            self.const_coef_var = l.declared_MBVars['test_const_coef']
            
    def test_cubic_drive_caller_creation_valid(self):
        """Test that CubicDriveCaller works with valid input"""
        # Create with required parameters
        cubic_drive = l.CubicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        self.assertIsInstance(cubic_drive, l.CubicDriveCaller)
        self.assertEqual(cubic_drive.const_coef, 1.0)
        self.assertEqual(cubic_drive.linear_coef, 2.0)
        self.assertEqual(cubic_drive.parabolic_coef, 3.0)
        self.assertEqual(cubic_drive.cubic_coef, 4.0)
        
        # Create with specific idx
        cubic_drive = l.CubicDriveCaller(
            idx=10,
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        self.assertIsInstance(cubic_drive, l.CubicDriveCaller)
        self.assertEqual(cubic_drive.idx, 10)

    def test_cubic_drive_caller_with_mbvars(self):
        """Test CubicDriveCaller with MBVar objects for parameters"""
        # Create with MBVar for const_coef
        cubic_drive = l.CubicDriveCaller(
            const_coef=self.const_coef_var,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        self.assertIsInstance(cubic_drive, l.CubicDriveCaller)
        self.assertEqual(cubic_drive.const_coef, self.const_coef_var)

    def test_cubic_drive_caller_str_representation(self):
        """Test the string representation of CubicDriveCaller"""
        # Test without idx
        cubic_drive = l.CubicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        expected_str = "cubic, 1.0, 2.0, 3.0, 4.0"
        self.assertEqual(str(cubic_drive), expected_str)
        
        # Test with idx
        cubic_drive = l.CubicDriveCaller(
            idx=10,
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        expected_str = "drive caller: 10, cubic, 1.0, 2.0, 3.0, 4.0"
        self.assertEqual(str(cubic_drive), expected_str)
        
        # Test with MBVar
        cubic_drive = l.CubicDriveCaller(
            const_coef=self.const_coef_var,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        expected_str = f"cubic, {self.const_coef_var}, 2.0, 3.0, 4.0"
        self.assertEqual(str(cubic_drive), expected_str)

    def test_cubic_drive_caller_drive_type(self):
        """Test the drive_type method of CubicDriveCaller"""
        cubic_drive = l.CubicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0,
            cubic_coef=4.0
        )
        self.assertEqual(cubic_drive.drive_type(), 'cubic')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_cubic_drive_caller_missing_required_field(self):
        """Test creating a CubicDriveCaller instance missing a required field"""
        # Missing const_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                linear_coef=2.0,
                parabolic_coef=3.0,
                cubic_coef=4.0
            )
        
        # Missing linear_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                const_coef=1.0,
                parabolic_coef=3.0,
                cubic_coef=4.0
            )
        
        # Missing parabolic_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                const_coef=1.0,
                linear_coef=2.0,
                cubic_coef=4.0
            )
        
        # Missing cubic_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                const_coef=1.0,
                linear_coef=2.0,
                parabolic_coef=3.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_cubic_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for const_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                const_coef="invalid",
                linear_coef=2.0,
                parabolic_coef=3.0,
                cubic_coef=4.0
            )
        
        # Invalid type for linear_coef
        with self.assertRaises(Exception):
            l.CubicDriveCaller(
                const_coef=1.0,
                linear_coef="invalid",
                parabolic_coef=3.0,
                cubic_coef=4.0
            )


class TestDirectDriveCaller(unittest.TestCase):
    def test_direct_drive_caller_creation_valid(self):
        """Test that DirectDriveCaller works with valid input"""
        # Create without idx
        direct_drive = l.DirectDriveCaller()
        self.assertIsInstance(direct_drive, l.DirectDriveCaller)
        self.assertEqual(direct_drive.idx, None)
        
        # Create with specific idx
        direct_drive = l.DirectDriveCaller(idx=10)
        self.assertIsInstance(direct_drive, l.DirectDriveCaller)
        self.assertEqual(direct_drive.idx, 10)

    def test_direct_drive_caller_str_representation(self):
        """Test the string representation of DirectDriveCaller"""
        # Test without idx
        direct_drive = l.DirectDriveCaller()
        expected_str = "direct"
        self.assertEqual(str(direct_drive), expected_str)
        
        # Test with idx
        direct_drive = l.DirectDriveCaller(idx=10)
        expected_str = "drive caller: 10, direct"
        self.assertEqual(str(direct_drive), expected_str)

    def test_direct_drive_caller_drive_type(self):
        """Test the drive_type method of DirectDriveCaller"""
        direct_drive = l.DirectDriveCaller()
        self.assertEqual(direct_drive.drive_type(), 'direct')

class TestDiscreteFilterDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive caller for testing
        self.const_drive = l.ConstDriveCaller(const_value=1.5)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=2.0)
        
        # Create MBVar objects for testing
        if 'test_n_a' not in l.declared_MBVars:
            self.n_a_var = l.MBVar(name='test_n_a', var_type='integer', expression=2)
        else:
            self.n_a_var = l.declared_MBVars['test_n_a']
            
        if 'test_b_0' not in l.declared_MBVars:
            self.b_0_var = l.MBVar(name='test_b_0', var_type='real', expression=0.5)
        else:
            self.b_0_var = l.declared_MBVars['test_b_0']

    def test_discrete_filter_drive_caller_creation_valid(self):
        """Test that DiscreteFilterDriveCaller works with valid input"""
        # Create with integer values for n_a and n_b
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        self.assertIsInstance(discrete_filter, l.DiscreteFilterDriveCaller)
        self.assertEqual(discrete_filter.n_a, 2)
        self.assertEqual(discrete_filter.a, [0.1, 0.2])
        self.assertEqual(discrete_filter.b_0, 0.5)
        self.assertEqual(discrete_filter.n_b, 3)
        self.assertEqual(discrete_filter.b, [0.3, 0.4, 0.5])
        self.assertEqual(discrete_filter.input_drive, self.const_drive)
        
        # Create with MBVar for n_a
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=self.n_a_var,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        self.assertEqual(discrete_filter.n_a, self.n_a_var)
        
        # Create with MBVar for b_0
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=self.b_0_var,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        self.assertEqual(discrete_filter.b_0, self.b_0_var)
        
        # Create with specific idx
        discrete_filter = l.DiscreteFilterDriveCaller(
            idx=10,
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        self.assertEqual(discrete_filter.idx, 10)
        
        # Create with reference drive
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive_with_idx
        )
        self.assertEqual(discrete_filter.input_drive, self.const_drive_with_idx)

    def test_discrete_filter_drive_caller_str_representation(self):
        """Test the string representation of DiscreteFilterDriveCaller"""
        # Test without idx
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        expected_str = "discrete filter,\n\t2, 0.1, 0.2,\n\t0.5,\n\t3, 0.3, 0.4, 0.5,\n\tconst, 1.5"
        self.assertEqual(str(discrete_filter), expected_str)
        
        # Test with idx
        discrete_filter = l.DiscreteFilterDriveCaller(
            idx=10,
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        expected_str = "drive caller: 10, discrete filter,\n\t2, 0.1, 0.2,\n\t0.5,\n\t3, 0.3, 0.4, 0.5,\n\tconst, 1.5"
        self.assertEqual(str(discrete_filter), expected_str)
        
        # Test with reference drive
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive_with_idx
        )
        expected_str = "discrete filter,\n\t2, 0.1, 0.2,\n\t0.5,\n\t3, 0.3, 0.4, 0.5,\n\treference, 5"
        self.assertEqual(str(discrete_filter), expected_str)
        
        # Test with MBVar parameters
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=self.n_a_var,
            a=[0.1, 0.2],
            b_0=self.b_0_var,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        expected_str = f"discrete filter,\n\t{self.n_a_var}, 0.1, 0.2,\n\t{self.b_0_var},\n\t3, 0.3, 0.4, 0.5,\n\tconst, 1.5"
        self.assertEqual(str(discrete_filter), expected_str)

    def test_discrete_filter_drive_caller_drive_type(self):
        """Test the drive_type method of DiscreteFilterDriveCaller"""
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=self.const_drive
        )
        self.assertEqual(discrete_filter.drive_type(), 'discrete filter')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_discrete_filter_drive_caller_missing_required_field(self):
        """Test creating a DiscreteFilterDriveCaller instance missing a required field"""
        # Missing n_a
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Missing a
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Missing b_0
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Missing n_b
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0=0.5,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Missing b
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                input_drive=self.const_drive
            )
        
        # Missing input_drive
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5]
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_discrete_filter_drive_caller_coefficient_length_validation(self):
        """Test validation of coefficient array lengths"""
        # a list length doesn't match n_a
        with self.assertRaises(ValueError) as context:
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1],  # Should have 2 elements
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        self.assertIn("Length of 'a' list", str(context.exception))
        
        # b list length doesn't match n_b
        with self.assertRaises(ValueError) as context:
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4],  # Should have 3 elements
                input_drive=self.const_drive
            )
        self.assertIn("Length of 'b' list", str(context.exception))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_discrete_filter_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for n_a
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a="invalid",
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Invalid type for a
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a="invalid",  # Should be a list
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Invalid type for b_0
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0="invalid",
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive=self.const_drive
            )
        
        # Invalid type for input_drive
        with self.assertRaises(Exception):
            l.DiscreteFilterDriveCaller(
                n_a=2,
                a=[0.1, 0.2],
                b_0=0.5,
                n_b=3,
                b=[0.3, 0.4, 0.5],
                input_drive="not a drive"  # Should be a DriveCaller
            )

    def test_discrete_filter_drive_caller_nested(self):
        """Test nesting DiscreteFilterDriveCaller with other drive callers"""
        # Create an ArrayDriveCaller to use as input
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive, self.const_drive_with_idx])
        
        # Use it as input to a DiscreteFilterDriveCaller
        discrete_filter = l.DiscreteFilterDriveCaller(
            n_a=2,
            a=[0.1, 0.2],
            b_0=0.5,
            n_b=3,
            b=[0.3, 0.4, 0.5],
            input_drive=array_drive
        )
        
        self.assertIsInstance(discrete_filter, l.DiscreteFilterDriveCaller)
        self.assertEqual(discrete_filter.input_drive, array_drive)
        
        # Check string representation with nested drive
        expected_str = "discrete filter,\n\t2, 0.1, 0.2,\n\t0.5,\n\t3, 0.3, 0.4, 0.5,\n\tarray, 2,\n\tconst, 1.5,\n\treference, 5"
        self.assertEqual(str(discrete_filter), expected_str)

class TestDofDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive = l.ConstDriveCaller(const_value=1.5)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=2.0)
        
        # Create NodeDof instances for testing
        self.node_dof = l.NodeDof(node_label=1, node_type="structural", dof_number=2)
        self.node_dof_with_order = l.NodeDof(
            node_label=2, 
            node_type="structural", 
            dof_number=3, 
            dof_order="differential"
        )

    def test_dof_drive_caller_creation_valid(self):
        """Test that DofDriveCaller works with valid input"""
        # Create with required parameters
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        self.assertIsInstance(dof_drive, l.DofDriveCaller)
        self.assertEqual(dof_drive.driving_dof, self.node_dof)
        self.assertEqual(dof_drive.func_drive, self.const_drive)
        
        # Create with specific idx
        dof_drive = l.DofDriveCaller(
            idx=10,
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        self.assertIsInstance(dof_drive, l.DofDriveCaller)
        self.assertEqual(dof_drive.idx, 10)
        
        # Create with reference drive
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive_with_idx
        )
        self.assertEqual(dof_drive.func_drive, self.const_drive_with_idx)
        
        # Create with a different NodeDof that includes dof_order
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof_with_order,
            func_drive=self.const_drive
        )
        self.assertEqual(dof_drive.driving_dof, self.node_dof_with_order)

    def test_dof_drive_caller_str_representation(self):
        """Test the string representation of DofDriveCaller"""
        # Test without idx
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        expected_str = "dof,\n\t1, structural, 2,\n\tconst, 1.5"
        self.assertEqual(str(dof_drive), expected_str)
        
        # Test with idx
        dof_drive = l.DofDriveCaller(
            idx=10,
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        expected_str = "drive caller: 10, dof,\n\t1, structural, 2,\n\tconst, 1.5"
        self.assertEqual(str(dof_drive), expected_str)
        
        # Test with reference drive
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive_with_idx
        )
        expected_str = "dof,\n\t1, structural, 2,\n\treference, 5"
        self.assertEqual(str(dof_drive), expected_str)
        
        # Test with NodeDof that includes dof_order
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof_with_order,
            func_drive=self.const_drive
        )
        expected_str = "dof,\n\t2, structural, 3, differential,\n\tconst, 1.5"
        self.assertEqual(str(dof_drive), expected_str)

    def test_dof_drive_caller_drive_type(self):
        """Test the drive_type method of DofDriveCaller"""
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        self.assertEqual(dof_drive.drive_type(), 'dof')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_dof_drive_caller_missing_required_field(self):
        """Test creating a DofDriveCaller instance missing a required field"""
        # Missing driving_dof
        with self.assertRaises(Exception):
            l.DofDriveCaller(
                func_drive=self.const_drive
            )
        
        # Missing func_drive
        with self.assertRaises(Exception):
            l.DofDriveCaller(
                driving_dof=self.node_dof
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_dof_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for driving_dof
        with self.assertRaises(Exception):
            l.DofDriveCaller(
                driving_dof="not a NodeDof",
                func_drive=self.const_drive
            )
        
        # Invalid type for func_drive
        with self.assertRaises(Exception):
            l.DofDriveCaller(
                driving_dof=self.node_dof,
                func_drive="not a drive"
            )

    def test_dof_drive_caller_nested(self):
        """Test nesting DofDriveCaller with other drive callers"""
        # Create an ArrayDriveCaller to use as func_drive
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive, self.const_drive_with_idx])
        
        # Use it as func_drive in a DofDriveCaller
        dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=array_drive
        )
        
        self.assertIsInstance(dof_drive, l.DofDriveCaller)
        self.assertEqual(dof_drive.func_drive, array_drive)
        
        # Check string representation with nested drive
        expected_str = "dof,\n\t1, structural, 2,\n\tarray, 2,\n\tconst, 1.5,\n\treference, 5"
        self.assertEqual(str(dof_drive), expected_str)

    def test_dof_drive_caller_complex_nesting(self):
        """Test complex nesting with DofDriveCaller"""
        # Create a DofDriveCaller to be used as func_drive in another DofDriveCaller
        inner_dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof,
            func_drive=self.const_drive
        )
        
        # Create a DofDriveCaller that uses another DofDriveCaller as its func_drive
        outer_dof_drive = l.DofDriveCaller(
            driving_dof=self.node_dof_with_order,
            func_drive=inner_dof_drive
        )
        
        self.assertIsInstance(outer_dof_drive, l.DofDriveCaller)
        self.assertEqual(outer_dof_drive.driving_dof, self.node_dof_with_order)
        self.assertEqual(outer_dof_drive.func_drive, inner_dof_drive)
        
        # Check the string representation with nested dof drive
        expected_str = "dof,\n\t2, structural, 3, differential,\n\tdof,\n\t1, structural, 2,\n\tconst, 1.5"
        self.assertEqual(str(outer_dof_drive), expected_str)

class TestDoubleRampDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create MBVar objects for testing
        if 'test_slope' not in l.declared_MBVars:
            self.slope_var = l.MBVar(name='test_slope', var_type='real', expression=2.5)
        else:
            self.slope_var = l.declared_MBVars['test_slope']
            
        if 'test_time' not in l.declared_MBVars:
            self.time_var = l.MBVar(name='test_time', var_type='real', expression=3.0)
        else:
            self.time_var = l.declared_MBVars['test_time']

    def test_double_ramp_drive_caller_creation_valid(self):
        """Test that DoubleRampDriveCaller works with valid input"""
        # Create with all required parameters (a_initial_time is optional with default)
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertIsInstance(double_ramp_drive, l.DoubleRampDriveCaller)
        self.assertEqual(double_ramp_drive.a_slope, 1.0)
        self.assertEqual(double_ramp_drive.a_initial_time, 0.0)  # Default value
        self.assertEqual(double_ramp_drive.a_final_time, 2.0)
        self.assertEqual(double_ramp_drive.d_slope, 0.5)
        self.assertEqual(double_ramp_drive.d_initial_time, 3.0)
        self.assertEqual(double_ramp_drive.d_final_time, 4.0)
        self.assertEqual(double_ramp_drive.initial_value, 0.0)
        
        # Create with explicit a_initial_time
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_initial_time=0.5,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertEqual(double_ramp_drive.a_initial_time, 0.5)
        
        # Create with 'forever' for d_final_time
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time='forever',
            initial_value=0.0
        )
        self.assertEqual(double_ramp_drive.d_final_time, 'forever')
        
        # Create with specific idx
        double_ramp_drive = l.DoubleRampDriveCaller(
            idx=10,
            a_slope=1.0,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertEqual(double_ramp_drive.idx, 10)

    def test_double_ramp_drive_caller_with_mbvars(self):
        """Test DoubleRampDriveCaller with MBVar objects for parameters"""
        # Create with MBVar for slopes
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=self.slope_var,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertIsInstance(double_ramp_drive, l.DoubleRampDriveCaller)
        self.assertEqual(double_ramp_drive.a_slope, self.slope_var)
        
        # Create with MBVar for times
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_final_time=self.time_var,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertEqual(double_ramp_drive.a_final_time, self.time_var)

    def test_double_ramp_drive_caller_default_a_initial_time(self):
        """Test the default behavior of a_initial_time if not provided"""
        # Test that a warning is issued when a_initial_time is not provided
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            double_ramp_drive = l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
            
            # Verify a warning was raised
            self.assertTrue(len(w) > 0)
            self.assertTrue(any("<a_initial_time> is not set" in str(warning.message) for warning in w))
            
            # Verify the default value was set
            self.assertEqual(double_ramp_drive.a_initial_time, 0.0)

    def test_double_ramp_drive_caller_default_warning(self):
        """Test that a warning is issued when a_initial_time isn't provided"""
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered
            warnings.simplefilter("always")
            
            # Create a DoubleRampDriveCaller without specifying a_initial_time
            double_ramp_drive = l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=5.0,
                d_slope=2.0,
                d_initial_time=10.0,
                d_final_time=15.0,
                initial_value=0.0
            )
            
            # Verify a warning was raised
            self.assertEqual(len(w), 1)
            self.assertIn("<a_initial_time> is not set, assuming 0.0.", str(w[0].message))

        # No warning when a_initial_time is explicitly provided
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            double_ramp_drive = l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_initial_time=1.0,  # Explicitly provided
                a_final_time=5.0,
                d_slope=2.0,
                d_initial_time=10.0,
                d_final_time=15.0,
                initial_value=0.0
            )
            
            self.assertEqual(len(w), 0)  # No warnings

    def test_double_ramp_drive_caller_str_representation(self):
        """Test the string representation of DoubleRampDriveCaller"""
        # Test without idx
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_initial_time=0.5,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        expected_str = "double ramp,\n\t1.0, 0.5, 2.0,\n\t0.5, 3.0, 4.0,\n\t0.0"
        self.assertEqual(str(double_ramp_drive), expected_str)
        
        # Test with idx
        double_ramp_drive = l.DoubleRampDriveCaller(
            idx=10,
            a_slope=1.0,
            a_initial_time=0.5,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        expected_str = "drive caller: 10, double ramp,\n\t1.0, 0.5, 2.0,\n\t0.5, 3.0, 4.0,\n\t0.0"
        self.assertEqual(str(double_ramp_drive), expected_str)
        
        # Test with 'forever'
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_initial_time=0.5,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time='forever',
            initial_value=0.0
        )
        expected_str = "double ramp,\n\t1.0, 0.5, 2.0,\n\t0.5, 3.0, forever,\n\t0.0"
        self.assertEqual(str(double_ramp_drive), expected_str)
        
        # Test with MBVar parameters
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=self.slope_var,
            a_final_time=self.time_var,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        expected_str = f"double ramp,\n\t{self.slope_var}, 0.0, {self.time_var},\n\t0.5, 3.0, 4.0,\n\t0.0"
        self.assertEqual(str(double_ramp_drive), expected_str)

    def test_double_ramp_drive_caller_drive_type(self):
        """Test the drive_type method of DoubleRampDriveCaller"""
        double_ramp_drive = l.DoubleRampDriveCaller(
            a_slope=1.0,
            a_final_time=2.0,
            d_slope=0.5,
            d_initial_time=3.0,
            d_final_time=4.0,
            initial_value=0.0
        )
        self.assertEqual(double_ramp_drive.drive_type(), 'double ramp')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_double_ramp_drive_caller_missing_required_field(self):
        """Test creating a DoubleRampDriveCaller instance missing a required field"""
        # Missing a_slope
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Missing a_final_time
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Missing d_slope
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Missing d_initial_time
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_slope=0.5,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Missing d_final_time
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                initial_value=0.0
            )
        
        # Missing initial_value
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_double_ramp_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for a_slope
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope="invalid",
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Invalid type for a_initial_time
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_initial_time="invalid",
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time=4.0,
                initial_value=0.0
            )
        
        # Invalid type for d_final_time (should be float, MBVar, or 'forever')
        with self.assertRaises(Exception):
            l.DoubleRampDriveCaller(
                a_slope=1.0,
                a_final_time=2.0,
                d_slope=0.5,
                d_initial_time=3.0,
                d_final_time="invalid",
                initial_value=0.0
            )

class TestDoubleStepDriveCaller(unittest.TestCase):
    def setUp(self):
        # Reset warnings to make sure we capture them in tests
        warnings.resetwarnings()
        # Setup common values for testing
        self.initial_time = 1.0
        self.final_time = 5.0
        self.step_value = 10.0
        self.initial_value = 2.0
    
    def test_double_step_drive_caller_creation_valid(self):
        """Test creating a DoubleStepDriveCaller with valid parameters"""
        drive = l.DoubleStepDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time,
            step_value=self.step_value,
            initial_value=self.initial_value
        )
        
        # Verify all properties are set correctly
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.final_time, self.final_time)
        self.assertEqual(drive.step_value, self.step_value)
        self.assertEqual(drive.initial_value, self.initial_value)
        self.assertEqual(drive.drive_type(), "double step")
    
    def test_double_step_drive_caller_default_values(self):
        """Test that default values are set correctly with warnings"""
        with warnings.catch_warnings(record=True) as w:
            # Create drive without initial_time and initial_value
            drive = l.DoubleStepDriveCaller(
                final_time=self.final_time,
                step_value=self.step_value
            )
            
            # Check default values
            self.assertEqual(drive.initial_time, 0.0)
            self.assertEqual(drive.initial_value, 0.0)
            
            # Verify warnings were raised
            self.assertEqual(len(w), 2)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertTrue("<initial_time> is not set, assuming 0.0." in str(w[0].message))
            self.assertTrue(issubclass(w[1].category, UserWarning))
            self.assertTrue("<initial_value> is not set, assuming 0.0." in str(w[1].message))
    
    def test_double_step_drive_caller_str_representation(self):
        """Test string representation of the drive caller"""
        # Create a drive caller with an index
        drive = l.DoubleStepDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            final_time=self.final_time,
            step_value=self.step_value,
            initial_value=self.initial_value
        )
        
        # Expected string representation
        expected_str = "drive caller: 5, double step,\n\t1.0, 5.0,\n\t10.0, 2.0"
        
        # Test string representation
        self.assertEqual(str(drive), expected_str)
        
        # Create a drive caller without an index
        drive_no_idx = l.DoubleStepDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time,
            step_value=self.step_value,
            initial_value=self.initial_value
        )
        
        # Expected string representation
        expected_str_no_idx = "double step,\n\t1.0, 5.0,\n\t10.0, 2.0"
        
        # Test string representation
        self.assertEqual(str(drive_no_idx), expected_str_no_idx)
    
    def test_double_step_drive_caller_with_mbvars(self):
        """Test creating a DoubleStepDriveCaller with MBVar objects"""
        # Create MBVar objects
        if 'init_time' not in l.declared_MBVars:
            initial_time_var = l.MBVar(name='init_time', var_type='real', expression=1.5)
        else:
            initial_time_var = l.declared_MBVars['init_time']
        if 'final_time' not in l.declared_MBVars:
            final_time_var = l.MBVar(name='final_time', var_type='real', expression=6.0)
        else:
            final_time_var = l.declared_MBVars['final_time']
        if 'step_val' not in l.declared_MBVars:
            step_value_var = l.MBVar(name='step_val', var_type='real', expression=12.5)
        else:
            step_value_var = l.declared_MBVars['step_val']
        if 'init_val' not in l.declared_MBVars:
            initial_value_var = l.MBVar(name='init_val', var_type='real', expression=3.5)
        else:
            initial_value_var = l.declared_MBVars['init_val']
        
        # Create drive with MBVar objects
        drive = l.DoubleStepDriveCaller(
            initial_time=initial_time_var,
            final_time=final_time_var,
            step_value=step_value_var,
            initial_value=initial_value_var
        )
        
        # Check that MBVar references are stored correctly
        self.assertEqual(drive.initial_time, initial_time_var)
        self.assertEqual(drive.final_time, final_time_var)
        self.assertEqual(drive.step_value, step_value_var)
        self.assertEqual(drive.initial_value, initial_value_var)
        
        # Check string representation with MBVars
        expected_str = "double step,\n\tinit_time, final_time,\n\tstep_val, init_val"
        self.assertEqual(str(drive), expected_str)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_double_step_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing final_time
        with self.assertRaises(ValueError):
            l.DoubleStepDriveCaller(
                step_value=self.step_value
            )
        
        # Missing step_value
        with self.assertRaises(ValueError):
            l.DoubleStepDriveCaller(
                final_time=self.final_time
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_double_step_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # final_time as string
        with self.assertRaises(ValueError):
            l.DoubleStepDriveCaller(
                final_time="invalid",
                step_value=self.step_value
            )
        
        # step_value as string
        with self.assertRaises(ValueError):
            l.DoubleStepDriveCaller(
                final_time=self.final_time,
                step_value="invalid"
            )
        
        # initial_value as list
        with self.assertRaises(ValueError):
            l.DoubleStepDriveCaller(
                final_time=self.final_time,
                step_value=self.step_value,
                initial_value=[1, 2, 3]
            )

class TestDriveDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive1 = l.ConstDriveCaller(const_value=1.5)
        self.const_drive2 = l.ConstDriveCaller(const_value=2.5)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=3.0)
        self.direct_drive = l.DirectDriveCaller()
        
    def test_drive_drive_caller_creation_valid(self):
        """Test that DriveDriveCaller works with valid input"""
        # Create with required parameters
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        self.assertIsInstance(drive_drive, l.DriveDriveCaller)
        self.assertEqual(drive_drive.drive_caller1, self.const_drive1)
        self.assertEqual(drive_drive.drive_caller2, self.const_drive2)
        
        # Create with specific idx
        drive_drive = l.DriveDriveCaller(
            idx=10,
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        self.assertIsInstance(drive_drive, l.DriveDriveCaller)
        self.assertEqual(drive_drive.idx, 10)
        
        # Create with a mix of drive types
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.direct_drive
        )
        self.assertEqual(drive_drive.drive_caller1, self.const_drive1)
        self.assertEqual(drive_drive.drive_caller2, self.direct_drive)

    def test_drive_drive_caller_with_reference_drives(self):
        """Test DriveDriveCaller with reference drives (drives with idx)"""
        # Create with one reference drive
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive_with_idx,
            drive_caller2=self.const_drive2
        )
        self.assertEqual(drive_drive.drive_caller1, self.const_drive_with_idx)
        
        # Create with both drives as references
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive_with_idx,
            drive_caller2=self.const_drive_with_idx
        )
        self.assertEqual(drive_drive.drive_caller1, self.const_drive_with_idx)
        self.assertEqual(drive_drive.drive_caller2, self.const_drive_with_idx)

    def test_drive_drive_caller_str_representation(self):
        """Test the string representation of DriveDriveCaller"""
        # Test without idx, with inline drives
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        expected_str = "drive,\n\tconst, 1.5,\n\tconst, 2.5"
        self.assertEqual(str(drive_drive), expected_str)
        
        # Test with idx
        drive_drive = l.DriveDriveCaller(
            idx=10,
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        expected_str = "drive caller: 10, drive,\n\tconst, 1.5,\n\tconst, 2.5"
        self.assertEqual(str(drive_drive), expected_str)
        
        # Test with reference drive for drive_caller1
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive_with_idx,
            drive_caller2=self.const_drive2
        )
        expected_str = "drive,\n\treference, 5,\n\tconst, 2.5"
        self.assertEqual(str(drive_drive), expected_str)
        
        # Test with reference drive for drive_caller2
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive_with_idx
        )
        expected_str = "drive,\n\tconst, 1.5,\n\treference, 5"
        self.assertEqual(str(drive_drive), expected_str)
        
        # Test with reference drives for both
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive_with_idx,
            drive_caller2=self.const_drive_with_idx
        )
        expected_str = "drive,\n\treference, 5,\n\treference, 5"
        self.assertEqual(str(drive_drive), expected_str)

    def test_drive_drive_caller_drive_type(self):
        """Test the drive_type method of DriveDriveCaller"""
        drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        self.assertEqual(drive_drive.drive_type(), 'drive')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_drive_drive_caller_missing_required_field(self):
        """Test creating a DriveDriveCaller instance missing a required field"""
        # Missing drive_caller1
        with self.assertRaises(Exception):
            l.DriveDriveCaller(
                drive_caller2=self.const_drive2
            )
        
        # Missing drive_caller2
        with self.assertRaises(Exception):
            l.DriveDriveCaller(
                drive_caller1=self.const_drive1
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_drive_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for drive_caller1
        with self.assertRaises(Exception):
            l.DriveDriveCaller(
                drive_caller1="not a drive",
                drive_caller2=self.const_drive2
            )
        
        # Invalid type for drive_caller2
        with self.assertRaises(Exception):
            l.DriveDriveCaller(
                drive_caller1=self.const_drive1,
                drive_caller2="not a drive"
            )

    def test_drive_drive_caller_nested(self):
        """Test nesting DriveDriveCallers"""
        # Create a DriveDriveCaller to be used in another DriveDriveCaller
        inner_drive_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.const_drive2
        )
        
        # Use it as drive_caller1 in another DriveDriveCaller
        outer_drive_drive = l.DriveDriveCaller(
            drive_caller1=inner_drive_drive,
            drive_caller2=self.const_drive_with_idx
        )
        
        self.assertIsInstance(outer_drive_drive, l.DriveDriveCaller)
        self.assertEqual(outer_drive_drive.drive_caller1, inner_drive_drive)
        self.assertEqual(outer_drive_drive.drive_caller2, self.const_drive_with_idx)
        
        # Check string representation with nested drive
        expected_str = "drive,\n\tdrive,\n\tconst, 1.5,\n\tconst, 2.5,\n\treference, 5"
        self.assertEqual(str(outer_drive_drive), expected_str)

    def test_drive_drive_caller_complex_nesting(self):
        """Test complex nesting with DriveDriveCaller"""
        # Create two nested DriveDriveCallers
        inner_drive1 = l.DriveDriveCaller(
            drive_caller1=self.const_drive1,
            drive_caller2=self.direct_drive
        )
        
        inner_drive2 = l.DriveDriveCaller(
            drive_caller1=self.const_drive2,
            drive_caller2=self.const_drive_with_idx
        )
        
        # Combine them in an outer DriveDriveCaller
        outer_drive = l.DriveDriveCaller(
            drive_caller1=inner_drive1,
            drive_caller2=inner_drive2
        )
        
        self.assertIsInstance(outer_drive, l.DriveDriveCaller)
        self.assertEqual(outer_drive.drive_caller1, inner_drive1)
        self.assertEqual(outer_drive.drive_caller2, inner_drive2)
        
        # Check the string representation with complex nesting
        expected_str = "drive,\n\tdrive,\n\tconst, 1.5,\n\tdirect,\n\tdrive,\n\tconst, 2.5,\n\treference, 5"
        self.assertEqual(str(outer_drive), expected_str)

class TestElementDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive = l.ConstDriveCaller(const_value=1.5)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=2.0)
        
        # Create concrete Element2 instances for testing
        self.angular_acceleration = l.AngularAcceleration(
            idx=1,
            node_label=101,
            relative_direction=[0.0, 0.0, 1.0],  # Unit vector in z direction
            acceleration=self.const_drive
        )
        
        self.angular_velocity = l.AngularVelocity(
            idx=2,
            node_label=102,
            relative_direction=[0.0, 1.0, 0.0],  # Unit vector in y direction
            velocity=self.const_drive
        )

    def test_element_drive_caller_creation_valid(self):
        """Test that ElementDriveCaller works with valid input"""
        # Create with DriveCaller for func_drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        self.assertIsInstance(element_drive, l.ElementDriveCaller)
        self.assertEqual(element_drive.element, self.angular_acceleration)
        self.assertEqual(element_drive.private_data, "test data")
        self.assertEqual(element_drive.func_drive, self.const_drive)
        
        # Create with 'direct' for func_drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive='direct'
        )
        self.assertEqual(element_drive.func_drive, 'direct')
        
        # Create with specific idx
        element_drive = l.ElementDriveCaller(
            idx=10,
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        self.assertEqual(element_drive.idx, 10)
        
        # Create with a different element
        element_drive = l.ElementDriveCaller(
            element=self.angular_velocity,
            private_data="test data",
            func_drive=self.const_drive
        )
        self.assertEqual(element_drive.element, self.angular_velocity)

    def test_element_drive_caller_str_representation(self):
        """Test the string representation of ElementDriveCaller"""
        # Test with drive_caller for func_drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        expected_str = 'element, 1, joint, string, "test data", const, 1.5'
        self.assertEqual(str(element_drive), expected_str)
        
        # Test with idx
        element_drive = l.ElementDriveCaller(
            idx=10,
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        expected_str = 'drive caller: 10, element, 1, joint, string, "test data", const, 1.5'
        self.assertEqual(str(element_drive), expected_str)
        
        # Test with 'direct' for func_drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive='direct'
        )
        expected_str = 'element, 1, joint, string, "test data", direct'
        self.assertEqual(str(element_drive), expected_str)
        
        # Test with reference drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive_with_idx
        )
        expected_str = 'element, 1, joint, string, "test data", reference, 5'
        self.assertEqual(str(element_drive), expected_str)
        
        # Test with different element type
        element_drive = l.ElementDriveCaller(
            element=self.angular_velocity,
            private_data="test data",
            func_drive=self.const_drive
        )
        expected_str = 'element, 2, joint, string, "test data", const, 1.5'
        self.assertEqual(str(element_drive), expected_str)

    def test_element_drive_caller_drive_type(self):
        """Test the drive_type method of ElementDriveCaller"""
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        self.assertEqual(element_drive.drive_type(), 'element')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_element_drive_caller_missing_required_field(self):
        """Test creating an ElementDriveCaller instance missing a required field"""
        # Missing element
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                private_data="test data",
                func_drive=self.const_drive
            )
        
        # Missing private_data
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                element=self.angular_acceleration,
                func_drive=self.const_drive
            )
        
        # Missing func_drive
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                element=self.angular_acceleration,
                private_data="test data"
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_element_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for element
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                element="not an element",
                private_data="test data",
                func_drive=self.const_drive
            )
        
        # Invalid type for private_data
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                element=self.angular_acceleration,
                private_data=123,  # Should be a string
                func_drive=self.const_drive
            )
        
        # Invalid type for func_drive (not a DriveCaller or 'direct')
        with self.assertRaises(Exception):
            l.ElementDriveCaller(
                element=self.angular_acceleration,
                private_data="test data",
                func_drive="invalid"  # Should be a DriveCaller or 'direct'
            )

    def test_element_drive_caller_nested(self):
        """Test nesting ElementDriveCaller with other drive callers"""
        # Create an ElementDriveCaller to be used as input to another drive
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=self.const_drive
        )
        
        # Create a DriveDriveCaller that uses the ElementDriveCaller
        drive_drive = l.DriveDriveCaller(
            drive_caller1=element_drive,
            drive_caller2=self.const_drive_with_idx
        )
        
        self.assertIsInstance(drive_drive, l.DriveDriveCaller)
        self.assertEqual(drive_drive.drive_caller1, element_drive)
        
        # Check string representation with nested drive
        expected_str = 'drive,\n\telement, 1, joint, string, "test data", const, 1.5,\n\treference, 5'
        self.assertEqual(str(drive_drive), expected_str)

    def test_element_drive_caller_complex_nesting(self):
        """Test complex nesting with ElementDriveCaller"""
        # Create an ElementDriveCaller with a nested drive
        nested_drive = l.DriveDriveCaller(
            drive_caller1=self.const_drive,
            drive_caller2=self.const_drive_with_idx
        )
        
        element_drive = l.ElementDriveCaller(
            element=self.angular_acceleration,
            private_data="test data",
            func_drive=nested_drive
        )
        
        self.assertIsInstance(element_drive, l.ElementDriveCaller)
        self.assertEqual(element_drive.func_drive, nested_drive)
        
        # Check the string representation with complex nesting
        expected_str = 'element, 1, joint, string, "test data", drive,\n\tconst, 1.5,\n\treference, 5'
        self.assertEqual(str(element_drive), expected_str)

class TestExponentialDriveCaller(unittest.TestCase):
    def setUp(self):
        # Create standard values for tests
        self.amplitude = 2.5
        self.time_constant = 1.0
        self.initial_time = 0.0
        self.initial_value = 1.0
        
    def test_exponential_drive_caller_creation_valid(self):
        # Test basic creation with all parameters
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant,
            initial_time=self.initial_time,
            initial_value=self.initial_value
        )
        
        self.assertEqual(drive.amplitude_value, self.amplitude)
        self.assertEqual(drive.time_constant_value, self.time_constant)
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.initial_value, self.initial_value)
        
        # Test creation without optional parameters
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant
        )
        
        self.assertEqual(drive.amplitude_value, self.amplitude)
        self.assertEqual(drive.time_constant_value, self.time_constant)
        self.assertEqual(drive.initial_time, 0.0)
        self.assertEqual(drive.initial_value, 0.0)
        
        # Test with idx parameter
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant,
            initial_time=self.initial_time,
            initial_value=self.initial_value,
            idx=5
        )
        
        self.assertEqual(drive.idx, 5)
        
    def test_exponential_drive_caller_with_mbvars(self):
        # Test with MBVars
        if 'test_amplitude' not in l.declared_MBVars:
            amplitude_var = l.MBVar("test_amplitude", "real", 2.5)
        else:
            amplitude_var = l.declared_MBVars['test_amplitude']
        if 'test_time_constant' not in l.declared_MBVars:
            time_constant_var = l.MBVar("test_time_constant", "real", 1.0)
        else:
            time_constant_var = l.declared_MBVars['test_time_constant']
        if 'test_initial_time' not in l.declared_MBVars:
            initial_time_var = l.MBVar("test_initial_time", "real", 0.5)
        else:
            initial_time_var = l.declared_MBVars['test_initial_time']
        if 'test_initial_value' not in l.declared_MBVars:
            initial_value_var = l.MBVar("test_initial_value", "real", 1.5)
        else:
            initial_value_var = l.declared_MBVars['test_initial_value']
        
        drive = l.ExponentialDriveCaller(
            amplitude_value=amplitude_var,
            time_constant_value=time_constant_var,
            initial_time=initial_time_var,
            initial_value=initial_value_var
        )
        
        self.assertEqual(drive.amplitude_value, amplitude_var)
        self.assertEqual(drive.time_constant_value, time_constant_var)
        self.assertEqual(drive.initial_time, initial_time_var)
        self.assertEqual(drive.initial_value, initial_value_var)
        
    def test_exponential_drive_caller_default_warning(self):
        # Test warnings for default parameters
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered
            warnings.simplefilter("always")
            
            drive = l.ExponentialDriveCaller(
                amplitude_value=self.amplitude,
                time_constant_value=self.time_constant
            )
            
            # Check that two warnings were generated
            self.assertEqual(len(w), 2)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertTrue("<initial_time> is not set, assuming 0.0." in str(w[0].message))
            self.assertTrue(issubclass(w[1].category, UserWarning))
            self.assertTrue("<initial_value> is not set, assuming 0.0." in str(w[1].message))
            
    def test_exponential_drive_caller_str_representation(self):
        # Test string representation without idx
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant,
            initial_time=self.initial_time,
            initial_value=self.initial_value
        )
        
        expected_str = "exponential, 2.5, 1.0, 0.0, 1.0"
        self.assertEqual(str(drive), expected_str)
        
        # Test string representation with idx
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant,
            initial_time=self.initial_time,
            initial_value=self.initial_value,
            idx=5
        )
        
        expected_str = "drive caller: 5, exponential, 2.5, 1.0, 0.0, 1.0"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVars
        amplitude_var = l.MBVar("test_amplitude", "real", 2.5)
        time_constant_var = l.MBVar("test_time_constant", "real", 1.0)
        
        drive = l.ExponentialDriveCaller(
            amplitude_value=amplitude_var,
            time_constant_value=time_constant_var,
            initial_time=self.initial_time,
            initial_value=self.initial_value
        )
        
        expected_str = "exponential, test_amplitude, test_time_constant, 0.0, 1.0"
        self.assertEqual(str(drive), expected_str)
        
    def test_exponential_drive_caller_drive_type(self):
        drive = l.ExponentialDriveCaller(
            amplitude_value=self.amplitude,
            time_constant_value=self.time_constant
        )
        
        self.assertEqual(drive.drive_type(), "exponential")
        
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_exponential_drive_caller_missing_required_field(self):
        # Test missing amplitude_value
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                time_constant_value=self.time_constant
            )
            
        # Test missing time_constant_value
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                amplitude_value=self.amplitude
            )
            
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_exponential_drive_caller_invalid_types(self):
        # Test invalid type for amplitude_value
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                amplitude_value="invalid",
                time_constant_value=self.time_constant
            )
            
        # Test invalid type for time_constant_value
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                amplitude_value=self.amplitude,
                time_constant_value="invalid"
            )
            
        # Test invalid type for initial_time
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                amplitude_value=self.amplitude,
                time_constant_value=self.time_constant,
                initial_time="invalid"
            )
            
        # Test invalid type for initial_value
        with self.assertRaises(Exception):
            l.ExponentialDriveCaller(
                amplitude_value=self.amplitude,
                time_constant_value=self.time_constant,
                initial_value="invalid"
            )

class TestFourierSeriesDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset warnings to make sure we capture them in tests
        warnings.simplefilter("always")
        
        # Create standard values for tests
        self.initial_time = 0.0
        self.angular_velocity = 2.0
        self.number_of_terms = 2
        self.a_0 = 1.0
        self.coefficients = [0.5, 0.6, 0.3, 0.2]  # a_1, b_1, a_2, b_2
        self.number_of_cycles = "forever"
        self.initial_value = 0.0
        
    def test_fourier_series_drive_caller_creation_valid(self):
        """Test that FourierSeriesDriveCaller works with valid input"""
        # Test with all parameters specified
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.angular_velocity, self.angular_velocity)
        self.assertEqual(drive.number_of_terms, self.number_of_terms)
        self.assertEqual(drive.a_0, self.a_0)
        self.assertEqual(drive.coefficients, self.coefficients)
        self.assertEqual(drive.number_of_cycles, self.number_of_cycles)
        self.assertEqual(drive.initial_value, self.initial_value)
        
        # Test with idx
        drive = l.FourierSeriesDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.idx, 5)
        
        # Test with "one" for number_of_cycles
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles="one",
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.number_of_cycles, "one")
        
        # Test with integer for number_of_cycles
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=3,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.number_of_cycles, 3)

    def test_fourier_series_drive_caller_with_mbvars(self):
        """Test with MBVar inputs"""
        # Create MBVars for testing
        if 'init_time' not in l.declared_MBVars:
            init_time_var = l.MBVar(name="init_time", var_type="real", expression=1.5)
        else:
            init_time_var = l.declared_MBVars['init_time']
        if 'ang_vel' not in l.declared_MBVars:
            ang_vel_var = l.MBVar(name="ang_vel", var_type="real", expression=3.14)
        else:
            ang_vel_var = l.declared_MBVars['ang_vel']
        if 'terms' not in l.declared_MBVars:
            terms_var = l.MBVar(name="terms", var_type="integer", expression=2)
        else:
            terms_var = l.declared_MBVars['terms']
        if 'a0' not in l.declared_MBVars:
            a0_var = l.MBVar(name="a0", var_type="real", expression=2.0)
        else:
            a0_var = l.declared_MBVars['a0']
        if 'cycles' not in l.declared_MBVars:
            cycles_var = l.MBVar(name="cycles", var_type="integer", expression=4)
        else:
            cycles_var = l.declared_MBVars['cycles']
        if 'init_val' not in l.declared_MBVars:
            init_val_var = l.MBVar(name="init_val", var_type="real", expression=0.5)
        else:
            init_val_var = l.declared_MBVars['init_val']
        
        # Test with MBVar inputs
        drive = l.FourierSeriesDriveCaller(
            initial_time=init_time_var,
            angular_velocity=ang_vel_var,
            number_of_terms=terms_var,
            a_0=a0_var,
            coefficients=self.coefficients,
            number_of_cycles=cycles_var,
            initial_value=init_val_var
        )
        
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.initial_time, init_time_var)
        self.assertEqual(drive.angular_velocity, ang_vel_var)
        self.assertEqual(drive.number_of_terms, terms_var)
        self.assertEqual(drive.a_0, a0_var)
        self.assertEqual(drive.coefficients, self.coefficients)
        self.assertEqual(drive.number_of_cycles, cycles_var)
        self.assertEqual(drive.initial_value, init_val_var)

    def test_fourier_series_drive_caller_default_warning(self):
        """Test warnings for default parameters"""
        # Test warning for default initial_time
        with warnings.catch_warnings(record=True) as w:
            drive = l.FourierSeriesDriveCaller(
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,#
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(drive.initial_time, 0.0)
        
        # Test warning for default initial_value
        with warnings.catch_warnings(record=True) as w:
            drive = l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles
            )
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(drive.initial_value, 0.0)
        
        # Test both warnings together
        with warnings.catch_warnings(record=True) as w:
            drive = l.FourierSeriesDriveCaller(
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles
            )
            warning_messages = [str(warning.message) for warning in w]
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in msg for msg in warning_messages))
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in msg for msg in warning_messages))
            self.assertEqual(drive.initial_time, 0.0)
            self.assertEqual(drive.initial_value, 0.0)

    def test_fourier_series_drive_caller_str_representation(self):
        """Test string representation"""
        # Test string representation without idx
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        expected_str = (
            "fourier series, 0.0, 2.0, 2,\n"
            "\t1.0,\n"
            "\t0.5, 0.6,\n"
            "\t0.3, 0.2,\n"
            "\tforever, 0.0"
        )
        self.assertEqual(str(drive), expected_str)
        
        # Test string representation with idx
        drive = l.FourierSeriesDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        expected_str = (
            "drive caller: 5, fourier series, 0.0, 2.0, 2,\n"
            "\t1.0,\n"
            "\t0.5, 0.6,\n"
            "\t0.3, 0.2,\n"
            "\tforever, 0.0"
        )
        self.assertEqual(str(drive), expected_str)
        
        # Test string representation with "one" for number_of_cycles
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles="one",
            initial_value=self.initial_value
        )
        expected_str = (
            "fourier series, 0.0, 2.0, 2,\n"
            "\t1.0,\n"
            "\t0.5, 0.6,\n"
            "\t0.3, 0.2,\n"
            "\tone, 0.0"
        )
        self.assertEqual(str(drive), expected_str)

    def test_fourier_series_drive_caller_drive_type(self):
        """Test the drive_type method"""
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=self.number_of_terms,
            a_0=self.a_0,
            coefficients=self.coefficients,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertEqual(drive.drive_type(), "fourier series")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_fourier_series_drive_caller_missing_required_field(self):
        """Test creating a FourierSeriesDriveCaller instance missing a required field"""
        # Missing angular_velocity
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Missing number_of_terms
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Missing a_0
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Missing coefficients
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Missing number_of_cycles
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                initial_value=self.initial_value
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_fourier_series_drive_caller_invalid_types(self):
        """Test with invalid field types"""
        # Test invalid type for angular_velocity
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity="not a number",  # Invalid type
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Test invalid type for number_of_terms
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms="not an integer",  # Invalid type
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Test invalid type for coefficients
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients="not a list",  # Invalid type
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Test invalid MBVar type for initial_time (should be real)
        with self.assertRaises(TypeError):
            if 'invalid_string_mbvar' not in l.declared_MBVars:
                invalid_string_mbvar = l.MBVar(name="invalid_string_mbvar", var_type="string", expression="string value")
            else:
                invalid_string_mbvar = l.declared_MBVars['invalid_string_mbvar']
            l.FourierSeriesDriveCaller(
                initial_time=invalid_string_mbvar,  # Invalid MBVar type
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Test invalid MBVar type for number_of_terms (should be integer)
        with self.assertRaises(TypeError):
            if 'invalid_real_mbvar' not in l.declared_MBVars:
                invalid_real_mbvar = l.MBVar(name="invalid_real_mbvar", var_type="real", expression=2.5)
            else:
                invalid_real_mbvar = l.declared_MBVars['invalid_real_mbvar']
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=invalid_real_mbvar,  # Invalid MBVar type
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )

        # Test invalid type for number_of_cycles
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=self.coefficients,
                number_of_cycles=2.5,  # Invalid type
                initial_value=self.initial_value
            )

        # Test invalid type for a_0
        with self.assertRaises(Exception):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0='invalid',  # Invalid type
                coefficients=self.coefficients,
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_coefficients_validation(self):
        """Test validation of coefficients list length"""
        # Test with incorrect number of coefficients
        with self.assertRaises(ValueError):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=2,
                a_0=self.a_0,
                coefficients=[0.5, 0.6, 0.3],  # Should be 4 elements for 2 terms
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )
            
        # Test with correct number of coefficients for different number_of_terms
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=1,
            a_0=self.a_0,
            coefficients=[0.5, 0.6],  # Correct for 1 term (a_1, b_1)
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        
        # Test with 3 terms
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=3,
            a_0=self.a_0,
            coefficients=[0.5, 0.6, 0.3, 0.2, 0.1, 0.05],  # 6 elements for 3 terms
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)

    def test_fourier_series_with_mbvar_coefficients(self):
        """Test validation of coefficients when using MBVars"""
        # Create MBVars for coefficients
        if 'coef_a1' not in l.declared_MBVars:
            coef_a1 = l.MBVar(name="coef_a1", var_type="real", expression=0.75)
        else:
            coef_a1 = l.declared_MBVars['coef_a1']
            
        if 'coef_b1' not in l.declared_MBVars:
            coef_b1 = l.MBVar(name="coef_b1", var_type="real", expression=0.25)
        else:
            coef_b1 = l.declared_MBVars['coef_b1']
        
        # Test with MBVars in the coefficients list
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=1,
            a_0=self.a_0,
            coefficients=[coef_a1, coef_b1],  # Using MBVars as coefficients
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        self.assertEqual(drive.coefficients[0], coef_a1)
        self.assertEqual(drive.coefficients[1], coef_b1)
        
        # Test with mixed regular values and MBVars
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=2,
            a_0=self.a_0,
            coefficients=[coef_a1, coef_b1, 0.3, 0.2],  # Mix of MBVars and floats
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        
    def test_fourier_series_with_mbvar_number_of_terms(self):
        """Test validation of coefficients when number_of_terms is an MBVar"""
        # Create MBVar for number_of_terms
        if 'terms_var' not in l.declared_MBVars:
            terms_var = l.MBVar(name="terms_var", var_type="integer", expression=2)
        else:
            terms_var = l.declared_MBVars['terms_var']
        
        # Test with correct number of coefficients for MBVar terms
        drive = l.FourierSeriesDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            number_of_terms=terms_var,
            a_0=self.a_0,
            coefficients=[0.5, 0.6, 0.3, 0.2],  # 4 coefficients for 2 terms
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertIsInstance(drive, l.FourierSeriesDriveCaller)
        
        # Test with incorrect number of coefficients
        with self.assertRaises(ValueError):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=terms_var,
                a_0=self.a_0,
                coefficients=[0.5, 0.6, 0.3],  # Only 3 coefficients instead of 4
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )

    def test_fourier_series_with_invalid_mbvar_coefficients(self):
        """Test validation fails with invalid MBVar types in coefficients"""
        # Create invalid MBVar (string type instead of real)
        if 'invalid_coef' not in l.declared_MBVars:
            invalid_coef = l.MBVar(name="invalid_coef", var_type="integer", expression=10)
        else:
            invalid_coef = l.declared_MBVars['invalid_coef']
        
        # Test with invalid MBVar type in coefficients
        with self.assertRaises(TypeError):
            l.FourierSeriesDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                number_of_terms=self.number_of_terms,
                a_0=self.a_0,
                coefficients=[0.5, 0.6, invalid_coef, 0.2],  # Invalid MBVar type
                number_of_cycles=self.number_of_cycles,
                initial_value=self.initial_value
            )

class TestFrequencySweepDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive1 = l.ConstDriveCaller(const_value=2.0)
        self.const_drive2 = l.ConstDriveCaller(const_value=3.0)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=4.0)
        
        # Create MBVar objects for testing
        if 'test_initial_time' not in l.declared_MBVars:
            self.initial_time_var = l.MBVar(name='test_initial_time', var_type='real', expression=1.5)
        else:
            self.initial_time_var = l.declared_MBVars['test_initial_time']
            
        if 'test_final_value' not in l.declared_MBVars:
            self.final_value_var = l.MBVar(name='test_final_value', var_type='real', expression=7.5)
        else:
            self.final_value_var = l.declared_MBVars['test_final_value']

    def test_frequency_sweep_drive_caller_creation_valid(self):
        """Test that FrequencySweepDriveCaller works with valid input"""
        # Create with all required parameters (initial_time and initial_value are optional with defaults)
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=5.0
        )
        self.assertIsInstance(freq_sweep_drive, l.FrequencySweepDriveCaller)
        self.assertEqual(freq_sweep_drive.initial_time, 0.0)  # Default value
        self.assertEqual(freq_sweep_drive.angular_velocity_drive, self.const_drive1)
        self.assertEqual(freq_sweep_drive.amplitude_drive, self.const_drive2)
        self.assertEqual(freq_sweep_drive.initial_value, 0.0)  # Default value
        self.assertEqual(freq_sweep_drive.final_time, 10.0)
        self.assertEqual(freq_sweep_drive.final_value, 5.0)
        
        # Create with explicit initial_time and initial_value
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=2.0,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time=10.0,
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.initial_time, 2.0)
        self.assertEqual(freq_sweep_drive.initial_value, 1.0)
        
        # Create with 'forever' for final_time
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time='forever',
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.final_time, 'forever')
        
        # Create with specific idx
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            idx=10,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.idx, 10)
        
        # Create with drive callers that have idx
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=self.const_drive_with_idx,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.angular_velocity_drive, self.const_drive_with_idx)

    def test_frequency_sweep_drive_caller_with_mbvars(self):
        """Test FrequencySweepDriveCaller with MBVar objects for parameters"""
        # Create with MBVar for initial_time
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=self.initial_time_var,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.initial_time, self.initial_time_var)
        
        # Create with MBVar for final_value
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=self.final_value_var
        )
        self.assertEqual(freq_sweep_drive.final_value, self.final_value_var)

    def test_frequency_sweep_drive_caller_default_warning(self):
        """Test that a warning is issued when initial_time or initial_value isn't provided"""
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered
            warnings.simplefilter("always")
            
            # Create a FrequencySweepDriveCaller without specifying initial_time or initial_value
            freq_sweep_drive = l.FrequencySweepDriveCaller(
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_time=10.0,
                final_value=5.0
            )
            
            # Verify warnings were raised
            self.assertEqual(len(w), 2)
            warning_messages = [str(warning.message) for warning in w]
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in msg for msg in warning_messages))
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in msg for msg in warning_messages))

        # No warnings when initial_time and initial_value are explicitly provided
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            freq_sweep_drive = l.FrequencySweepDriveCaller(
                initial_time=2.0,
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                initial_value=1.0,
                final_time=10.0,
                final_value=5.0
            )
            
            self.assertEqual(len(w), 0)  # No warnings

    def test_frequency_sweep_drive_caller_str_representation(self):
        """Test the string representation of FrequencySweepDriveCaller"""
        # Test without idx
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=2.0,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time=10.0,
            final_value=5.0
        )
        expected_str = "frequency sweep, 2.0,\n\tconst, 2.0,\n\tconst, 3.0\n\t1.0, 10.0, 5.0"
        self.assertEqual(str(freq_sweep_drive), expected_str)
        
        # Test with idx
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            idx=10,
            initial_time=2.0,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time=10.0,
            final_value=5.0
        )
        expected_str = "drive caller: 10, frequency sweep, 2.0,\n\tconst, 2.0,\n\tconst, 3.0\n\t1.0, 10.0, 5.0"
        self.assertEqual(str(freq_sweep_drive), expected_str)
        
        # Test with drive callers that have idx
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=2.0,
            angular_velocity_drive=self.const_drive_with_idx,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time=10.0,
            final_value=5.0
        )
        expected_str = "frequency sweep, 2.0,\n\treference, 5,\n\tconst, 3.0\n\t1.0, 10.0, 5.0"
        self.assertEqual(str(freq_sweep_drive), expected_str)
        
        # Test with MBVar parameters
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=self.initial_time_var,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time=10.0,
            final_value=self.final_value_var
        )
        expected_str = f"frequency sweep, {self.initial_time_var},\n\tconst, 2.0,\n\tconst, 3.0\n\t1.0, 10.0, {self.final_value_var}"
        self.assertEqual(str(freq_sweep_drive), expected_str)
        
        # Test with 'forever'
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            initial_time=2.0,
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            initial_value=1.0,
            final_time='forever',
            final_value=5.0
        )
        expected_str = "frequency sweep, 2.0,\n\tconst, 2.0,\n\tconst, 3.0\n\t1.0, forever, 5.0"
        self.assertEqual(str(freq_sweep_drive), expected_str)

    def test_frequency_sweep_drive_caller_drive_type(self):
        """Test the drive_type method of FrequencySweepDriveCaller"""
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=self.const_drive1,
            amplitude_drive=self.const_drive2,
            final_time=10.0,
            final_value=5.0
        )
        self.assertEqual(freq_sweep_drive.drive_type(), 'frequency sweep')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_frequency_sweep_drive_caller_missing_required_field(self):
        """Test creating a FrequencySweepDriveCaller instance missing a required field"""
        # Missing angular_velocity_drive
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                amplitude_drive=self.const_drive2,
                final_time=10.0,
                final_value=5.0
            )
        
        # Missing amplitude_drive
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                angular_velocity_drive=self.const_drive1,
                final_time=10.0,
                final_value=5.0
            )
        
        # Missing final_time
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_value=5.0
            )
        
        # Missing final_value
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_time=10.0,
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_frequency_sweep_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for initial_time
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                initial_time="invalid",
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_time=10.0,
                final_value=5.0
            )
        
        # Invalid type for angular_velocity_drive
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                angular_velocity_drive="not a drive",
                amplitude_drive=self.const_drive2,
                final_time=10.0,
                final_value=5.0
            )
        
        # Invalid type for final_time (should be float, MBVar, or 'forever')
        with self.assertRaises(Exception):
            l.FrequencySweepDriveCaller(
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_time="invalid",
                final_value=5.0
            )

    def test_frequency_sweep_drive_caller_nested(self):
        """Test nesting FrequencySweepDriveCaller with other drive callers"""
        # Create array drives to use as angular_velocity_drive and amplitude_drive
        array_drive1 = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive2])
        array_drive2 = l.ArrayDriveCaller(drives=[self.const_drive2, self.const_drive_with_idx])
        
        # Use them in a FrequencySweepDriveCaller
        freq_sweep_drive = l.FrequencySweepDriveCaller(
            angular_velocity_drive=array_drive1,
            amplitude_drive=array_drive2,
            final_time=10.0,
            final_value=5.0
        )
        
        self.assertIsInstance(freq_sweep_drive, l.FrequencySweepDriveCaller)
        self.assertEqual(freq_sweep_drive.angular_velocity_drive, array_drive1)
        self.assertEqual(freq_sweep_drive.amplitude_drive, array_drive2)
        
        # Check string representation with nested drives
        expected_str = "frequency sweep, 0.0,\n\tarray, 2,\n\tconst, 2.0,\n\tconst, 3.0,\n\tarray, 2,\n\tconst, 3.0,\n\treference, 5\n\t0.0, 10.0, 5.0"
        self.assertEqual(str(freq_sweep_drive), expected_str)

    def test_field_validator_for_real_mbvar(self):
        """Test the field validator for real MBVar fields"""
        # Create an MBVar with a non-real type
        if 'test_non_real' not in l.declared_MBVars:
            non_real_var = l.MBVar(name='test_non_real', var_type='integer', expression=5)
        else:
            non_real_var = l.declared_MBVars['test_non_real']
        
        # Test that using a non-real MBVar for initial_time raises an error
        with self.assertRaises(TypeError) as context:
            l.FrequencySweepDriveCaller(
                initial_time=non_real_var,  # This should be a real MBVar or float
                angular_velocity_drive=self.const_drive1,
                amplitude_drive=self.const_drive2,
                final_time=10.0,
                final_value=5.0
            )
        self.assertIn("Field must be an MBVar of type real or a float", str(context.exception))

class TestGiNaCDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset warnings to make sure we capture them in tests
        warnings.resetwarnings()
        
        # Create sample expressions for testing
        self.expression_str = "x^2 + sin(x)"
        self.symbol_str = "x"
        
        # Create MBVars for testing
        if pydantic is not None:  # Only create if validation is available
            if 'expr_var' not in l.declared_MBVars:
                self.expr_var=l.MBVar("expr_var", "string", "x^2 + 2*x + 1")
            else:
                self.expr_var=l.declared_MBVars['expr_var'] = l.MBVar("expr_var", "string", "x^2 + 2*x + 1")
            if 'sym_var' not in l.declared_MBVars:
                self.sym_var=l.MBVar("sym_var", "string", "x")
            else:
                self.sym_var=l.declared_MBVars['sym_var'] = l.MBVar("sym_var", "string", "x")

    def test_ginac_drive_caller_creation_valid(self):
        """Test that GiNaCDriveCaller works with valid input"""
        # Create with just expression (string)
        ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str)
        self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
        self.assertEqual(ginac_drive.expression, self.expression_str)
        self.assertIsNone(ginac_drive.symbol)
        
        # Create with expression and symbol (both strings)
        ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str, symbol=self.symbol_str)
        self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
        self.assertEqual(ginac_drive.expression, self.expression_str)
        self.assertEqual(ginac_drive.symbol, self.symbol_str)
        
        # Create with specific idx
        ginac_drive = l.GiNaCDriveCaller(idx=10, expression=self.expression_str)
        self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
        self.assertEqual(ginac_drive.idx, 10)
        self.assertEqual(ginac_drive.expression, self.expression_str)
        
        if pydantic is not None:  # Only test with MBVars if validation is available
            # Create with expression as MBVar
            ginac_drive = l.GiNaCDriveCaller(expression=self.expr_var)
            self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
            self.assertEqual(ginac_drive.expression, self.expr_var)
            
            # Create with symbol as MBVar
            ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str, symbol=self.sym_var)
            self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
            self.assertEqual(ginac_drive.symbol, self.sym_var)
            
            # Create with both expression and symbol as MBVars
            ginac_drive = l.GiNaCDriveCaller(expression=self.expr_var, symbol=self.sym_var)
            self.assertIsInstance(ginac_drive, l.GiNaCDriveCaller)
            self.assertEqual(ginac_drive.expression, self.expr_var)
            self.assertEqual(ginac_drive.symbol, self.sym_var)

    def test_ginac_drive_caller_str_representation(self):
        """Test the string representation of GiNaCDriveCaller"""
        # Test without idx or symbol
        ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str)
        expected_str = f'ginac, "{self.expression_str}"'
        self.assertEqual(str(ginac_drive), expected_str)
        
        # Test with idx but no symbol
        ginac_drive = l.GiNaCDriveCaller(idx=10, expression=self.expression_str)
        expected_str = f'drive caller: 10, ginac, "{self.expression_str}"'
        self.assertEqual(str(ginac_drive), expected_str)
        
        # Test with symbol but no idx
        ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str, symbol=self.symbol_str)
        expected_str = f'ginac, symbol, "{self.symbol_str}", "{self.expression_str}"'
        self.assertEqual(str(ginac_drive), expected_str)
        
        # Test with both idx and symbol
        ginac_drive = l.GiNaCDriveCaller(idx=10, expression=self.expression_str, symbol=self.symbol_str)
        expected_str = f'drive caller: 10, ginac, symbol, "{self.symbol_str}", "{self.expression_str}"'
        self.assertEqual(str(ginac_drive), expected_str)
        
        if pydantic is not None:  # Only test with MBVars if validation is available
            # Test with expression as MBVar
            ginac_drive = l.GiNaCDriveCaller(expression=self.expr_var)
            expected_str = f'ginac, {self.expr_var}'
            self.assertEqual(str(ginac_drive), expected_str)
            
            # Test with symbol as MBVar
            ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str, symbol=self.sym_var)
            expected_str = f'ginac, symbol, {self.sym_var}, "{self.expression_str}"'
            self.assertEqual(str(ginac_drive), expected_str)
            
            # Test with both expression and symbol as MBVars
            ginac_drive = l.GiNaCDriveCaller(expression=self.expr_var, symbol=self.sym_var)
            expected_str = f'ginac, symbol, {self.sym_var}, {self.expr_var}'
            self.assertEqual(str(ginac_drive), expected_str)

    def test_ginac_drive_caller_drive_type(self):
        """Test the drive_type method returns the correct string"""
        ginac_drive = l.GiNaCDriveCaller(expression=self.expression_str)
        self.assertEqual(ginac_drive.drive_type(), "ginac")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_ginac_drive_caller_missing_required_field(self):
        """Test creating a GiNaCDriveCaller instance missing a required field"""
        # Test missing expression
        with self.assertRaises(Exception):
            l.GiNaCDriveCaller()  # Missing expression field
        
        # Test with only idx (missing expression)
        with self.assertRaises(Exception):
            l.GiNaCDriveCaller(idx=10)  # Missing expression field
        
        # Test with only symbol (missing expression)
        with self.assertRaises(Exception):
            l.GiNaCDriveCaller(symbol=self.symbol_str)  # Missing expression field

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_ginac_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Test invalid type for expression (not string or string MBVar)
        with self.assertRaises(TypeError):
            if 'invalid_expr' not in l.declared_MBVars:
                invalid_expr = l.MBVar("invalid_expr", "integer", 42)
            else:
                invalid_expr = l.declared_MBVars['invalid_expr']
            l.GiNaCDriveCaller(expression=invalid_expr)  # Should be string MBVar
        
        with self.assertRaises(Exception):
            l.GiNaCDriveCaller(expression=42)  # Should be string
        
        # Test invalid type for symbol (not string or string MBVar)
        with self.assertRaises(TypeError):
            if 'invalid_sym' not in l.declared_MBVars:
                invalid_sym = l.MBVar("invalid_sym", "integer", 42)
            else:
                invalid_sym = l.declared_MBVars['invalid_sym']
            l.GiNaCDriveCaller(expression=self.expression_str, symbol=invalid_sym)  # Should be string MBVar
        
        with self.assertRaises(Exception):
            l.GiNaCDriveCaller(expression=self.expression_str, symbol=42)  # Should be string

    def test_ginac_drive_caller_complex_expressions(self):
        """Test with more complex mathematical expressions"""
        # Test complex expression
        complex_expr = "sin(x)^2 + cos(x)^2 + tan(x/2) + sqrt(abs(x))"
        ginac_drive = l.GiNaCDriveCaller(expression=complex_expr)
        self.assertEqual(ginac_drive.expression, complex_expr)
        expected_str = f'ginac, "{complex_expr}"'
        self.assertEqual(str(ginac_drive), expected_str)
        
        # Test expression with different variable
        ginac_drive = l.GiNaCDriveCaller(expression="a*t^2 + b*t + c", symbol="t")
        expected_str = 'ginac, symbol, "t", "a*t^2 + b*t + c"'
        self.assertEqual(str(ginac_drive), expected_str)

class TestLinearDriveCaller(unittest.TestCase):
    def setUp(self):
        # Setup common values for testing
        self.const_coef = 2.5
        self.slope_coef = 1.5

        # Create MBVar objects
        if 'const_coef_var' not in l.declared_MBVars:
            self.const_coef_var = l.MBVar(name='const_coef_var', var_type='real', expression=3.0)
        else:
            self.const_coef_var = l.declared_MBVars['const_coef_var']
            
        if 'slope_coef_var' not in l.declared_MBVars:
            self.slope_coef_var = l.MBVar(name='slope_coef_var', var_type='real', expression=2.0)
        else:
            self.slope_coef_var = l.declared_MBVars['slope_coef_var']
    
    def test_linear_drive_caller_creation_valid(self):
        """Test creating a LinearDriveCaller with valid parameters"""
        # Basic creation with all parameters
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef,
            slope_coef=self.slope_coef
        )
        
        # Verify properties are set correctly
        self.assertEqual(drive.const_coef, self.const_coef)
        self.assertEqual(drive.slope_coef, self.slope_coef)
        self.assertEqual(drive.drive_type(), "linear")
        
        # Test with idx parameter
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef,
            slope_coef=self.slope_coef,
            idx=5
        )
        
        self.assertEqual(drive.idx, 5)
        
    def test_linear_drive_caller_with_mbvars(self):
        """Test creating a LinearDriveCaller with MBVar objects"""        
        # Create with MBVar objects
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef_var,
            slope_coef=self.slope_coef_var
        )
        
        # Check that MBVar references are stored correctly
        self.assertEqual(drive.const_coef, self.const_coef_var)
        self.assertEqual(drive.slope_coef, self.slope_coef_var)
        
        # Create with mixed types (MBVar and float)
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef_var,
            slope_coef=self.slope_coef
        )
        
        self.assertEqual(drive.const_coef, self.const_coef_var)
        self.assertEqual(drive.slope_coef, self.slope_coef)
        
    def test_linear_drive_caller_str_representation(self):
        """Test the string representation of LinearDriveCaller"""
        # Test without idx
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef,
            slope_coef=self.slope_coef
        )
        
        expected_str = f"linear, {self.const_coef}, {self.slope_coef}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef,
            slope_coef=self.slope_coef,
            idx=5
        )
        
        expected_str = f"drive caller: 5, linear, {self.const_coef}, {self.slope_coef}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVars        
        drive = l.LinearDriveCaller(
            const_coef=self.const_coef_var,
            slope_coef=self.slope_coef_var
        )
        
        expected_str = f"linear, {self.const_coef_var}, {self.slope_coef_var}"
        self.assertEqual(str(drive), expected_str)
        
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing const_coef
        with self.assertRaises(Exception):
            l.LinearDriveCaller(
                slope_coef=self.slope_coef
            )
        
        # Missing slope_coef
        with self.assertRaises(Exception):
            l.LinearDriveCaller(
                const_coef=self.const_coef
            )
            
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # const_coef as string
        with self.assertRaises(ValueError):
            l.LinearDriveCaller(
                const_coef="invalid",
                slope_coef=self.slope_coef
            )
        
        # slope_coef as string
        with self.assertRaises(ValueError):
            l.LinearDriveCaller(
                const_coef=self.const_coef,
                slope_coef="invalid"
            )
        
        # Create an MBVar with a non-real type
        if 'non_real_var' not in l.declared_MBVars:
            non_real_var = l.MBVar(name='non_real_var', var_type='integer', expression=5)
        else:
            non_real_var = l.declared_MBVars['non_real_var']
        
        # Test with non-real MBVar for const_coef
        with self.assertRaises(TypeError) as context:
            l.LinearDriveCaller(
                const_coef=non_real_var,
                slope_coef=self.slope_coef
            )
        self.assertIn("coefficients must be real numbers or MBVars of type real", str(context.exception))
        
        # Test with non-real MBVar for slope_coef
        with self.assertRaises(TypeError) as context:
            l.LinearDriveCaller(
                const_coef=self.const_coef,
                slope_coef=non_real_var
            )
        self.assertIn("coefficients must be real numbers or MBVars of type real", str(context.exception))
        
    def test_linear_drive_caller_in_other_drives(self):
        """Test using LinearDriveCaller as input to other drive callers"""
        # Create a LinearDriveCaller
        linear_drive = l.LinearDriveCaller(
            const_coef=self.const_coef,
            slope_coef=self.slope_coef
        )
        
        # Use it in a MultDriveCaller
        const_drive = l.ConstDriveCaller(const_value=3.0)
        mult_drive = l.MultDriveCaller(
            drive_1=linear_drive,
            drive_2=const_drive
        )
        
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, linear_drive)
        self.assertEqual(mult_drive.drive_2, const_drive)
        
        # Check string representation
        expected_str = f"mult,\n\tlinear, {self.const_coef}, {self.slope_coef},\n\tconst, 3.0"
        self.assertEqual(str(mult_drive), expected_str)

class TestMeterDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset warnings to make sure we capture them in tests
        warnings.resetwarnings()
        # Setup common values for testing
        self.initial_time = 1.0
        self.final_time = 10.0
        self.steps_between_spikes = 5
    
    def test_meter_drive_caller_creation_valid(self):
        """Test creating a MeterDriveCaller with valid parameters"""
        # Create with all parameters
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time,
            steps_between_spikes=self.steps_between_spikes
        )
        
        # Verify all properties are set correctly
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.final_time, self.final_time)
        self.assertEqual(drive.steps_between_spikes, self.steps_between_spikes)
        self.assertEqual(drive.drive_type(), "meter")
        
        # Test without optional steps_between_spikes
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time
        )
        
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.final_time, self.final_time)
        self.assertIsNone(drive.steps_between_spikes)
        
        # Test with 'forever' for final_time
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time='forever',
            steps_between_spikes=self.steps_between_spikes
        )
        
        self.assertEqual(drive.final_time, 'forever')
        
        # Test with idx parameter
        drive = l.MeterDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            final_time=self.final_time,
            steps_between_spikes=self.steps_between_spikes
        )
        
        self.assertEqual(drive.idx, 5)
    
    def test_meter_drive_caller_default_values(self):
        """Test that default values are set correctly with warnings"""
        with warnings.catch_warnings(record=True) as w:
            # Create drive without initial_time
            drive = l.MeterDriveCaller(
                final_time=self.final_time,
                steps_between_spikes=self.steps_between_spikes
            )
            
            # Check default value
            self.assertEqual(drive.initial_time, 0.0)
            
            # Verify warnings were raised
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertTrue("<initial_time> is not set, assuming 0.0." in str(w[0].message))
    
    def test_meter_drive_caller_with_mbvars(self):
        """Test creating a MeterDriveCaller with MBVar objects"""
        # Create MBVar objects
        if 'meter_init_time' not in l.declared_MBVars:
            initial_time_var = l.MBVar(name='meter_init_time', var_type='real', expression=1.5)
        else:
            initial_time_var = l.declared_MBVars['meter_init_time']
            
        if 'meter_final_time' not in l.declared_MBVars:
            final_time_var = l.MBVar(name='meter_final_time', var_type='real', expression=12.0)
        else:
            final_time_var = l.declared_MBVars['meter_final_time']
            
        if 'meter_steps' not in l.declared_MBVars:
            steps_var = l.MBVar(name='meter_steps', var_type='integer', expression=8)
        else:
            steps_var = l.declared_MBVars['meter_steps']
        
        # Create with all MBVar objects
        drive = l.MeterDriveCaller(
            initial_time=initial_time_var,
            final_time=final_time_var,
            steps_between_spikes=steps_var
        )
        
        # Check that MBVar references are stored correctly
        self.assertEqual(drive.initial_time, initial_time_var)
        self.assertEqual(drive.final_time, final_time_var)
        self.assertEqual(drive.steps_between_spikes, steps_var)
        
        # Create with mixed parameters (some MBVars, some values)
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=final_time_var,
            steps_between_spikes=self.steps_between_spikes
        )
        
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.final_time, final_time_var)
        self.assertEqual(drive.steps_between_spikes, self.steps_between_spikes)
    
    def test_meter_drive_caller_str_representation(self):
        """Test the string representation of MeterDriveCaller"""
        # Test with all parameters
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time,
            steps_between_spikes=self.steps_between_spikes
        )
        
        expected_str = f"meter, {self.initial_time}, {self.final_time}, steps, {self.steps_between_spikes}"
        self.assertEqual(str(drive), expected_str)
        
        # Test without steps_between_spikes
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time
        )
        
        expected_str = f"meter, {self.initial_time}, {self.final_time}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.MeterDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            final_time=self.final_time,
            steps_between_spikes=self.steps_between_spikes
        )
        
        expected_str = f"drive caller: 5, meter, {self.initial_time}, {self.final_time}, steps, {self.steps_between_spikes}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with 'forever'
        drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time='forever',
            steps_between_spikes=self.steps_between_spikes
        )
        
        expected_str = f"meter, {self.initial_time}, forever, steps, {self.steps_between_spikes}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVars
        initial_time_var = l.MBVar(name='meter_init_time', var_type='real', expression=1.5)
        final_time_var = l.MBVar(name='meter_final_time', var_type='real', expression=12.0)
        
        drive = l.MeterDriveCaller(
            initial_time=initial_time_var,
            final_time=final_time_var,
            steps_between_spikes=self.steps_between_spikes
        )
        
        expected_str = f"meter, {initial_time_var}, {final_time_var}, steps, {self.steps_between_spikes}"
        self.assertEqual(str(drive), expected_str)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_meter_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing final_time
        with self.assertRaises(Exception):
            l.MeterDriveCaller(
                initial_time=self.initial_time
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_meter_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # Invalid type for initial_time
        with self.assertRaises(Exception):
            l.MeterDriveCaller(
                initial_time="invalid",
                final_time=self.final_time
            )
        
        # Invalid string for final_time (not 'forever')
        with self.assertRaises(Exception):
            l.MeterDriveCaller(
                initial_time=self.initial_time,
                final_time="invalid"
            )
        
        # Invalid (non-positive) value for steps_between_spikes
        with self.assertRaises(ValueError) as context:
            l.MeterDriveCaller(
                initial_time=self.initial_time,
                final_time=self.final_time,
                steps_between_spikes=0  # Should be positive
            )
        self.assertIn("steps_between_spikes must be a positive integer", str(context.exception))
        
        # Create MBVars with wrong types
        if 'non_real_var' not in l.declared_MBVars:
            non_real_var = l.MBVar(name='non_real_var', var_type='integer', expression=5)
        else:
            non_real_var = l.declared_MBVars['non_real_var']
            
        if 'non_integer_var' not in l.declared_MBVars:
            non_integer_var = l.MBVar(name='non_integer_var', var_type='real', expression=5.0)
        else:
            non_integer_var = l.declared_MBVars['non_integer_var']
        
        # Test with non-real MBVar for initial_time
        with self.assertRaises(TypeError) as context:
            l.MeterDriveCaller(
                initial_time=non_real_var,
                final_time=self.final_time
            )
        self.assertIn("Field must be a real number or an MBVar of type real", str(context.exception))
        
        # Test with non-real MBVar for final_time
        with self.assertRaises(TypeError) as context:
            l.MeterDriveCaller(
                initial_time=self.initial_time,
                final_time=non_real_var
            )
        self.assertIn("Field must be a real number or an MBVar of type real", str(context.exception))
        
        # Test with non-integer MBVar for steps_between_spikes
        with self.assertRaises(TypeError) as context:
            l.MeterDriveCaller(
                initial_time=self.initial_time,
                final_time=self.final_time,
                steps_between_spikes=non_integer_var
            )
        self.assertIn("steps_between_spikes must be an integer or an MBVar of type integer", str(context.exception))
    
    def test_meter_drive_caller_in_other_drives(self):
        """Test using MeterDriveCaller as input to other drive callers"""
        # Create a MeterDriveCaller
        meter_drive = l.MeterDriveCaller(
            initial_time=self.initial_time,
            final_time=self.final_time,
            steps_between_spikes=self.steps_between_spikes
        )
        
        # Use it in a MultDriveCaller
        const_drive = l.ConstDriveCaller(const_value=2.0)
        mult_drive = l.MultDriveCaller(
            drive_1=meter_drive,
            drive_2=const_drive
        )
        
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, meter_drive)
        self.assertEqual(mult_drive.drive_2, const_drive)
        
        # Check string representation
        expected_mult_str = f"mult,\n\tmeter, {self.initial_time}, {self.final_time}, steps, {self.steps_between_spikes},\n\tconst, 2.0"
        self.assertEqual(str(mult_drive), expected_mult_str)
        
        # Use it in an ArrayDriveCaller
        array_drive = l.ArrayDriveCaller(drives=[meter_drive, const_drive])
        
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(len(array_drive.drives), 2)
        self.assertEqual(array_drive.drives[0], meter_drive)
        
        # Check string representation
        expected_array_str = f"array, 2,\n\tmeter, {self.initial_time}, {self.final_time}, steps, {self.steps_between_spikes},\n\tconst, 2.0"
        self.assertEqual(str(array_drive), expected_array_str)

class TestMultDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create sample drive callers for testing
        self.const_drive1 = l.ConstDriveCaller(const_value=3.0)
        self.const_drive2 = l.ConstDriveCaller(const_value=4.0)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=5.0)

    def test_mult_drive_caller_creation_valid(self):
        """Test that MultDriveCaller works with valid input"""
        # Create with two drives without idx
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive2)
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, self.const_drive1)
        self.assertEqual(mult_drive.drive_2, self.const_drive2)
        
        # Create with one drive having idx
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive_with_idx, drive_2=self.const_drive2)
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, self.const_drive_with_idx)
        self.assertEqual(mult_drive.drive_2, self.const_drive2)
        
        # Create with both drives having idx
        other_drive_with_idx = l.ConstDriveCaller(idx=6, const_value=6.0)
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive_with_idx, drive_2=other_drive_with_idx)
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, self.const_drive_with_idx)
        self.assertEqual(mult_drive.drive_2, other_drive_with_idx)
        
        # Create with specific idx for the mult drive
        mult_drive = l.MultDriveCaller(idx=10, drive_1=self.const_drive1, drive_2=self.const_drive2)
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.idx, 10)
        self.assertEqual(mult_drive.drive_1, self.const_drive1)
        self.assertEqual(mult_drive.drive_2, self.const_drive2)

    def test_mult_drive_caller_str_representation(self):
        """Test the string representation of MultDriveCaller"""
        # Test without idx
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive2)
        expected_str = "mult,\n\tconst, 3.0,\n\tconst, 4.0"
        self.assertEqual(str(mult_drive), expected_str)
        
        # Test with idx
        mult_drive = l.MultDriveCaller(idx=10, drive_1=self.const_drive1, drive_2=self.const_drive2)
        expected_str = "drive caller: 10, mult,\n\tconst, 3.0,\n\tconst, 4.0"
        self.assertEqual(str(mult_drive), expected_str)
        
        # Test with drive_1 having idx
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive_with_idx, drive_2=self.const_drive2)
        expected_str = "mult,\n\treference, 5,\n\tconst, 4.0"
        self.assertEqual(str(mult_drive), expected_str)
        
        # Test with drive_2 having idx
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive_with_idx)
        expected_str = "mult,\n\tconst, 3.0,\n\treference, 5"
        self.assertEqual(str(mult_drive), expected_str)
        
        # Test with both drives having idx
        other_drive_with_idx = l.ConstDriveCaller(idx=6, const_value=6.0)
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive_with_idx, drive_2=other_drive_with_idx)
        expected_str = "mult,\n\treference, 5,\n\treference, 6"
        self.assertEqual(str(mult_drive), expected_str)

    def test_mult_drive_caller_drive_type(self):
        """Test the drive_type method of MultDriveCaller"""
        mult_drive = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive2)
        self.assertEqual(mult_drive.drive_type(), 'mult')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_mult_drive_caller_missing_required_field(self):
        """Test creating a MultDriveCaller instance missing a required field"""
        # Missing drive_1
        with self.assertRaises(Exception):
            l.MultDriveCaller(drive_2=self.const_drive2)
        
        # Missing drive_2
        with self.assertRaises(Exception):
            l.MultDriveCaller(drive_1=self.const_drive1)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_mult_drive_caller_invalid_types(self):
        """Test invalid types for fields"""
        # Invalid type for drive_1
        with self.assertRaises(Exception):
            l.MultDriveCaller(drive_1="not a drive", drive_2=self.const_drive2)
        
        # Invalid type for drive_2
        with self.assertRaises(Exception):
            l.MultDriveCaller(drive_1=self.const_drive1, drive_2="not a drive")

    def test_mult_drive_caller_nested(self):
        """Test nesting MultDriveCaller with other drive callers"""
        # Create a nested mult drive caller
        inner_mult = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive2)
        outer_mult = l.MultDriveCaller(drive_1=inner_mult, drive_2=self.const_drive_with_idx)
        
        self.assertIsInstance(outer_mult, l.MultDriveCaller)
        self.assertEqual(outer_mult.drive_1, inner_mult)
        self.assertEqual(outer_mult.drive_2, self.const_drive_with_idx)
        
        # Check string representation with nested drive
        expected_str = "mult,\n\tmult,\n\tconst, 3.0,\n\tconst, 4.0,\n\treference, 5"
        self.assertEqual(str(outer_mult), expected_str)

    def test_mult_drive_caller_complex_nesting(self):
        """Test more complex nesting with different drive caller types"""
        # Create an array drive caller
        array_drive = l.ArrayDriveCaller(drives=[self.const_drive1, self.const_drive2])
        
        # Use it in a mult drive caller
        mult_drive = l.MultDriveCaller(drive_1=array_drive, drive_2=self.const_drive_with_idx)
        
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, array_drive)
        self.assertEqual(mult_drive.drive_2, self.const_drive_with_idx)
        
        # Check string representation with complex nesting
        expected_str = "mult,\n\tarray, 2,\n\tconst, 3.0,\n\tconst, 4.0,\n\treference, 5"
        self.assertEqual(str(mult_drive), expected_str)
        
        # Create an even more complex nested structure
        inner_mult = l.MultDriveCaller(drive_1=self.const_drive1, drive_2=self.const_drive2)
        middle_mult = l.MultDriveCaller(drive_1=inner_mult, drive_2=array_drive)
        outer_mult = l.MultDriveCaller(drive_1=middle_mult, drive_2=self.const_drive_with_idx)
        
        self.assertIsInstance(outer_mult, l.MultDriveCaller)
        
        # Check string representation with complex nesting
        expected_str = "mult,\n\tmult,\n\tmult,\n\tconst, 3.0,\n\tconst, 4.0,\n\tarray, 2,\n\tconst, 3.0,\n\tconst, 4.0,\n\treference, 5"
        self.assertEqual(str(outer_mult), expected_str)

class TestNullDriveCaller(unittest.TestCase):
    def test_null_drive_caller_creation_and_representation(self):
        """Test creating a NullDriveCaller and its string representation"""
        # Create without idx
        null_drive = l.NullDriveCaller()
        self.assertIsInstance(null_drive, l.NullDriveCaller)
        self.assertEqual(str(null_drive), 'null')
        self.assertEqual(null_drive.drive_type(), 'null')
        
        # Create with idx
        null_drive = l.NullDriveCaller(idx=5)
        self.assertIsInstance(null_drive, l.NullDriveCaller)
        self.assertEqual(str(null_drive), 'drive caller: 5, null')
        self.assertEqual(null_drive.idx, 5)

class TestParabolicDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create MBVar objects for testing if needed
        if 'const_coef_var' not in l.declared_MBVars:
            self.const_coef_var = l.MBVar(name='const_coef_var', var_type='real', expression=1.5)
        else:
            self.const_coef_var = l.declared_MBVars['const_coef_var']
    
    def test_parabolic_drive_caller_creation_valid(self):
        """Test that ParabolicDriveCaller works with valid numeric inputs"""
        # Create with all parameters
        parabolic_drive = l.ParabolicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        self.assertIsInstance(parabolic_drive, l.ParabolicDriveCaller)
        self.assertEqual(parabolic_drive.const_coef, 1.0)
        self.assertEqual(parabolic_drive.linear_coef, 2.0)
        self.assertEqual(parabolic_drive.parabolic_coef, 3.0)
        
        # Create with idx
        parabolic_drive = l.ParabolicDriveCaller(
            idx=10,
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        self.assertEqual(parabolic_drive.idx, 10)
    
    def test_parabolic_drive_caller_with_mbvars(self):
        """Test ParabolicDriveCaller with MBVar objects"""
        parabolic_drive = l.ParabolicDriveCaller(
            const_coef=self.const_coef_var,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        self.assertEqual(parabolic_drive.const_coef, self.const_coef_var)
    
    def test_parabolic_drive_caller_str_representation(self):
        """Test string representation of ParabolicDriveCaller"""
        # Test without idx
        parabolic_drive = l.ParabolicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        expected_str = "parabolic, 1.0, 2.0, 3.0"
        self.assertEqual(str(parabolic_drive), expected_str)
        
        # Test with idx
        parabolic_drive = l.ParabolicDriveCaller(
            idx=10,
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        expected_str = "drive caller: 10, parabolic, 1.0, 2.0, 3.0"
        self.assertEqual(str(parabolic_drive), expected_str)
        
        # Test with MBVar
        parabolic_drive = l.ParabolicDriveCaller(
            const_coef=self.const_coef_var,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        expected_str = f"parabolic, {self.const_coef_var}, 2.0, 3.0"
        self.assertEqual(str(parabolic_drive), expected_str)
    
    def test_parabolic_drive_caller_drive_type(self):
        """Test the drive_type method"""
        parabolic_drive = l.ParabolicDriveCaller(
            const_coef=1.0,
            linear_coef=2.0,
            parabolic_coef=3.0
        )
        self.assertEqual(parabolic_drive.drive_type(), "parabolic")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_parabolic_drive_caller_missing_required_field(self):
        """Test creating a ParabolicDriveCaller missing a required field"""
        # Missing const_coef
        with self.assertRaises(Exception):
            l.ParabolicDriveCaller(
                linear_coef=2.0,
                parabolic_coef=3.0
            )
        
        # Missing linear_coef
        with self.assertRaises(Exception):
            l.ParabolicDriveCaller(
                const_coef=1.0,
                parabolic_coef=3.0
            )
        
        # Missing parabolic_coef
        with self.assertRaises(Exception):
            l.ParabolicDriveCaller(
                const_coef=1.0,
                linear_coef=2.0
            )

    unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_parabolic_drive_caller_invalid_types(self):
        """Test invalid types for ParabolicDriveCaller fields"""
        # Invalid type for const_coef
        with self.assertRaises(Exception):
            l.ParabolicDriveCaller(
                const_coef="invalid string",
                linear_coef=2.0,
                parabolic_coef=3.0
            )
        
        # Invalid type for linear_coef
        with self.assertRaises(Exception):
            l.ParabolicDriveCaller(
                const_coef=1.0,
                linear_coef=[1, 2, 3],  # List is invalid
                parabolic_coef=3.0
            )
        
        # Invalid MBVar type for parabolic_coef (using string MBVar)
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.ParabolicDriveCaller(
                const_coef=1.0,
                linear_coef=2.0,
                parabolic_coef=string_var  # String MBVar is invalid for numeric field
            )

class TestPeriodicDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create sample drive callers for testing
        self.const_drive = l.ConstDriveCaller(const_value=5.0)
        self.const_drive_with_idx = l.ConstDriveCaller(idx=5, const_value=10.0)
    
    def test_periodic_drive_caller_creation_valid(self):
        """Test that PeriodicDriveCaller works with valid inputs"""
        # Create with all parameters
        periodic_drive = l.PeriodicDriveCaller(
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        self.assertIsInstance(periodic_drive, l.PeriodicDriveCaller)
        self.assertEqual(periodic_drive.initial_time, 1.0)
        self.assertEqual(periodic_drive.period, 2.0)
        self.assertEqual(periodic_drive.func_drive, self.const_drive)
        
        # Create with default initial_time
        periodic_drive = l.PeriodicDriveCaller(
            period=2.0,
            func_drive=self.const_drive
        )
        self.assertEqual(periodic_drive.initial_time, 0.0)
        
        # Create with idx
        periodic_drive = l.PeriodicDriveCaller(
            idx=10,
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        self.assertEqual(periodic_drive.idx, 10)
    
    def test_periodic_drive_caller_default_warning(self):
        """Test warning when initial_time is not specified"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            periodic_drive = l.PeriodicDriveCaller(
                period=2.0,
                func_drive=self.const_drive
            )
            
            # Verify a warning was raised
            self.assertTrue(len(w) > 0)
            self.assertTrue(any("<initial_time> is not set" in str(warning.message) for warning in w))
    
    def test_periodic_drive_caller_str_representation(self):
        """Test string representation of PeriodicDriveCaller"""
        # Test with regular drive
        periodic_drive = l.PeriodicDriveCaller(
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        expected_str = "periodic, 1.0, 2.0,\n\tconst, 5.0"
        self.assertEqual(str(periodic_drive), expected_str)
        
        # Test with idx
        periodic_drive = l.PeriodicDriveCaller(
            idx=10,
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        expected_str = "drive caller: 10, periodic, 1.0, 2.0,\n\tconst, 5.0"
        self.assertEqual(str(periodic_drive), expected_str)
        
        # Test with referenced drive
        periodic_drive = l.PeriodicDriveCaller(
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive_with_idx
        )
        expected_str = "periodic, 1.0, 2.0,\n\treference, 5"
        self.assertEqual(str(periodic_drive), expected_str)
    
    def test_periodic_drive_caller_drive_type(self):
        """Test the drive_type method"""
        periodic_drive = l.PeriodicDriveCaller(
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        self.assertEqual(periodic_drive.drive_type(), "periodic")
    
    def test_periodic_drive_caller_nested(self):
        """Test nesting PeriodicDriveCaller with other drives"""
        # Create a periodic drive
        periodic_drive = l.PeriodicDriveCaller(
            initial_time=1.0,
            period=2.0,
            func_drive=self.const_drive
        )
        
        # Create a drive that uses the periodic drive
        outer_drive = l.DriveDriveCaller(
            drive_caller1=periodic_drive,
            drive_caller2=self.const_drive_with_idx
        )
        
        self.assertIsInstance(outer_drive, l.DriveDriveCaller)
        self.assertEqual(outer_drive.drive_caller1, periodic_drive)
        
        # Check string representation
        expected_str = "drive,\n\tperiodic, 1.0, 2.0,\n\tconst, 5.0,\n\treference, 5"
        self.assertEqual(str(outer_drive), expected_str)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_periodic_drive_caller_missing_required_field(self):
        """Test creating a PeriodicDriveCaller missing a required field"""
        # Missing period
        with self.assertRaises(Exception):
            l.PeriodicDriveCaller(
                initial_time=1.0,
                func_drive=self.const_drive
            )
        
        # Missing func_drive
        with self.assertRaises(Exception):
            l.PeriodicDriveCaller(
                initial_time=1.0,
                period=2.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_periodic_drive_caller_invalid_types(self):
        """Test invalid types for PeriodicDriveCaller fields"""
        # Invalid type for initial_time
        with self.assertRaises(Exception):
            l.PeriodicDriveCaller(
                initial_time="invalid string",
                period=2.0,
                func_drive=self.const_drive
            )
        
        # Invalid type for period
        with self.assertRaises(Exception):
            l.PeriodicDriveCaller(
                initial_time=1.0,
                period=[1, 2, 3],  # List is invalid
                func_drive=self.const_drive
            )
        
        # Invalid type for func_drive (not a DriveCaller)
        with self.assertRaises(Exception):
            l.PeriodicDriveCaller(
                initial_time=1.0,
                period=2.0,
                func_drive="not a drive caller"
            )
        
        # Invalid MBVar type for initial_time (using string MBVar)
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.PeriodicDriveCaller(
                initial_time=string_var,  # String MBVar is invalid for numeric field
                period=2.0,
                func_drive=self.const_drive
            )

class TestPiecewiseLinearDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Sample points and values for testing
        self.num_points = 3
        self.points_values = [(0.0, 0.0), (1.0, 2.0), (2.0, 1.0)]
        
        # Create MBVar objects for testing if needed
        if 'points_var' not in l.declared_MBVars:
            self.points_var = l.MBVar(name='points_var', var_type='integer', expression=3)
        else:
            self.points_var = l.declared_MBVars['points_var']
    
    def test_piecewise_linear_drive_caller_creation_valid(self):
        """Test creating a PiecewiseLinearDriveCaller with valid parameters"""
        # Create with integer num_points
        drive = l.PiecewiseLinearDriveCaller(
            num_points=self.num_points,
            points_values=self.points_values
        )
        self.assertIsInstance(drive, l.PiecewiseLinearDriveCaller)
        self.assertEqual(drive.num_points, self.num_points)
        self.assertEqual(drive.points_values, self.points_values)
        self.assertEqual(drive.drive_type(), "piecewise linear")
        
        # Create with MBVar for num_points
        drive = l.PiecewiseLinearDriveCaller(
            num_points=self.points_var,
            points_values=self.points_values
        )
        self.assertIsInstance(drive, l.PiecewiseLinearDriveCaller)
        self.assertEqual(drive.num_points, self.points_var)
        
        # Create with idx
        drive = l.PiecewiseLinearDriveCaller(
            idx=10,
            num_points=self.num_points,
            points_values=self.points_values
        )
        self.assertEqual(drive.idx, 10)
    
    def test_piecewise_linear_drive_caller_with_mbvars_in_points(self):
        """Test using MBVar objects in the points_values list"""
        # Create MBVars for points and values
        if 'point1' not in l.declared_MBVars:
            point1_var = l.MBVar(name='point1', var_type='real', expression=1.5)
        else:
            point1_var = l.declared_MBVars['point1']
            
        if 'value1' not in l.declared_MBVars:
            value1_var = l.MBVar(name='value1', var_type='real', expression=3.0)
        else:
            value1_var = l.declared_MBVars['value1']
        
        # Create points_values with MBVars
        points_values_with_mbvars = [(0.0, 0.0), (point1_var, value1_var), (2.0, 1.0)]
        
        # Create drive with MBVars in points_values
        drive = l.PiecewiseLinearDriveCaller(
            num_points=self.num_points,
            points_values=points_values_with_mbvars
        )
        self.assertIsInstance(drive, l.PiecewiseLinearDriveCaller)
        self.assertEqual(drive.points_values[1][0], point1_var)
        self.assertEqual(drive.points_values[1][1], value1_var)
    
    def test_piecewise_linear_drive_caller_str_representation(self):
        """Test string representation of the drive caller"""
        # Create a drive
        drive = l.PiecewiseLinearDriveCaller(
            num_points=self.num_points,
            points_values=self.points_values
        )
        
        # Expected string representation
        expected_str = "piecewise linear, 3,\n\t0.0, 0.0,\n\t1.0, 2.0,\n\t2.0, 1.0"
        
        # Test string representation
        self.assertEqual(str(drive), expected_str)
        
        # Create with idx
        drive = l.PiecewiseLinearDriveCaller(
            idx=10,
            num_points=self.num_points,
            points_values=self.points_values
        )
        
        # Expected string with idx
        expected_str_with_idx = "drive caller: 10, piecewise linear, 3,\n\t0.0, 0.0,\n\t1.0, 2.0,\n\t2.0, 1.0"
        
        self.assertEqual(str(drive), expected_str_with_idx)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_piecewise_linear_drive_caller_validation(self):
        """Test validation of num_points against length of points_values"""
        # Test with mismatched num_points and points_values length
        with self.assertRaises(ValueError) as context:
            l.PiecewiseLinearDriveCaller(
                num_points=4,  # Doesn't match the 3 points below
                points_values=self.points_values
            )
        self.assertIn("number of (point, value) pairs", str(context.exception))
        
        # Test with too few points
        with self.assertRaises(ValueError):
            l.PiecewiseLinearDriveCaller(
                num_points=1,
                points_values=[(0.0, 0.0)]  # Need at least 2 points for interpolation
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_piecewise_linear_drive_caller_missing_required_field(self):
        """Test validation when required fields are missing"""
        # Missing num_points
        with self.assertRaises(Exception):
            l.PiecewiseLinearDriveCaller(
                points_values=self.points_values
            )
        
        # Missing points_values
        with self.assertRaises(Exception):
            l.PiecewiseLinearDriveCaller(
                num_points=self.num_points
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_piecewise_linear_drive_caller_invalid_types(self):
        """Test with invalid field types"""
        # Invalid type for num_points
        with self.assertRaises(Exception):
            l.PiecewiseLinearDriveCaller(
                num_points="not an integer",
                points_values=self.points_values
            )
        
        # Invalid type for points_values
        with self.assertRaises(Exception):
            l.PiecewiseLinearDriveCaller(
                num_points=self.num_points,
                points_values="not a list"
            )
        
        # # Works with lists too
        # # Invalid structure in points_values (not tuples)
        # with self.assertRaises(Exception):
        #     l.PiecewiseLinearDriveCaller(
        #         num_points=3,
        #         points_values=[[0.0, 0.0], [1.0, 2.0], [2.0, 1.0]]  # Lists instead of tuples
        #     )
        
        # Invalid values in tuples
        with self.assertRaises(Exception):
            l.PiecewiseLinearDriveCaller(
                num_points=3,
                points_values=[(0.0, 0.0), ("string", 2.0), (2.0, 1.0)]  # String instead of number
            )

class TestPostponedDriveCaller(unittest.TestCase):
    def test_postponed_drive_caller_creation_valid(self):
        """Test creating a PostponedDriveCaller with valid parameters"""
        # Create without idx but with required label
        drive = l.PostponedDriveCaller(label=42)
        self.assertIsInstance(drive, l.PostponedDriveCaller)
        self.assertEqual(drive.label, 42)
        
        # Create with idx and label
        drive = l.PostponedDriveCaller(idx=10, label=42)
        self.assertIsInstance(drive, l.PostponedDriveCaller)
        self.assertEqual(drive.idx, 10)
        self.assertEqual(drive.label, 42)
    
    def test_postponed_drive_caller_with_mbvar_label(self):
        """Test using MBVar as label"""
        if 'post_label' not in l.declared_MBVars:
            label_var = l.MBVar(name='post_label', var_type='integer', expression=5)
        else:
            label_var = l.declared_MBVars['post_label']
        
        drive = l.PostponedDriveCaller(label=label_var)
        self.assertIsInstance(drive, l.PostponedDriveCaller)
        self.assertEqual(drive.label, label_var)
    
    def test_postponed_drive_caller_str_representation(self):
        """Test string representation of the drive caller"""
        # Create drive without idx
        drive = l.PostponedDriveCaller(label=42)
        expected_str = "postponed, 42"
        self.assertEqual(str(drive), expected_str)
        
        # Create drive with idx
        drive = l.PostponedDriveCaller(idx=10, label=42)
        expected_str = "drive caller: 10, postponed, 42"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVar label
        if 'post_label' not in l.declared_MBVars:
            label_var = l.MBVar(name='post_label', var_type='integer', expression=5)
        else:
            label_var = l.declared_MBVars['post_label']
        
        drive = l.PostponedDriveCaller(label=label_var)
        expected_str = "postponed, post_label"
        self.assertEqual(str(drive), expected_str)
    
    def test_postponed_drive_caller_drive_type(self):
        """Test the drive_type method"""
        drive = l.PostponedDriveCaller(label=42)
        self.assertEqual(drive.drive_type(), "postponed")
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_postponed_drive_caller_missing_required_field(self):
        """Test creating a PostponedDriveCaller missing the required label field"""
        with self.assertRaises(Exception):
            l.PostponedDriveCaller()  # Missing required label
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_postponed_drive_caller_invalid_types(self):
        """Test with invalid field types"""
        # Invalid type for idx
        with self.assertRaises(Exception):
            l.PostponedDriveCaller(idx="not an integer", label=42)
            
        # Invalid type for label (should be int or MBVar)
        with self.assertRaises(Exception):
            l.PostponedDriveCaller(label="not an integer")

class TestRampDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create MBVar objects for testing
        if 'slope_var' not in l.declared_MBVars:
            self.slope_var = l.MBVar(name='slope_var', var_type='real', expression=2.5)
        else:
            self.slope_var = l.declared_MBVars['slope_var']
    
    def test_ramp_drive_caller_creation_valid(self):
        """Test that RampDriveCaller works with valid inputs"""
        # Create with all parameters
        ramp_drive = l.RampDriveCaller(
            slope=1.0,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        self.assertIsInstance(ramp_drive, l.RampDriveCaller)
        self.assertEqual(ramp_drive.slope, 1.0)
        self.assertEqual(ramp_drive.initial_time, 0.0)
        self.assertEqual(ramp_drive.final_time, 10.0)
        self.assertEqual(ramp_drive.initial_value, 5.0)
        
        # Create with 'forever' as final_time
        ramp_drive = l.RampDriveCaller(
            slope=1.0,
            initial_time=0.0,
            final_time='forever',
            initial_value=5.0
        )
        self.assertIsInstance(ramp_drive, l.RampDriveCaller)
        self.assertEqual(ramp_drive.final_time, 'forever')
        
        # Create with idx
        ramp_drive = l.RampDriveCaller(
            idx=10,
            slope=1.0,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        self.assertEqual(ramp_drive.idx, 10)
    
    def test_ramp_drive_caller_with_mbvars(self):
        """Test RampDriveCaller with MBVar objects"""
        ramp_drive = l.RampDriveCaller(
            slope=self.slope_var,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        self.assertEqual(ramp_drive.slope, self.slope_var)
    
    def test_ramp_drive_caller_default_warning(self):
        """Test warnings for default parameters"""
        # Test warning for default initial_time
        with warnings.catch_warnings(record=True) as w:
            ramp_drive = l.RampDriveCaller(
                slope=1.0,
                final_time=10.0
            )
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(ramp_drive.initial_time, 0.0)
        
        # Test warning for default initial_value
        with warnings.catch_warnings(record=True) as w:
            ramp_drive = l.RampDriveCaller(
                slope=1.0,
                final_time=10.0,
                initial_time=0.0
            )
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(ramp_drive.initial_value, 0.0)
    
    def test_ramp_drive_caller_str_representation(self):
        """Test string representation of RampDriveCaller"""
        # Test without idx
        ramp_drive = l.RampDriveCaller(
            slope=1.0,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        expected_str = "ramp, 1.0, 0.0, 10.0, 5.0"
        self.assertEqual(str(ramp_drive), expected_str)
        
        # Test with idx
        ramp_drive = l.RampDriveCaller(
            idx=10,
            slope=1.0,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        expected_str = "drive caller: 10, ramp, 1.0, 0.0, 10.0, 5.0"
        self.assertEqual(str(ramp_drive), expected_str)
        
        # Test with 'forever' and MBVar
        ramp_drive = l.RampDriveCaller(
            slope=self.slope_var,
            initial_time=0.0,
            final_time='forever',
            initial_value=5.0
        )
        expected_str = f"ramp, {self.slope_var}, 0.0, forever, 5.0"
        self.assertEqual(str(ramp_drive), expected_str)
    
    def test_ramp_drive_caller_drive_type(self):
        """Test the drive_type method"""
        ramp_drive = l.RampDriveCaller(
            slope=1.0,
            initial_time=0.0,
            final_time=10.0,
            initial_value=5.0
        )
        self.assertEqual(ramp_drive.drive_type(), "ramp")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_ramp_drive_caller_missing_required_field(self):
        """Test creating a RampDriveCaller missing a required field"""
        # Missing slope
        with self.assertRaises(Exception):
            l.RampDriveCaller(
                initial_time=0.0,
                final_time=10.0,
                initial_value=5.0
            )
        
        # Missing final_time
        with self.assertRaises(Exception):
            l.RampDriveCaller(
                slope=1.0,
                initial_time=0.0,
                initial_value=5.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_ramp_drive_caller_invalid_types(self):
        """Test invalid types for RampDriveCaller fields"""
        # Invalid type for slope
        with self.assertRaises(Exception):
            l.RampDriveCaller(
                slope="invalid string",
                initial_time=0.0,
                final_time=10.0,
                initial_value=5.0
            )
        
        # Invalid MBVar type for initial_value (using string MBVar)
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.RampDriveCaller(
                slope=1.0,
                initial_time=0.0,
                final_time=10.0,
                initial_value=string_var  # String MBVar is invalid for numeric field
            )

class TestRandomDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create MBVar objects for testing
        if 'amplitude_var' not in l.declared_MBVars:
            self.amplitude_var = l.MBVar(name='amplitude_var', var_type='real', expression=2.5)
        else:
            self.amplitude_var = l.declared_MBVars['amplitude_var']
            
        if 'steps_var' not in l.declared_MBVars:
            self.steps_var = l.MBVar(name='steps_var', var_type='integer', expression=5)
        else:
            self.steps_var = l.declared_MBVars['steps_var']
    
    def test_random_drive_caller_creation_valid(self):
        """Test that RandomDriveCaller works with valid inputs"""
        # Create with required parameters
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0
        )
        self.assertIsInstance(random_drive, l.RandomDriveCaller)
        self.assertEqual(random_drive.amplitude_value, 1.0)
        self.assertEqual(random_drive.mean_value, 0.0)
        self.assertEqual(random_drive.initial_time, 0.0)
        self.assertEqual(random_drive.final_time, 10.0)
        self.assertIsNone(random_drive.steps_to_hold_value)
        self.assertIsNone(random_drive.seed_value)
        
        # Create with all parameters
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0,
            steps_to_hold_value=3,
            seed_value=42
        )
        self.assertIsInstance(random_drive, l.RandomDriveCaller)
        self.assertEqual(random_drive.steps_to_hold_value, 3)
        self.assertEqual(random_drive.seed_value, 42)
        
        # Create with 'time' as seed_value
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0,
            seed_value='time'
        )
        self.assertEqual(random_drive.seed_value, 'time')
        
        # Create with 'forever' as final_time
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time='forever'
        )
        self.assertEqual(random_drive.final_time, 'forever')
    
    def test_random_drive_caller_with_mbvars(self):
        """Test RandomDriveCaller with MBVar objects"""
        random_drive = l.RandomDriveCaller(
            amplitude_value=self.amplitude_var,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0,
            steps_to_hold_value=self.steps_var
        )
        self.assertEqual(random_drive.amplitude_value, self.amplitude_var)
        self.assertEqual(random_drive.steps_to_hold_value, self.steps_var)
    
    def test_random_drive_caller_default_warning(self):
        """Test warning for default initial_time"""
        with warnings.catch_warnings(record=True) as w:
            random_drive = l.RandomDriveCaller(
                amplitude_value=1.0,
                mean_value=0.0,
                final_time=10.0
            )
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(random_drive.initial_time, 0.0)
    
    def test_random_drive_caller_str_representation(self):
        """Test string representation of RandomDriveCaller"""
        # Test basic parameters
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0
        )
        expected_str = "random, 1.0, 0.0, 0.0, 10.0"
        self.assertEqual(str(random_drive), expected_str)
        
        # Test with all parameters
        random_drive = l.RandomDriveCaller(
            idx=5,
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0,
            steps_to_hold_value=3,
            seed_value=42
        )
        expected_str = "drive caller: 5, random, 1.0, 0.0, 0.0, 10.0, steps, 3, seed, 42"
        self.assertEqual(str(random_drive), expected_str)
        
        # Test with 'time' and 'forever'
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time='forever',
            seed_value='time'
        )
        expected_str = "random, 1.0, 0.0, 0.0, forever, seed, time"
        self.assertEqual(str(random_drive), expected_str)
    
    def test_random_drive_caller_drive_type(self):
        """Test the drive_type method"""
        random_drive = l.RandomDriveCaller(
            amplitude_value=1.0,
            mean_value=0.0,
            initial_time=0.0,
            final_time=10.0
        )
        self.assertEqual(random_drive.drive_type(), "random")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_random_drive_caller_missing_required_field(self):
        """Test creating a RandomDriveCaller missing a required field"""
        # Missing amplitude_value
        with self.assertRaises(Exception):
            l.RandomDriveCaller(
                mean_value=0.0,
                initial_time=0.0,
                final_time=10.0
            )
        
        # Missing mean_value
        with self.assertRaises(Exception):
            l.RandomDriveCaller(
                amplitude_value=1.0,
                initial_time=0.0,
                final_time=10.0
            )
        
        # Missing final_time
        with self.assertRaises(Exception):
            l.RandomDriveCaller(
                amplitude_value=1.0,
                mean_value=0.0,
                initial_time=0.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_random_drive_caller_invalid_types(self):
        """Test invalid types for RandomDriveCaller fields"""
        # Invalid type for steps_to_hold_value
        with self.assertRaises(Exception):
            l.RandomDriveCaller(
                amplitude_value=1.0,
                mean_value=0.0,
                initial_time=0.0,
                final_time=10.0,
                steps_to_hold_value="not an integer"
            )
        
        # Invalid value for steps_to_hold_value (must be positive)
        with self.assertRaises(ValueError):
            l.RandomDriveCaller(
                amplitude_value=1.0,
                mean_value=0.0,
                initial_time=0.0,
                final_time=10.0,
                steps_to_hold_value=0
            )
        
        # Invalid type for seed_value
        with self.assertRaises(Exception):
            l.RandomDriveCaller(
                amplitude_value=1.0,
                mean_value=0.0,
                initial_time=0.0,
                final_time=10.0,
                seed_value=3.14  # Should be integer or 'time'
            )


class TestSampleAndHoldDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create sample drive callers for testing
        self.function_drive = l.ConstDriveCaller(const_value=5.0)
        self.trigger_drive = l.ConstDriveCaller(const_value=1.0)
        self.function_drive_with_idx = l.ConstDriveCaller(idx=10, const_value=5.0)
        self.trigger_drive_with_idx = l.ConstDriveCaller(idx=20, const_value=1.0)
        
        # Create MBVar for testing
        if 'initial_val_var' not in l.declared_MBVars:
            self.initial_val_var = l.MBVar(name='initial_val_var', var_type='real', expression=3.0)
        else:
            self.initial_val_var = l.declared_MBVars['initial_val_var']
    
    def test_sample_and_hold_drive_caller_creation_valid(self):
        """Test creating a SampleAndHoldDriveCaller with valid parameters"""
        # Create with required parameters
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive
        )
        self.assertIsInstance(drive, l.SampleAndHoldDriveCaller)
        self.assertEqual(drive.function, self.function_drive)
        self.assertEqual(drive.trigger, self.trigger_drive)
        self.assertIsNone(drive.initial_value)
        
        # Create with initial_value
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive,
            initial_value=2.0
        )
        self.assertIsInstance(drive, l.SampleAndHoldDriveCaller)
        self.assertEqual(drive.initial_value, 2.0)
        
        # Create with idx
        drive = l.SampleAndHoldDriveCaller(
            idx=5,
            function=self.function_drive,
            trigger=self.trigger_drive
        )
        self.assertEqual(drive.idx, 5)
    
    def test_sample_and_hold_drive_caller_with_mbvars(self):
        """Test SampleAndHoldDriveCaller with MBVar for initial_value"""
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive,
            initial_value=self.initial_val_var
        )
        self.assertEqual(drive.initial_value, self.initial_val_var)
    
    def test_sample_and_hold_drive_caller_with_reference_drives(self):
        """Test SampleAndHoldDriveCaller with drives that have idx (reference)"""
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive_with_idx,
            trigger=self.trigger_drive_with_idx
        )
        self.assertEqual(drive.function, self.function_drive_with_idx)
        self.assertEqual(drive.trigger, self.trigger_drive_with_idx)
    
    def test_sample_and_hold_drive_caller_str_representation(self):
        """Test string representation of SampleAndHoldDriveCaller"""
        # Test with regular drives
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive
        )
        expected_str = "sample and hold,\n\tconst, 5.0,\n\tconst, 1.0"
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.SampleAndHoldDriveCaller(
            idx=5,
            function=self.function_drive,
            trigger=self.trigger_drive
        )
        expected_str = "drive caller: 5, sample and hold,\n\tconst, 5.0,\n\tconst, 1.0"
        self.assertEqual(str(drive), expected_str)
        
        # Test with reference drives and initial_value
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive_with_idx,
            trigger=self.trigger_drive_with_idx,
            initial_value=2.0
        )
        expected_str = "sample and hold,\n\treference, 10,\n\treference, 20, initial value, 2.0"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVar initial_value
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive,
            initial_value=self.initial_val_var
        )
        expected_str = f"sample and hold,\n\tconst, 5.0,\n\tconst, 1.0, initial value, {self.initial_val_var}"
        self.assertEqual(str(drive), expected_str)
    
    def test_sample_and_hold_drive_caller_drive_type(self):
        """Test the drive_type method"""
        drive = l.SampleAndHoldDriveCaller(
            function=self.function_drive,
            trigger=self.trigger_drive
        )
        self.assertEqual(drive.drive_type(), "sample and hold")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_sample_and_hold_drive_caller_missing_required_field(self):
        """Test creating a SampleAndHoldDriveCaller missing a required field"""
        # Missing function
        with self.assertRaises(Exception):
            l.SampleAndHoldDriveCaller(
                trigger=self.trigger_drive
            )
        
        # Missing trigger
        with self.assertRaises(Exception):
            l.SampleAndHoldDriveCaller(
                function=self.function_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_sample_and_hold_drive_caller_invalid_types(self):
        """Test invalid types for SampleAndHoldDriveCaller fields"""
        # Invalid type for function
        with self.assertRaises(Exception):
            l.SampleAndHoldDriveCaller(
                function="not a drive caller",
                trigger=self.trigger_drive
            )
        
        # Invalid type for trigger
        with self.assertRaises(Exception):
            l.SampleAndHoldDriveCaller(
                function=self.function_drive,
                trigger="not a drive caller"
            )
        
        # Invalid MBVar type for initial_value
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.SampleAndHoldDriveCaller(
                function=self.function_drive,
                trigger=self.trigger_drive,
                initial_value=string_var  # String MBVar is invalid for numeric field
            )
            
class TestSineDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset warnings to make sure we capture them in tests
        warnings.resetwarnings()
        # Setup common values for testing
        self.initial_time = 1.0
        self.angular_velocity = 2.0
        self.amplitude = 3.0
        self.number_of_cycles = 5
        self.initial_value = 0.5
    
    def test_sine_drive_caller_creation_valid(self):
        """Test creating a SineDriveCaller with valid parameters"""
        # Create with all parameters
        drive = l.SineDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        
        # Verify all properties are set correctly
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.angular_velocity, self.angular_velocity)
        self.assertEqual(drive.amplitude, self.amplitude)
        self.assertEqual(drive.number_of_cycles, self.number_of_cycles)
        self.assertEqual(drive.initial_value, self.initial_value)
        self.assertEqual(drive.drive_type(), "sine")
        
        # Test with string values for number_of_cycles
        for cycle_value in ['half', 'one', 'forever']:
            drive = l.SineDriveCaller(
                initial_time=self.initial_time,
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles=cycle_value,
                initial_value=self.initial_value
            )
            self.assertEqual(drive.number_of_cycles, cycle_value)
        
        # Test with idx parameter
        drive = l.SineDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        self.assertEqual(drive.idx, 5)
    
    def test_sine_drive_caller_default_values(self):
        """Test that default values are set correctly with warnings"""
        with warnings.catch_warnings(record=True) as w:
            # Create drive without initial_time and initial_value
            drive = l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
            
            # Check default values
            self.assertEqual(drive.initial_time, 0.0)
            self.assertEqual(drive.initial_value, 0.0)
            
            # Verify warnings were raised
            self.assertEqual(len(w), 2)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertTrue("<initial_time> is not set, assuming 0.0." in str(w[0].message))
            self.assertTrue(issubclass(w[1].category, UserWarning))
            self.assertTrue("<initial_value> is not set, assuming 0.0." in str(w[1].message))
    
    def test_sine_drive_caller_with_mbvars(self):
        """Test creating a SineDriveCaller with MBVar objects"""
        # Create MBVar objects
        if 'init_time' not in l.declared_MBVars:
            initial_time_var = l.MBVar(name='init_time', var_type='real', expression=1.5)
        else:
            initial_time_var = l.declared_MBVars['init_time']
            
        if 'omega' not in l.declared_MBVars:
            angular_velocity_var = l.MBVar(name='omega', var_type='real', expression=3.5)
        else:
            angular_velocity_var = l.declared_MBVars['omega']
            
        if 'amp' not in l.declared_MBVars:
            amplitude_var = l.MBVar(name='amp', var_type='real', expression=4.0)
        else:
            amplitude_var = l.declared_MBVars['amp']
            
        if 'cycles' not in l.declared_MBVars:
            number_of_cycles_var = l.MBVar(name='cycles', var_type='integer', expression=10)
        else:
            number_of_cycles_var = l.declared_MBVars['cycles']
            
        if 'init_val' not in l.declared_MBVars:
            initial_value_var = l.MBVar(name='init_val', var_type='real', expression=0.75)
        else:
            initial_value_var = l.declared_MBVars['init_val']
        
        # Create with all MBVar objects
        drive = l.SineDriveCaller(
            initial_time=initial_time_var,
            angular_velocity=angular_velocity_var,
            amplitude=amplitude_var,
            number_of_cycles=number_of_cycles_var,
            initial_value=initial_value_var
        )
        
        # Check that MBVar references are stored correctly
        self.assertEqual(drive.initial_time, initial_time_var)
        self.assertEqual(drive.angular_velocity, angular_velocity_var)
        self.assertEqual(drive.amplitude, amplitude_var)
        self.assertEqual(drive.number_of_cycles, number_of_cycles_var)
        self.assertEqual(drive.initial_value, initial_value_var)
        
        # Create with mixed parameters (some MBVars, some values)
        drive = l.SineDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=angular_velocity_var,
            amplitude=self.amplitude,
            number_of_cycles='forever',
            initial_value=initial_value_var
        )
        
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.angular_velocity, angular_velocity_var)
        self.assertEqual(drive.amplitude, self.amplitude)
        self.assertEqual(drive.number_of_cycles, 'forever')
        self.assertEqual(drive.initial_value, initial_value_var)
    
    def test_sine_drive_caller_str_representation(self):
        """Test the string representation of SineDriveCaller"""
        # Test without idx
        drive = l.SineDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        
        expected_str = f"sine, {self.initial_time}, {self.angular_velocity}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.SineDriveCaller(
            idx=5,
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        
        expected_str = f"drive caller: 5, sine, {self.initial_time}, {self.angular_velocity}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with string number_of_cycles
        drive = l.SineDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles='forever',
            initial_value=self.initial_value
        )
        
        expected_str = f"sine, {self.initial_time}, {self.angular_velocity}, {self.amplitude}, forever, {self.initial_value}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with MBVars
        initial_time_var = l.MBVar(name='init_time', var_type='real', expression=1.5)
        angular_velocity_var = l.MBVar(name='omega', var_type='real', expression=3.5)
        
        drive = l.SineDriveCaller(
            initial_time=initial_time_var,
            angular_velocity=angular_velocity_var,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        
        expected_str = f"sine, {initial_time_var}, {angular_velocity_var}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value}"
        self.assertEqual(str(drive), expected_str)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_sine_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing angular_velocity
        with self.assertRaises(Exception):
            l.SineDriveCaller(
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        
        # Missing amplitude
        with self.assertRaises(Exception):
            l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                number_of_cycles=self.number_of_cycles
            )
        
        # Missing number_of_cycles
        with self.assertRaises(Exception):
            l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_sine_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # Invalid type for angular_velocity
        with self.assertRaises(Exception):
            l.SineDriveCaller(
                angular_velocity="invalid",
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        
        # Invalid type for amplitude
        with self.assertRaises(Exception):
            l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude="invalid",
                number_of_cycles=self.number_of_cycles
            )
        
        # Invalid string for number_of_cycles
        with self.assertRaises(ValueError):
            l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles="invalid"
            )
        
        # Create an MBVar with a non-real type
        if 'test_non_real' not in l.declared_MBVars:
            non_real_var = l.MBVar(name='test_non_real', var_type='integer', expression=5)
        else:
            non_real_var = l.declared_MBVars['test_non_real']
        
        # Create an MBVar with a non-integer type
        if 'test_non_integer' not in l.declared_MBVars:
            non_integer_var = l.MBVar(name='test_non_integer', var_type='real', expression=5.0)
        else:
            non_integer_var = l.declared_MBVars['test_non_integer']
        
        # Test with non-real MBVar for initial_time
        with self.assertRaises(TypeError) as context:
            l.SineDriveCaller(
                initial_time=non_real_var,
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        self.assertIn("Field must be an MBVar of type real or a float", str(context.exception))
        
        # Test with non-real MBVar for angular_velocity
        with self.assertRaises(TypeError) as context:
            l.SineDriveCaller(
                angular_velocity=non_real_var,
                amplitude=self.amplitude,
                number_of_cycles=self.number_of_cycles
            )
        self.assertIn("Field must be an MBVar of type real or a float", str(context.exception))
        
        # Test with non-integer MBVar for number_of_cycles
        with self.assertRaises(TypeError) as context:
            l.SineDriveCaller(
                angular_velocity=self.angular_velocity,
                amplitude=self.amplitude,
                number_of_cycles=non_integer_var
            )
        self.assertIn("number_of_cycles must be an integer", str(context.exception))
    
    def test_sine_drive_caller_in_other_drives(self):
        """Test using SineDriveCaller as input to other drive callers"""
        # Create a SineDriveCaller
        sine_drive = l.SineDriveCaller(
            initial_time=self.initial_time,
            angular_velocity=self.angular_velocity,
            amplitude=self.amplitude,
            number_of_cycles=self.number_of_cycles,
            initial_value=self.initial_value
        )
        
        # Use it in a MultDriveCaller
        const_drive = l.ConstDriveCaller(const_value=2.0)
        mult_drive = l.MultDriveCaller(
            drive_1=sine_drive,
            drive_2=const_drive
        )
        
        self.assertIsInstance(mult_drive, l.MultDriveCaller)
        self.assertEqual(mult_drive.drive_1, sine_drive)
        self.assertEqual(mult_drive.drive_2, const_drive)
        
        # Check string representation
        expected_mult_str = f"mult,\n\tsine, {self.initial_time}, {self.angular_velocity}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value},\n\tconst, 2.0"
        self.assertEqual(str(mult_drive), expected_mult_str)
        
        # Use it in an ArrayDriveCaller
        array_drive = l.ArrayDriveCaller(drives=[sine_drive, const_drive])
        
        self.assertIsInstance(array_drive, l.ArrayDriveCaller)
        self.assertEqual(len(array_drive.drives), 2)
        self.assertEqual(array_drive.drives[0], sine_drive)
        
        # Check string representation
        expected_array_str = f"array, 2,\n\tsine, {self.initial_time}, {self.angular_velocity}, {self.amplitude}, {self.number_of_cycles}, {self.initial_value},\n\tconst, 2.0"
        self.assertEqual(str(array_drive), expected_array_str)

class TestStepDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create MBVar objects for testing if needed
        if 'step_var' not in l.declared_MBVars:
            self.step_var = l.MBVar(name='step_var', var_type='real', expression=2.5)
        else:
            self.step_var = l.declared_MBVars['step_var']
    
    def test_step_drive_caller_creation_valid(self):
        """Test that StepDriveCaller works with valid inputs"""
        # Create with all parameters
        step_drive = l.StepDriveCaller(
            initial_time=1.0,
            step_value=5.0,
            initial_value=0.0
        )
        self.assertIsInstance(step_drive, l.StepDriveCaller)
        self.assertEqual(step_drive.initial_time, 1.0)
        self.assertEqual(step_drive.step_value, 5.0)
        self.assertEqual(step_drive.initial_value, 0.0)
        
        # Create with idx
        step_drive = l.StepDriveCaller(
            idx=10,
            initial_time=1.0,
            step_value=5.0,
            initial_value=0.0
        )
        self.assertEqual(step_drive.idx, 10)
    
    def test_step_drive_caller_with_mbvars(self):
        """Test StepDriveCaller with MBVar objects"""
        step_drive = l.StepDriveCaller(
            initial_time=1.0,
            step_value=self.step_var,
            initial_value=0.0
        )
        self.assertEqual(step_drive.step_value, self.step_var)
    
    def test_step_drive_caller_default_warning(self):
        """Test warnings for default parameters"""
        # Test warning for default initial_time
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            step_drive = l.StepDriveCaller(
                step_value=5.0
            )
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(step_drive.initial_time, 0.0)
        
        # Test warning for default initial_value
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            step_drive = l.StepDriveCaller(
                step_value=5.0,
                initial_time=1.0
            )
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(step_drive.initial_value, 0.0)
    
    def test_step_drive_caller_str_representation(self):
        """Test string representation of StepDriveCaller"""
        # Test without idx
        step_drive = l.StepDriveCaller(
            initial_time=1.0,
            step_value=5.0,
            initial_value=0.0
        )
        expected_str = "step, 1.0, 5.0, 0.0"
        self.assertEqual(str(step_drive), expected_str)
        
        # Test with idx
        step_drive = l.StepDriveCaller(
            idx=10,
            initial_time=1.0,
            step_value=5.0,
            initial_value=0.0
        )
        expected_str = "drive caller: 10, step, 1.0, 5.0, 0.0"
        self.assertEqual(str(step_drive), expected_str)
        
        # Test with MBVar
        step_drive = l.StepDriveCaller(
            initial_time=1.0,
            step_value=self.step_var,
            initial_value=0.0
        )
        expected_str = f"step, 1.0, {self.step_var}, 0.0"
        self.assertEqual(str(step_drive), expected_str)
    
    def test_step_drive_caller_drive_type(self):
        """Test the drive_type method"""
        step_drive = l.StepDriveCaller(
            initial_time=1.0,
            step_value=5.0,
            initial_value=0.0
        )
        self.assertEqual(step_drive.drive_type(), "step")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_step_drive_caller_missing_required_field(self):
        """Test creating a StepDriveCaller missing a required field"""
        # Missing step_value
        with self.assertRaises(Exception):
            l.StepDriveCaller(
                initial_time=1.0,
                initial_value=0.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_step_drive_caller_invalid_types(self):
        """Test invalid types for StepDriveCaller fields"""
        # Invalid type for step_value
        with self.assertRaises(Exception):
            l.StepDriveCaller(
                initial_time=1.0,
                step_value="invalid string",
                initial_value=0.0
            )
        
        # Invalid MBVar type (using string MBVar)
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.StepDriveCaller(
                initial_time=1.0,
                step_value=string_var,  # String MBVar is invalid for numeric field
                initial_value=0.0
            )

class TestStep5DriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset warnings to make sure we capture them in tests
        warnings.resetwarnings()
        # Setup common values for testing
        self.initial_time = 1.0
        self.initial_value = 0.5
        self.final_time = 5.0
        self.final_value = 2.0
    
    def test_step5_drive_caller_creation_valid(self):
        """Test creating a Step5DriveCaller with valid parameters"""
        # Create with all parameters
        drive = l.Step5DriveCaller(
            initial_time=self.initial_time,
            initial_value=self.initial_value,
            final_time=self.final_time,
            final_value=self.final_value
        )
        
        # Verify all properties are set correctly
        self.assertEqual(drive.initial_time, self.initial_time)
        self.assertEqual(drive.initial_value, self.initial_value)
        self.assertEqual(drive.final_time, self.final_time)
        self.assertEqual(drive.final_value, self.final_value)
        self.assertEqual(drive.drive_type(), "step5")
    
    def test_step5_drive_caller_default_values(self):
        """Test that default values are set correctly with warnings"""
        with warnings.catch_warnings(record=True) as w:
            # Create drive without initial_time and initial_value
            drive = l.Step5DriveCaller(
                final_time=self.final_time,
                final_value=self.final_value
            )
            
            # Check default values
            self.assertEqual(drive.initial_time, 0.0)
            self.assertEqual(drive.initial_value, 0.0)
            
            # Verify warnings were raised
            self.assertEqual(len(w), 2)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertTrue("<initial_time> is not set, assuming 0.0." in str(w[0].message))
            self.assertTrue(issubclass(w[1].category, UserWarning))
            self.assertTrue("<initial_value> is not set, assuming 0.0." in str(w[1].message))
    
    def test_step5_drive_caller_str_representation(self):
        """Test string representation of the drive caller"""
        # Test without idx
        drive = l.Step5DriveCaller(
            initial_time=self.initial_time,
            initial_value=self.initial_value,
            final_time=self.final_time,
            final_value=self.final_value
        )
        
        expected_str = f"step5, {self.initial_time}, {self.initial_value}, {self.final_time}, {self.final_value}"
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.Step5DriveCaller(
            idx=5,
            initial_time=self.initial_time,
            initial_value=self.initial_value,
            final_time=self.final_time,
            final_value=self.final_value
        )
        
        expected_str = f"drive caller: 5, step5, {self.initial_time}, {self.initial_value}, {self.final_time}, {self.final_value}"
        self.assertEqual(str(drive), expected_str)
    
    def test_step5_drive_caller_with_mbvars(self):
        """Test creating a Step5DriveCaller with MBVar objects"""
        # Create MBVar objects
        if 'init_time' not in l.declared_MBVars:
            initial_time_var = l.MBVar(name='init_time', var_type='real', expression=1.5)
        else:
            initial_time_var = l.declared_MBVars['init_time']
            
        if 'final_val' not in l.declared_MBVars:
            final_value_var = l.MBVar(name='final_val', var_type='real', expression=3.0)
        else:
            final_value_var = l.declared_MBVars['final_val']
        
        # Create with MBVar objects
        drive = l.Step5DriveCaller(
            initial_time=initial_time_var,
            initial_value=self.initial_value,
            final_time=self.final_time,
            final_value=final_value_var
        )
        
        # Check that MBVar references are stored correctly
        self.assertEqual(drive.initial_time, initial_time_var)
        self.assertEqual(drive.final_value, final_value_var)
        
        # Check string representation with MBVars
        expected_str = f"step5, {initial_time_var}, {self.initial_value}, {self.final_time}, {final_value_var}"
        self.assertEqual(str(drive), expected_str)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_step5_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing final_time
        with self.assertRaises(Exception):
            l.Step5DriveCaller(
                initial_time=self.initial_time,
                initial_value=self.initial_value,
                final_value=self.final_value
            )
        
        # Missing final_value
        with self.assertRaises(Exception):
            l.Step5DriveCaller(
                initial_time=self.initial_time,
                initial_value=self.initial_value,
                final_time=self.final_time
            )
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_step5_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # Invalid type for initial_time
        with self.assertRaises(Exception):
            l.Step5DriveCaller(
                initial_time="invalid",
                initial_value=self.initial_value,
                final_time=self.final_time,
                final_value=self.final_value
            )
        
        # Create an MBVar with a non-real type
        if 'test_non_real' not in l.declared_MBVars:
            non_real_var = l.MBVar(name='test_non_real', var_type='integer', expression=5)
        else:
            non_real_var = l.declared_MBVars['test_non_real']
        
        # Test with non-real MBVar for final_value
        with self.assertRaises(TypeError) as context:
            l.Step5DriveCaller(
                initial_time=self.initial_time,
                initial_value=self.initial_value,
                final_time=self.final_time,
                final_value=non_real_var
            )
        self.assertIn("must be a real number or an MBVar of type real", str(context.exception))


class TestStringDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Sample expressions for testing
        self.expression = "e^(-Time)*cos(2.*pi*Time)"
        self.var_expression = "e^(-Var)*cos(2.*pi*Time)"
    
    def test_string_drive_caller_creation_valid(self):
        """Test creating a StringDriveCaller with valid parameters"""
        # Basic creation with string expression
        drive = l.StringDriveCaller(expression=self.expression)
        
        # Verify properties are set correctly
        self.assertEqual(drive.expression, self.expression)
        self.assertEqual(drive.drive_type(), "string")
        
        # Create with idx
        drive = l.StringDriveCaller(idx=5, expression=self.expression)
        self.assertEqual(drive.idx, 5)

        # Test complex expression
        complex_expr = "1+cos(2*pi*Time)"
        drive = l.StringDriveCaller(expression=complex_expr)
        self.assertEqual(drive.expression, complex_expr)
        
        # Test expression with integer_eval
        eval_expr = "integer_eval(model::distance(CURR_NODE, CURR_NODE+1))"
        drive = l.StringDriveCaller(expression=eval_expr)
        self.assertEqual(drive.expression, eval_expr)
        
        # Test expression with multiple functions
        func_expr = "sin(Time)^2 + cos(Time)^2 + tan(Time/2) + sqrt(abs(Time))"
        drive = l.StringDriveCaller(expression=func_expr)
        self.assertEqual(drive.expression, func_expr)

    
    def test_string_drive_caller_with_mbvar(self):
        """Test creating a StringDriveCaller with an MBVar"""
        # Create MBVar for expression
        if 'expr_var' not in l.declared_MBVars:
            expr_var = l.MBVar(name='expr_var', var_type='string', expression=self.expression)
        else:
            expr_var = l.declared_MBVars['expr_var']
        
        # Create with MBVar
        drive = l.StringDriveCaller(expression=expr_var)
        
        # Check properties
        self.assertEqual(drive.expression, expr_var)
        
        # Check string representation
        expected_str = f'string, "{expr_var.expression}"'
        self.assertEqual(str(drive), expected_str)
        print(expected_str)
        print(str(drive))
    
    def test_string_drive_caller_str_representation(self):
        """Test string representation of the drive caller"""
        # Test without idx
        drive = l.StringDriveCaller(expression=self.expression)
        expected_str = f'string, "{self.expression}"'
        self.assertEqual(str(drive), expected_str)
        
        # Test with idx
        drive = l.StringDriveCaller(idx=10, expression=self.expression)
        expected_str = f'drive caller: 10, string, "{self.expression}"'
        self.assertEqual(str(drive), expected_str)

        # Test with variable expression
        drive = l.StringDriveCaller(expression=self.var_expression)
        expected_str = f'string, "{self.var_expression}"'
        self.assertEqual(str(drive), expected_str)
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_string_drive_caller_missing_required_field(self):
        """Test validation fails when required fields are missing"""
        # Missing expression
        with self.assertRaises(Exception):
            l.StringDriveCaller()
    
    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_string_drive_caller_invalid_types(self):
        """Test validation fails with invalid types"""
        # Invalid non-string expression
        with self.assertRaises(Exception):
            l.StringDriveCaller(expression=123)
        
        # Create an MBVar with a non-string type
        if 'non_string_var' not in l.declared_MBVars:
            non_string_var = l.MBVar(name='non_string_var', var_type='integer', expression=42)
        else:
            non_string_var = l.declared_MBVars['non_string_var']
        
        # Test with non-string MBVar for expression
        with self.assertRaises(TypeError) as context:
            l.StringDriveCaller(expression=non_string_var)
        self.assertIn("must be a string or an MBVar of type string", str(context.exception))
    
class TestTanhDriveCaller(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method"""
        # Create MBVar objects for testing if needed
        if 'amplitude_var' not in l.declared_MBVars:
            self.amplitude_var = l.MBVar(name='amplitude_var', var_type='real', expression=2.5)
        else:
            self.amplitude_var = l.declared_MBVars['amplitude_var']
    
    def test_tanh_drive_caller_creation_valid(self):
        """Test that TanhDriveCaller works with valid inputs"""
        # Create with all parameters
        tanh_drive = l.TanhDriveCaller(
            initial_time=1.0,
            amplitude=2.0,
            nd_slope=3.0,
            initial_value=0.0
        )
        self.assertIsInstance(tanh_drive, l.TanhDriveCaller)
        self.assertEqual(tanh_drive.initial_time, 1.0)
        self.assertEqual(tanh_drive.amplitude, 2.0)
        self.assertEqual(tanh_drive.nd_slope, 3.0)
        self.assertEqual(tanh_drive.initial_value, 0.0)
        
        # Create with idx
        tanh_drive = l.TanhDriveCaller(
            idx=10,
            initial_time=1.0,
            amplitude=2.0,
            nd_slope=3.0,
            initial_value=0.0
        )
        self.assertEqual(tanh_drive.idx, 10)
    
    def test_tanh_drive_caller_with_mbvars(self):
        """Test TanhDriveCaller with MBVar objects"""
        tanh_drive = l.TanhDriveCaller(
            initial_time=1.0,
            amplitude=self.amplitude_var,
            nd_slope=3.0,
            initial_value=0.0
        )
        self.assertEqual(tanh_drive.amplitude, self.amplitude_var)
    
    def test_tanh_drive_caller_default_warning(self):
        """Test warnings for default parameters"""
        # Test warning for default initial_time
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            tanh_drive = l.TanhDriveCaller(
                amplitude=2.0,
                nd_slope=3.0
            )
            self.assertTrue(any("<initial_time> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(tanh_drive.initial_time, 0.0)
        
        # Test warning for default initial_value
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            tanh_drive = l.TanhDriveCaller(
                amplitude=2.0,
                nd_slope=3.0,
                initial_time=1.0
            )
            self.assertTrue(any("<initial_value> is not set, assuming 0.0." in str(warning.message) for warning in w))
            self.assertEqual(tanh_drive.initial_value, 0.0)
    
    def test_tanh_drive_caller_str_representation(self):
        """Test string representation of TanhDriveCaller"""
        # Test without idx
        tanh_drive = l.TanhDriveCaller(
            initial_time=1.0,
            amplitude=2.0,
            nd_slope=3.0,
            initial_value=0.0
        )
        expected_str = "tanh, 1.0, 2.0, 3.0, 0.0"
        self.assertEqual(str(tanh_drive), expected_str)
        
        # Test with idx
        tanh_drive = l.TanhDriveCaller(
            idx=10,
            initial_time=1.0,
            amplitude=2.0,
            nd_slope=3.0,
            initial_value=0.0
        )
        expected_str = "drive caller: 10, tanh, 1.0, 2.0, 3.0, 0.0"
        self.assertEqual(str(tanh_drive), expected_str)
        
        # Test with MBVar
        tanh_drive = l.TanhDriveCaller(
            initial_time=1.0,
            amplitude=self.amplitude_var,
            nd_slope=3.0,
            initial_value=0.0
        )
        expected_str = f"tanh, 1.0, {self.amplitude_var}, 3.0, 0.0"
        self.assertEqual(str(tanh_drive), expected_str)
    
    def test_tanh_drive_caller_drive_type(self):
        """Test the drive_type method"""
        tanh_drive = l.TanhDriveCaller(
            initial_time=1.0,
            amplitude=2.0,
            nd_slope=3.0,
            initial_value=0.0
        )
        self.assertEqual(tanh_drive.drive_type(), "tanh")



    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_tanh_drive_caller_missing_required_field(self):
        """Test creating a TanhDriveCaller missing a required field"""
        # Missing amplitude
        with self.assertRaises(Exception):
            l.TanhDriveCaller(
                initial_time=1.0,
                nd_slope=3.0,
                initial_value=0.0
            )
        
        # Missing nd_slope
        with self.assertRaises(Exception):
            l.TanhDriveCaller(
                initial_time=1.0,
                amplitude=2.0,
                initial_value=0.0
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_tanh_drive_caller_invalid_types(self):
        """Test invalid types for TanhDriveCaller fields"""
        # Invalid type for amplitude
        with self.assertRaises(Exception):
            l.TanhDriveCaller(
                initial_time=1.0,
                amplitude="invalid string",
                nd_slope=3.0,
                initial_value=0.0
            )
        
        # Invalid MBVar type (using string MBVar)
        if 'string_var' not in l.declared_MBVars:
            string_var = l.MBVar(name='string_var', var_type='string', expression="test")
        else:
            string_var = l.declared_MBVars['string_var']
            
        with self.assertRaises(TypeError):
            l.TanhDriveCaller(
                initial_time=1.0,
                amplitude=2.0,
                nd_slope=string_var,  # String MBVar is invalid for numeric field
                initial_value=0.0
            )

class TestTimeDriveCaller(unittest.TestCase):
    def test_time_drive_caller_creation_valid_and_str_representation(self):
        """Test that TimeDriveCaller works with valid inputs"""
        # Create without idx
        time_drive = l.TimeDriveCaller()
        self.assertIsInstance(time_drive, l.TimeDriveCaller)
        expected_str = "time"
        self.assertEqual(str(time_drive), expected_str)
        
        # Create with idx
        time_drive = l.TimeDriveCaller(idx=10)
        self.assertIsInstance(time_drive, l.TimeDriveCaller)
        self.assertEqual(time_drive.idx, 10)
        expected_str = "drive caller: 10, time"
        self.assertEqual(str(time_drive), expected_str)
        
        # Create with MBVar for idx
        if 'idx_var' not in l.declared_MBVars:
            idx_var = l.MBVar(name='idx_var', var_type='integer', expression=5)
        else:
            idx_var = l.declared_MBVars['idx_var']
        time_drive = l.TimeDriveCaller(idx=idx_var)
        self.assertEqual(time_drive.idx, idx_var)
        expected_str = f"drive caller: {idx_var}, time"
        self.assertEqual(str(time_drive), expected_str)

class TestTimestepDriveCaller(unittest.TestCase):
    def test_timestep_drive_caller_creation_valid_and_str_representation(self):
        """Test that TimestepDriveCaller works with valid inputs"""
        # Create without idx
        timestep_drive = l.TimestepDriveCaller()
        self.assertIsInstance(timestep_drive, l.TimestepDriveCaller)
        expected_str = "timestep"
        self.assertEqual(str(timestep_drive), expected_str)

        # Create with idx
        timestep_drive = l.TimestepDriveCaller(idx=10)
        self.assertIsInstance(timestep_drive, l.TimestepDriveCaller)
        self.assertEqual(timestep_drive.idx, 10)
        expected_str = "drive caller: 10, timestep"
        self.assertEqual(str(timestep_drive), expected_str)

class TestUnitDriveCaller(unittest.TestCase):
    def test_unit_drive_caller_creation_valid_and_str_representation(self):
        """Test that UnitDriveCaller works with valid inputs"""
        # Create without idx
        unit_drive = l.UnitDriveCaller()
        self.assertIsInstance(unit_drive, l.UnitDriveCaller)
        expected_str = "unit"
        self.assertEqual(str(unit_drive), expected_str)
        
        # Create with idx
        unit_drive = l.UnitDriveCaller(idx=10)
        self.assertIsInstance(unit_drive, l.UnitDriveCaller)
        self.assertEqual(unit_drive.idx, 10)
        expected_str = "drive caller: 10, unit"
        self.assertEqual(str(unit_drive), expected_str)
        
        # Create with MBVar for idx
        if 'idx_var' not in l.declared_MBVars:
            idx_var = l.MBVar(name='idx_var', var_type='integer', expression=5)
        else:
            idx_var = l.declared_MBVars['idx_var']
        unit_drive = l.UnitDriveCaller(idx=idx_var)
        self.assertEqual(unit_drive.idx, idx_var)
        expected_str = f"drive caller: {idx_var}, unit"
        self.assertEqual(str(unit_drive), expected_str)

class TestLinearElastic(unittest.TestCase):
    def setUp(self):
        self.scalar_law = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e9)
        self.vector_3d_law = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW, stiffness=1e9)
        self.vector_6d_law = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW, stiffness=1e9)

    def test_const_law_name(self):
        self.assertEqual(self.scalar_law.const_law_header(), 'linear elastic')
        self.assertEqual(self.vector_3d_law.const_law_header(), 'linear elastic isotropic')
        self.assertEqual(self.vector_6d_law.const_law_header(), 'linear elastic isotropic')

    def test_str(self):
        self.assertEqual(str(self.scalar_law), f'{self.scalar_law.const_law_header()}, 1000000000.0')
        self.assertEqual(str(self.vector_3d_law), f'{self.vector_3d_law.const_law_header()}, 1000000000.0')
        self.assertEqual(str(self.vector_6d_law), f'{self.vector_6d_law.const_law_header()}, 1000000000.0')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_law_type(self):
        with self.assertRaises(Exception):
            l.LinearElastic(law_type='INVALID_LAW_TYPE', stiffness=1e9)

    def test_different_stiffness_values(self):
        small_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e-9)
        large_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e12)
        zero_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=0)
        self.assertEqual(small_stiffness.stiffness, 1e-9)
        self.assertEqual(large_stiffness.stiffness, 1e12)
        self.assertEqual(zero_stiffness.stiffness, 0)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_missing_arguments(self):
        with self.assertRaises(Exception):
            l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW)
        with self.assertRaises(Exception):
            l.LinearElastic(stiffness=1e9)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_stiffness_type(self):
        with self.assertRaises(Exception):
            l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness="invalid")
        with self.assertRaises(Exception):
            l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=None)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_extra_arguments(self):
        with self.assertRaises(Exception):
            l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e9, foo=1.0)

    def test_str_with_different_stiffness(self):
        small_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e-9)
        large_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=1e12)
        zero_stiffness = l.LinearElastic(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, stiffness=0)
        self.assertEqual(str(small_stiffness), f'{small_stiffness.const_law_header()}, 1e-09')
        self.assertEqual(str(large_stiffness), f'{large_stiffness.const_law_header()}, 1000000000000.0')
        self.assertEqual(str(zero_stiffness), f'{zero_stiffness.const_law_header()}, 0.0')

class TestLinearViscousGeneric(unittest.TestCase):
    def setUp(self):
        self.scalar_law = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e9)
        self.vector_3d_law = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW, viscosity=[[1e9, 0, 0], [0, 1e9, 0], [0, 0, 1e9]])
        self.vector_6d_law = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW, viscosity=[[1e9, 0, 0, 0, 0, 0], [0, 1e9, 0, 0, 0, 0], [0, 0, 1e9, 0, 0, 0], [0, 0, 0, 1e9, 0, 0], [0, 0, 0, 0, 1e9, 0], [0, 0, 0, 0, 0, 1e9]])

    def test_name(self):
        self.assertEqual(self.scalar_law.law_type, l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW)
        self.assertEqual(self.vector_3d_law.law_type, l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW)
        self.assertEqual(self.vector_6d_law.law_type, l.ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW)
        
    def test_const_law_name(self):
        self.assertEqual(self.scalar_law.const_law_name(), 'linear viscous generic')
        self.assertEqual(self.vector_3d_law.const_law_name(), 'linear viscous generic')
        self.assertEqual(self.vector_6d_law.const_law_name(), 'linear viscous generic')

    def test_str(self):
        self.assertEqual(str(self.scalar_law), f'{self.scalar_law.const_law_header()}, 1000000000.0')
        self.assertEqual(str(self.vector_3d_law), f'{self.vector_3d_law.const_law_header()},\n\t1000000000.0, 0.0, 0.0,\n\t0.0, 1000000000.0, 0.0,\n\t0.0, 0.0, 1000000000.0')
        self.assertEqual(str(self.vector_6d_law), f'{self.vector_6d_law.const_law_header()},\n\t1000000000.0, 0.0, 0.0, 0.0, 0.0, 0.0,\n\t0.0, 1000000000.0, 0.0, 0.0, 0.0, 0.0,\n\t0.0, 0.0, 1000000000.0, 0.0, 0.0, 0.0,\n\t0.0, 0.0, 0.0, 1000000000.0, 0.0, 0.0,\n\t0.0, 0.0, 0.0, 0.0, 1000000000.0, 0.0,\n\t0.0, 0.0, 0.0, 0.0, 0.0, 1000000000.0')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_law_type(self):
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(law_type='INVALID_LAW_TYPE', viscosity=1e9)

    def test_different_viscosity_values(self):
        small_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e-9)
        large_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e12)
        zero_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=0)
        self.assertEqual(small_viscosity.viscosity, 1e-9)
        self.assertEqual(large_viscosity.viscosity, 1e12)
        self.assertEqual(zero_viscosity.viscosity, 0)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_missing_arguments(self):
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW)
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(viscosity=1e9)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_viscosity_type(self):
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity="invalid")
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=None)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_extra_arguments(self):
        with self.assertRaises(Exception):
            l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e9, foo=1.0)

    def test_str_with_different_viscosity(self):
        small_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e-9)
        large_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=1e12)
        zero_viscosity = l.LinearViscousGeneric(law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW, viscosity=0)
        self.assertEqual(str(small_viscosity), f'{small_viscosity.const_law_header()}, 1e-09')
        self.assertEqual(str(large_viscosity), f'{large_viscosity.const_law_header()}, 1000000000000.0')
        self.assertEqual(str(zero_viscosity), f'{zero_viscosity.const_law_header()}, 0.0')

class TestLinearViscoelasticGeneric(unittest.TestCase):
    def test_valid_initialization_with_viscosity(self):
        stiffness = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        viscosity = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        law = l.LinearViscoelasticGeneric(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=stiffness,
            viscosity=viscosity
        )
        self.assertEqual(law.const_law_name(), 'linear viscoelastic generic')
        self.assertEqual(law.stiffness, stiffness)
        self.assertEqual(law.viscosity, viscosity)
        self.assertIsNone(law.factor)

    def test_valid_initialization_with_factor(self):
        stiffness = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        factor = 0.5
        law = l.LinearViscoelasticGeneric(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=stiffness,
            factor=factor
        )
        self.assertEqual(law.const_law_name(), 'linear viscoelastic generic')
        self.assertEqual(law.stiffness, stiffness)
        self.assertEqual(law.factor, factor)
        self.assertIsNone(law.viscosity)

    def test_str_representation_with_viscosity(self):
        stiffness = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        viscosity = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
        law = l.LinearViscoelasticGeneric(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=stiffness,
            viscosity=viscosity
        )
        expected_str = f'{law.const_law_header()},\n\t1.0, 0.0, 0.0,\n\t0.0, 1.0, 0.0,\n\t0.0, 0.0, 1.0,\n\t2.0, 0.0, 0.0,\n\t0.0, 2.0, 0.0,\n\t0.0, 0.0, 2.0'
        self.assertEqual(str(law), expected_str)

    def test_str_representation_with_factor(self):
        stiffness = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        factor = 0.5
        law = l.LinearViscoelasticGeneric(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=stiffness,
            factor=factor
        )
        expected_str = f'{law.const_law_header()},\n\t1.0, 0.0, 0.0,\n\t0.0, 1.0, 0.0,\n\t0.0, 0.0, 1.0, proportional, 0.5'
        self.assertEqual(str(law), expected_str)

    def test_invalid_stiffness_matrix(self):
        with self.assertRaises(ValueError):
            invalid_stiffness = [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]
            l.LinearViscoelasticGeneric(
                law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
                stiffness=invalid_stiffness,
                factor=0.5
            )

    def test_invalid_viscosity_matrix(self):
        with self.assertRaises(ValueError):
            invalid_viscosity = [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]
            l.LinearViscoelasticGeneric(
                law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
                stiffness=[[1.0]],
                viscosity=invalid_viscosity
            )

    def test_missing_factor_or_viscosity(self):
        with self.assertRaises(ValueError):
            l.LinearViscoelasticGeneric(
                law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
                stiffness=[[1.0]]
            )

    def test_both_viscosity_and_factor_provided(self):
        with self.assertRaises(ValueError):
            l.LinearViscoelasticGeneric(
                law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
                stiffness=[[1.0]],
                viscosity=[[2.0]],
                factor=0.5
            )

class TestPosition2(unittest.TestCase):
    def test_initialization_with_list(self):
        pos = l.Position('', [1.0, 2.0, 3.0])
        self.assertEqual(pos.relative_position, [1.0, 2.0, 3.0])
        pos2 = l.Position2(reference='', relative_position=[1.0, 2.0, 3.0])
        self.assertEqual(pos2.relative_position, [1.0, 2.0, 3.0])
        self.assertEqual(pos.relative_position, pos2.relative_position)

    def test_initialization_with_non_list(self):
        pos = l.Position('', 1.0)
        self.assertEqual(pos.relative_position, [1.0])
        pos2 = l.Position2(reference='', relative_position=1.0)
        self.assertEqual(pos2.relative_position, [1.0])
        self.assertEqual(pos.relative_position, pos2.relative_position)

    def test_string_representation_with_empty_reference(self):
        pos = l.Position('', [1.0, 2.0, 3.0])
        self.assertEqual(str(pos), '1.0, 2.0, 3.0')
        pos2 = l.Position2(reference='', relative_position=[1.0, 2.0, 3.0])
        self.assertEqual(str(pos2), '1.0, 2.0, 3.0')
        self.assertEqual(str(pos), str(pos2))

    def test_string_representation_with_non_empty_reference(self):
        pos = l.Position('global', [1.0, 2.0, 3.0])
        self.assertEqual(str(pos), 'reference, global, 1.0, 2.0, 3.0')
        pos2 = l.Position2(reference='global', relative_position=[1.0, 2.0, 3.0])
        self.assertEqual(str(pos2), 'reference, global, 1.0, 2.0, 3.0')
        self.assertEqual(str(pos), str(pos2))

    def test_isnull(self):
        pos = l.Position('', [ l.null()])
        self.assertTrue(pos.isnull())
        pos2 = l.Position2(reference='', relative_position=[l.null()])
        self.assertTrue(pos2.isnull())
        self.assertEqual(str(pos), str(pos2))

    def test_iseye(self):
        pos = l.Position('', [l.eye()])
        self.assertTrue(pos.iseye())
        pos2 = l.Position2(reference='', relative_position=[l.eye()])
        self.assertTrue(pos2.iseye())
        self.assertEqual(str(pos), str(pos2))

class TestReference2(unittest.TestCase):
    def test_initialization(self):
        pos2 = l.Position2(reference='', relative_position=[1.0, 2.0, 3.0])
        orient2 = l.Position2(reference='', relative_position=[0.0, 0.0, 1.0])
        vel2 = l.Position2(reference='', relative_position=[0.0, 0.0, 0.0])
        angvel2 = l.Position2(reference='', relative_position=[0.1, 0.1, 0.1])
        ref2 = l.Reference2(idx=1, position=pos2, orientation=orient2, velocity=vel2, angular_velocity=angvel2)
        self.assertEqual(str(ref2), 'reference: 1, \n\t1.0, 2.0, 3.0,\n\t0.0, 0.0, 1.0,\n\t0.0, 0.0, 0.0,\n\t0.1, 0.1, 0.1;\n')

    def test_against_Reference(self):
        pos = l.Position('', [1.0, 2.0, 3.0])
        orient = l.Position('', [0.0, 0.0, 1.0])
        vel = l.Position('', [0.0, 0.0, 0.0])
        angvel = l.Position('', [0.1, 0.1, 0.1])
        ref = l.Reference(1, pos, orient, vel, angvel)
        pos2 = l.Position2(reference='', relative_position=[1.0, 2.0, 3.0])
        orient2 = l.Position2(reference='', relative_position=[0.0, 0.0, 1.0])
        vel2 = l.Position2(reference='', relative_position=[0.0, 0.0, 0.0])
        angvel2 = l.Position2(reference='', relative_position=[0.1, 0.1, 0.1])
        ref2 = l.Reference2(idx=1, position=pos2, orientation=orient2, velocity=vel2, angular_velocity=angvel2)
        self.assertEqual(str(ref), str(ref2))

class TestAngularAcceleration(unittest.TestCase):

    def test_valid_input(self):
        """Test that AngularAcceleration works with valid input"""
        # Valid instance of ConstDriveCaller
        const_drive = l.ConstDriveCaller(const_value=5.0)
        
        # Create an AngularAcceleration instance with valid inputs
        angular_accel = l.AngularAcceleration(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            acceleration=const_drive
        )
        
        expected_output = '''joint: 1, angular acceleration,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0;\n'''
        self.assertEqual(str(angular_accel), expected_output)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_relative_direction_length(self):
        """Test that AngularAcceleration raises an error for invalid relative_direction length"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # relative_direction must have exactly 3 elements, expect failure
        with self.assertRaises(Exception):
            l.AngularAcceleration(
                idx=1,
                node_label=1,
                relative_direction=[1, 0],  # Invalid length
                acceleration=const_drive
            )

    def test_invalid_relative_direction_magnitude(self):
        """Test that AngularAcceleration raises an error for non-unit vector relative_direction"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # relative_direction must be a unit vector (magnitude = 1), expect failure
        with self.assertRaises(ValueError):
            l.AngularAcceleration(
                idx=1,
                node_label=1,
                relative_direction=[2, 0, 0],  # Invalid magnitude (not a unit vector)
                acceleration=const_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_acceleration_type(self):
        """Test that AngularAcceleration raises an error for invalid acceleration type"""
        # Pass an invalid type for acceleration
        with self.assertRaises(Exception):
            l.AngularAcceleration(
                idx=1,
                node_label=1,
                relative_direction=[1, 0, 0],
                acceleration=5  # Invalid type, should be DriveCaller or its subclass
            )

    def test_optional_output(self):
        """Test that the 'output' field is optional and defaults to 'yes'"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # Create AngularAcceleration without specifying output
        angular_accel = l.AngularAcceleration(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            acceleration=const_drive
        )
        
        expected_output = '''joint: 1, angular acceleration,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0;\n'''
        self.assertEqual(str(angular_accel), expected_output)

    def test_custom_output(self):
        """Test that the 'output' field is properly set when customized"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # Create AngularAcceleration with custom output
        angular_accel = l.AngularAcceleration(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            acceleration=const_drive,
            output='no'
        )
        
        expected_output = '''joint: 1, angular acceleration,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0,\n\toutput, no;\n'''
        self.assertEqual(str(angular_accel), expected_output)

class TestAngularVelocity(unittest.TestCase):

    def test_valid_input(self):
        """Test that AngularVelocity works with valid input"""
        # Valid instance of ConstDriveCaller
        const_drive = l.ConstDriveCaller(const_value=5.0)
        
        # Create an AngularVelocity instance with valid inputs
        angular_vel = l.AngularVelocity(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            velocity=const_drive
        )
        
        expected_output = '''joint: 1, angular velocity,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0;\n'''
        self.assertEqual(str(angular_vel), expected_output)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_relative_direction_length(self):
        """Test that AngularVelocity raises an error for invalid relative_direction length"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # relative_direction must have exactly 3 elements, expect failure
        with self.assertRaises(Exception):
            l.AngularVelocity(
                idx=1,
                node_label=1,
                relative_direction=[1, 0],  # Invalid length
                velocity=const_drive
            )

    def test_invalid_relative_direction_magnitude(self):
        """Test that AngularVelocity raises an error for non-unit vector relative_direction"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # relative_direction must be a unit vector (magnitude = 1), expect failure
        with self.assertRaises(ValueError):
            l.AngularVelocity(
                idx=1,
                node_label=1,
                relative_direction=[2, 0, 0],  # Invalid magnitude (not a unit vector)
                velocity=const_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_velocity_type(self):
        """Test that AngularVelocity raises an error for invalid velocity type"""
        # Pass an invalid type for velocity
        with self.assertRaises(Exception):
            l.AngularVelocity(
                idx=1,
                node_label=1,
                relative_direction=[1, 0, 0],
                velocity=5  # Invalid type, should be DriveCaller or its subclass
            )

    def test_optional_output(self):
        """Test that the 'output' field is optional and defaults to 'yes'"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # Create AngularVelocity without specifying output
        angular_vel = l.AngularVelocity(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            velocity=const_drive
        )
        
        expected_output = '''joint: 1, angular velocity,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0;\n'''
        self.assertEqual(str(angular_vel), expected_output)

    def test_custom_output(self):
        """Test that the 'output' field is properly set when customized"""
        const_drive = l.ConstDriveCaller(const_value=5.0)

        # Create AngularVelocity with custom output
        angular_vel = l.AngularVelocity(
            idx=1,
            node_label=1,
            relative_direction=[1, 0, 0],
            velocity=const_drive,
            output='no'
        )
        
        expected_output = '''joint: 1, angular velocity,\n\t1, [1.0, 0.0, 0.0],\n\tconst, 5.0,\n\toutput, no;\n'''
        self.assertEqual(str(angular_vel), expected_output)

class TestAxialRotation(unittest.TestCase):

    def test_valid_input(self):
        """Test that AxialRotation works with valid input"""
        # Valid instances of Position2 and DriveCaller
        position1 = l.Position2(reference='global', relative_position=[0, 0, 0])
        orientation1 = l.Position2(reference='global', relative_position=[1, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        orientation2 = l.Position2(reference='global', relative_position=[0, 1, 0])
        drive_caller = l.ConstDriveCaller(const_value=5.0)
        
        # Create an AxialRotation instance with valid inputs
        axial_rot = l.AxialRotation(
            idx=1,
            node_1_label=1,
            position_1=position1,
            orientation_mat_1=orientation1,
            node_2_label=2,
            position_2=position2,
            orientation_mat_2=orientation2,
            angular_velocity=drive_caller
        )
        
        expected_output = (
            f"{axial_rot.element_header()}, axial rotation,\n"
            f"\t{axial_rot.node_1_label},\n"
            f"\t\tposition, {axial_rot.position_1},\n"
            f"\t\torientation, {axial_rot.orientation_mat_1},\n"
            f"\t{axial_rot.node_2_label},\n"
            f"\t\tposition, {axial_rot.position_2},\n"
            f"\t\torientation, {axial_rot.orientation_mat_2},\n"
            f"\t{axial_rot.angular_velocity}"
            f"{axial_rot.element_footer()}"
        )
        self.maxDiff=None
        self.assertEqual(str(axial_rot), expected_output)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_position_type(self):
        """Test that AxialRotation raises an error for invalid position type"""
        # Valid instances of Position2 and DriveCaller
        orientation1 = l.Position2(reference='global', relative_position=[1, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        orientation2 = l.Position2(reference='global', relative_position=[0, 1, 0])
        drive_caller = l.ConstDriveCaller(const_value=5.0)

        # position_1 must be of type Position2
        with self.assertRaises(Exception):
            l.AxialRotation(
                idx=1,
                node_1_label=1,
                position_1='invalid_position',  # Invalid type
                orientation_mat_1=orientation1,
                node_2_label=2,
                position_2=position2,
                orientation_mat_2=orientation2,
                angular_velocity=drive_caller
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_orientation_type(self):
        """Test that AxialRotation raises an error for invalid orientation type"""
        # Valid instances of Position2 and DriveCaller
        position1 = l.Position2(reference='global', relative_position=[0, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        drive_caller = l.ConstDriveCaller(const_value=5.0)

        # orientation_mat_1 must be of type Position2
        with self.assertRaises(Exception):
            l.AxialRotation(
                idx=1,
                node_1_label=1,
                position_1=position1,
                orientation_mat_1='invalid_orientation',  # Invalid type
                node_2_label=2,
                position_2=position2,
                orientation_mat_2=l.Position2(reference='global', relative_position=[0, 1, 0]),
                angular_velocity=drive_caller
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_angular_velocity_type(self):
        """Test that AxialRotation raises an error for invalid angular_velocity type"""
        # Valid instances of Position2
        position1 = l.Position2(reference='global', relative_position=[0, 0, 0])
        orientation1 = l.Position2(reference='global', relative_position=[1, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        orientation2 = l.Position2(reference='global', relative_position=[0, 1, 0])

        # angular_velocity must be of type DriveCaller or DriveCaller2
        with self.assertRaises(Exception):
            l.AxialRotation(
                idx=1,
                node_1_label=1,
                position_1=position1,
                orientation_mat_1=orientation1,
                node_2_label=2,
                position_2=position2,
                orientation_mat_2=orientation2,
                angular_velocity='invalid_velocity'  # Invalid type
            )

    def test_optional_output(self):
        """Test that the 'output' field is optional and defaults to 'yes'"""
        # Valid instances of Position2 and DriveCaller
        position1 = l.Position2(reference='global', relative_position=[0, 0, 0])
        orientation1 = l.Position2(reference='global', relative_position=[1, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        orientation2 = l.Position2(reference='global', relative_position=[0, 1, 0])
        drive_caller = l.ConstDriveCaller(const_value=5.0)

        # Create AxialRotation without specifying output
        axial_rot = l.AxialRotation(
            idx=1,
            node_1_label=1,
            position_1=position1,
            orientation_mat_1=orientation1,
            node_2_label=2,
            position_2=position2,
            orientation_mat_2=orientation2,
            angular_velocity=drive_caller
        )
        expected_output = (
            f"{axial_rot.element_header()}, axial rotation,\n"
            f"\t{axial_rot.node_1_label},\n"
            f"\t\tposition, {axial_rot.position_1},\n"
            f"\t\torientation, {axial_rot.orientation_mat_1},\n"
            f"\t{axial_rot.node_2_label},\n"
            f"\t\tposition, {axial_rot.position_2},\n"
            f"\t\torientation, {axial_rot.orientation_mat_2},\n"
            f"\t{axial_rot.angular_velocity}"
            f"{axial_rot.element_footer()}"
        )  
        self.maxDiff=None      
        self.assertEqual(str(axial_rot), expected_output)

    def test_custom_output(self):
        """Test that the 'output' field is properly set when customized"""
        # Valid instances of Position2 and DriveCaller
        position1 = l.Position2(reference='global', relative_position=[0, 0, 0])
        orientation1 = l.Position2(reference='global', relative_position=[1, 0, 0])
        position2 = l.Position2(reference='global', relative_position=[1, 1, 1])
        orientation2 = l.Position2(reference='global', relative_position=[0, 1, 0])
        drive_caller = l.ConstDriveCaller(const_value=5.0)

        # Create AxialRotation with custom output
        axial_rot = l.AxialRotation(
            idx=1,
            node_1_label=1,
            position_1=position1,
            orientation_mat_1=orientation1,
            node_2_label=2,
            position_2=position2,
            orientation_mat_2=orientation2,
            angular_velocity=drive_caller,
            output='no'
        )
        
        expected_output = (
            f"{axial_rot.element_header()}, axial rotation,\n"
            f"\t{axial_rot.node_1_label},\n"
            f"\t\tposition, {axial_rot.position_1},\n"
            f"\t\torientation, {axial_rot.orientation_mat_1},\n"
            f"\t{axial_rot.node_2_label},\n"
            f"\t\tposition, {axial_rot.position_2},\n"
            f"\t\torientation, {axial_rot.orientation_mat_2},\n"
            f"\t{axial_rot.angular_velocity}"
            f"{axial_rot.element_footer()}"
        )  
        self.maxDiff=None
        self.assertEqual(str(axial_rot), expected_output)

class TestBeamSlider(unittest.TestCase):
    def setUp(self):
        # Define Position2 instances
        self.position1 = l.Position2(
            relative_position=[[0.0, 0.0, 0.0]], 
            reference='global'
        )
        self.position2 = l.Position2(
            relative_position=[[1.0, 0.0, 0.0]], 
            reference='node'
        )
        self.position3 = l.Position2(
            relative_position=[[0.0, 1.0, 0.0]], 
            reference='other node'
        )

        # Define Constitutive Laws
        self.elastic_law = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=2000.0
        )

        # Define Beams
        self.beam = l.Beam(
            idx=1,
            nodes=[1, 2, 3],
            positions=[self.position1, self.position2, self.position3],
            orientations=[self.position1, self.position2, self.position3],
            const_laws_orientations=[self.position1, self.position2],
            const_laws=[self.elastic_law, self.elastic_law],
        )
    
    def test_valid_input(self):
        beam_slider = l.BeamSlider(
            idx=1,
            slider_node_label=1,
            position=self.position1,
            orientation=self.position2,
            slider_type='classic',
            beam_number=1,
            three_node_beam=self.beam,
            first_node_offset=self.position1,
            first_node_orientation=self.position2,
            mid_node_offset=self.position2,
            mid_node_orientation=self.position3,
            end_node_offset=self.position3,
            end_node_orientation=self.position1,
            initial_beam=self.beam,
            initial_node=None,
            smearing_factor=0.5
        )
        expected_str = (
            f"joint: {beam_slider.idx}, kinematic,\n"
            f"\t{beam_slider.slider_node_label},\n"
            f"\t\t{beam_slider.position},\n"
            f"\t\thinge, {beam_slider.orientation},\n"
            f"\ttype, {beam_slider.slider_type},\n"
            f"\t{beam_slider.beam_number},\n"
            f"\t\t{beam_slider.three_node_beam}"[:-2] + ",\n"  # Remove ';\n' and add ',\n'
            f"\t\t\t{beam_slider.first_node_offset},\n"
            f"\t\thinge, {beam_slider.first_node_orientation},\n"
            f"\t\t\t{beam_slider.mid_node_offset},\n"
            f"\t\thinge, {beam_slider.mid_node_orientation},\n"
            f"\t\t\t{beam_slider.end_node_offset},\n"
            f"\t\thinge, {beam_slider.end_node_orientation},\n"
            f"\tinitial beam, {beam_slider.initial_beam}"[:-2] + ",\n"  # Remove ';\n' and add ',\n'
            f"\tsmearing, {beam_slider.smearing_factor};\n"
        )
        self.maxDiff=None
        self.assertEqual(str(beam_slider), expected_str)
    
    def test_invalid_slider_type(self):
        with self.assertRaises(ValueError):
            l.BeamSlider(
                idx=1,
                slider_node_label=1,
                position=self.position1,
                orientation=self.position2,
                slider_type='invalid_type',
                beam_number=1,
                three_node_beam=self.beam,
                first_node_offset=self.position1,
                first_node_orientation=self.position2,
                mid_node_offset=self.position2,
                mid_node_orientation=self.position3,
                end_node_offset=self.position3,
                end_node_orientation=self.position1,
                initial_beam=self.beam,
                initial_node=None,
                smearing_factor=0.5
            )
    
    def test_optional_fields(self):
        beam_slider = l.BeamSlider(
            idx=1,
            slider_node_label=1,
            position=self.position1,
            orientation=None,
            slider_type=None,
            beam_number=1,
            three_node_beam=self.beam,
            first_node_offset=self.position1,
            first_node_orientation=None,
            mid_node_offset=self.position2,
            mid_node_orientation=None,
            end_node_offset=self.position3,
            end_node_orientation=None,
            initial_beam=None,
            initial_node=None,
            smearing_factor=None
        )
        expected_str = (
            f"{beam_slider.element_header()}, kinematic,\n"
            f"\t{beam_slider.slider_node_label},\n"
            f"\t\t{beam_slider.position},\n"
            f"\t{beam_slider.beam_number},\n"
            f"\t\t{beam_slider.three_node_beam}"[:-2] + ",\n"  # Remove ';\n' and add ',\n'
            f"\t\t\t{beam_slider.first_node_offset},\n"
            f"\t\t\t{beam_slider.mid_node_offset},\n"
            f"\t\t\t{beam_slider.end_node_offset}"
            f"{beam_slider.element_footer()}"
        )
        self.maxDiff = None
        self.assertEqual(str(beam_slider), expected_str)

    def test_invalid_mid_node_offset_type(self):
        with self.assertRaises(ValueError):
            l.BeamSlider(
                idx=1,
                slider_node_label=1,
                position=self.position1,
                orientation=self.position2,
                slider_type='classic',
                beam_number=1,
                three_node_beam=self.beam,
                first_node_offset=self.position1,
                first_node_orientation=self.position2,
                mid_node_offset='invalid_offset',
                mid_node_orientation=self.position3,
                end_node_offset=self.position3,
                end_node_orientation=self.position1,
                initial_beam=self.beam,
                initial_node=None,
                smearing_factor=0.5
            )

    def test_invalid_end_node_offset_type(self):
        with self.assertRaises(ValueError):
            l.BeamSlider(
                idx=1,
                slider_node_label=1,
                position=self.position1,
                orientation=self.position2,
                slider_type='classic',
                beam_number=1,
                three_node_beam=self.beam,
                first_node_offset=self.position1,
                first_node_orientation=self.position2,
                mid_node_offset=self.position2,
                mid_node_orientation=self.position3,
                end_node_offset='invalid_offset',
                end_node_orientation=self.position1,
                initial_beam=self.beam,
                initial_node=None,
                smearing_factor=0.5
            )

    def test_invalid_initial_node_type(self):
        with self.assertRaises(ValueError):
            l.BeamSlider(
                idx=1,
                slider_node_label=1,
                position=self.position1,
                orientation=self.position2,
                slider_type='classic',
                beam_number=1,
                three_node_beam=self.beam,
                first_node_offset=self.position1,
                first_node_orientation=self.position2,
                mid_node_offset=self.position2,
                mid_node_orientation=self.position3,
                end_node_offset=self.position3,
                end_node_orientation=self.position1,
                initial_beam=self.beam,
                initial_node='invalid_node',
                smearing_factor=0.5
            )

class TestBrake(unittest.TestCase):
    def setUp(self):
        self.position1 = l.Position2(relative_position=[[0.0, 0.0, 0.0]], reference='global')
        self.position2 = l.Position2(relative_position=[[1.0, 0.0, 0.0]], reference='node')
        self.normal_force = l.ConstDriveCaller(const_value=1000.0)

    def test_valid_brake(self):
        brake = l.Brake(
            idx=1,
            node_1_label=1,
            position_1=self.position1,
            node_2_label=2,
            position_2=self.position2,
            average_radius=0.5,
            friction_model="modlugre",
            shape_function="tanh",
            normal_force=self.normal_force
        )
        self.assertIsInstance(brake, l.Brake)

    def test_str_representation(self):
        brake = l.Brake(
            idx=1,
            node_1_label=1,
            position_1=self.position1,
            node_2_label=2,
            position_2=self.position2,
            average_radius=0.5,
            friction_model="modlugre",
            shape_function="tanh",
            normal_force=self.normal_force
        )
        expected_str = (
            "joint: 1, brake,\n"
            "\t1, reference, global, [0.0, 0.0, 0.0],\n"
            "\t2, reference, node, [1.0, 0.0, 0.0],\n"
            "\tfriction, 0.5,\n"
            "\t\tmodlugre,\n"
            "\t\ttanh,\n"
            "\tconst, 1000.0;\n"
        )
        self.assertEqual(str(brake), expected_str)

    def test_with_optional_fields(self):
        brake = l.Brake(
            idx=1,
            node_1_label=1,
            position_1=self.position1,
            orientation_mat_1=self.position2,
            node_2_label=2,
            position_2=self.position2,
            orientation_mat_2=self.position1,
            average_radius=0.5,
            preload=100,
            friction_model="modlugre",
            shape_function="tanh",
            normal_force=self.normal_force
        )
        self.assertIsInstance(brake, l.Brake)
        self.assertEqual(brake.preload, 100)

class TestCardanoPin(unittest.TestCase):

    def test_abstract_class(self):
        """Check that user can't create abstract classes, which are used only to share functionality (can't be part of MBDyn output)"""
        with self.assertRaises(TypeError):
            e = l.Element2()

    def setUp(self):
        self.node_label = 5
        self.relative_position = [0.0, 1.0, 2.0]
        self.absolute_pin_position = l.Position2(reference='global', relative_position=[3.0, 4.0, 5.0])

        # Optional values for testing with orientations
        self.relative_orientation_matrix = l.Position2(
            relative_position=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            reference=''
        )
        self.absolute_orientation_matrix = l.Position2(
            relative_position=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )

        # Initialize CardanoPin with required fields
        self.cardano_pin = l.CardanoPin(
            idx=1,
            output='yes',
            node_label=self.node_label,
            position=l.Position2(relative_position=self.relative_position, reference='global'),
            absolute_pin_position=self.absolute_pin_position
        )

    def test_initialization(self):
        # Test that the CardanoPin initializes correctly with provided values
        self.assertEqual(self.cardano_pin.node_label, self.node_label)
        self.assertIsInstance(self.cardano_pin.position, l.Position2)
        self.assertEqual(self.cardano_pin.position.relative_position, self.relative_position)
        self.assertEqual(self.cardano_pin.position.reference, 'global')
        self.assertEqual(self.cardano_pin.absolute_pin_position, self.absolute_pin_position)
        self.assertIsNone(self.cardano_pin.orientation_mat)
        self.assertIsNone(self.cardano_pin.absolute_pin_orientation_mat)
        self.assertEqual(self.cardano_pin.idx, 1)
        self.assertEqual(self.cardano_pin.output, 'yes')

    def test_str_representation_without_optional(self):
        # Test the string output when optional orientation matrices are not provided
        expected_str = (
            f'{self.cardano_pin.element_header()}, cardano pin,\n\t{self.node_label},'
            f'\n\t\tposition, {self.cardano_pin.position},'
            f'\n\tposition, {self.cardano_pin.absolute_pin_position}'
            f'{self.cardano_pin.element_footer()}'
        )
        self.assertEqual(str(self.cardano_pin), expected_str)

    def test_str_representation_with_optional(self):
        # Initialize with optional orientation matrices
        cardano_pin_with_orientation = l.CardanoPin(
            idx=2,
            output='no',
            node_label=self.node_label,
            position=l.Position2(relative_position=self.relative_position, reference='global'),
            orientation_mat=self.relative_orientation_matrix,
            absolute_pin_position=self.absolute_pin_position,
            absolute_pin_orientation_mat=self.absolute_orientation_matrix
        )

        expected_str = (
            f'{cardano_pin_with_orientation.element_header()}, cardano pin,\n\t{self.node_label},'
            f'\n\t\tposition, {cardano_pin_with_orientation.position},'
            f'\n\t\torientation, {self.relative_orientation_matrix},'
            f'\n\tposition, {cardano_pin_with_orientation.absolute_pin_position},'
            f'\n\torientation, {self.absolute_orientation_matrix}'
            f'{cardano_pin_with_orientation.element_footer()}'
        )
        self.assertEqual(str(cardano_pin_with_orientation), expected_str)

    def test_optional_none_handling(self):
        # Test to check that None is handled correctly for optional orientation matrices
        cardano_pin_without_orientation = l.CardanoPin(
            idx=3,
            node_label=self.node_label,
            position=l.Position2(relative_position=self.relative_position, reference=''),
            absolute_pin_position=self.absolute_pin_position
        )

        # Assert that orientation matrices remain None
        self.assertIsNone(cardano_pin_without_orientation.orientation_mat)
        self.assertIsNone(cardano_pin_without_orientation.absolute_pin_orientation_mat)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_node_label(self):
        # Test that passing invalid node_label raises the appropriate error
        with self.assertRaises(Exception):
            l.CardanoPin(
                idx=4,
                node_label="invalid_label",  # Invalid type for node_label
                position=l.Position2(relative_position=self.relative_position, reference='global'),
                absolute_pin_position=self.absolute_pin_position
            )

    def test_invalid_position(self):
        # Test invalid Position2 for position and absolute_pin_position
        with self.assertRaises(ValueError):
            l.CardanoPin(
                idx=5,
                node_label=self.node_label,
                position=l.Position2(relative_position='invalid_value', reference='global'),
                absolute_pin_position=self.absolute_pin_position
            )

    def test_isnull_function(self):
        # Test if the `isnull()` function works correctly in Position2
        null_position = l.Position2(relative_position=[l.null()], reference='')
        self.assertTrue(null_position.isnull())

    def test_iseye_function(self):
        # Test if the `iseye()` function works correctly in Position2
        eye_position = l.Position2(relative_position=[l.eye()], reference='')
        self.assertTrue(eye_position.iseye())

class TestCardanoRotation(unittest.TestCase):

    def setUp(self):
        self.node_1_label = 1
        self.node_2_label = 2

        # Optional values for testing with orientations
        self.orientation_matrix_1 = l.Position2(
            relative_position=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            reference=''
        )
        self.orientation_matrix_2 = l.Position2(
            relative_position=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )

        # Initialize CardanoRotation with required fields
        self.cardano_rotation = l.CardanoRotation(
            idx=1,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )

    def test_initialization(self):
        # Test that the CardanoRotation initializes correctly with provided values
        self.assertEqual(self.cardano_rotation.node_1_label, self.node_1_label)
        self.assertEqual(self.cardano_rotation.node_2_label, self.node_2_label)
        self.assertIsNone(self.cardano_rotation.orientation_mat_1)
        self.assertIsNone(self.cardano_rotation.orientation_mat_2)
        self.assertEqual(self.cardano_rotation.idx, 1)
        self.assertEqual(self.cardano_rotation.output, 'yes')

    def test_str_representation_without_optional(self):
        # Test the string output when optional orientation matrices are not provided
        expected_str = (
            f'{self.cardano_rotation.element_header()}, cardano rotation,\n\t{self.node_1_label},'
            f'\n\t{self.node_2_label}'
            f'{self.cardano_rotation.element_footer()}'
        )
        self.assertEqual(str(self.cardano_rotation), expected_str)

    def test_str_representation_with_optional(self):
        # Initialize with optional orientation matrices
        cardano_rotation_with_orientation = l.CardanoRotation(
            idx=2,
            output='no',
            node_1_label=self.node_1_label,
            orientation_mat_1=self.orientation_matrix_1,
            node_2_label=self.node_2_label,
            orientation_mat_2=self.orientation_matrix_2
        )

        expected_str = (
            f'{cardano_rotation_with_orientation.element_header()}, cardano rotation,\n\t{self.node_1_label},'
            f'\n\t\torientation, {self.orientation_matrix_1},'
            f'\n\t{self.node_2_label},'
            f'\n\t\torientation, {self.orientation_matrix_2}'
            f'{cardano_rotation_with_orientation.element_footer()}'
        )
        self.assertEqual(str(cardano_rotation_with_orientation), expected_str)

    def test_optional_none_handling(self):
        # Test to check that None is handled correctly for optional orientation matrices
        cardano_rotation_without_orientation = l.CardanoRotation(
            idx=3,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )

        # Assert that orientation matrices remain None
        self.assertIsNone(cardano_rotation_without_orientation.orientation_mat_1)
        self.assertIsNone(cardano_rotation_without_orientation.orientation_mat_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_invalid_node_labels(self):
        # Test that passing invalid node labels raises the appropriate error
        with self.assertRaises(Exception):
            l.CardanoRotation(
                idx=4,
                node_1_label="invalid_label",  # Invalid type for node_1_label
                node_2_label=self.node_2_label
            )

        with self.assertRaises(Exception):
            l.CardanoRotation(
                idx=5,
                node_1_label=self.node_1_label,
                node_2_label="invalid_label"  # Invalid type for node_2_label
            )

class TestDeformableAxial(unittest.TestCase):

    def setUp(self):
        self.node_1_label = 1
        self.node_2_label = 2

        # Optional values for testing with positions and orientations
        self.position_1 = l.Position2(
            relative_position=[1.0, 0.0, 0.0],
            reference=''
        )
        self.orientation_mat_1 = l.Position2(
            relative_position=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )
        self.position_2 = l.Position2(
            relative_position=[0.0, 1.0, 0.0],
            reference='global'
        )
        self.orientation_mat_2 = l.Position2(
            relative_position=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )

        # Example constitutive laws
        self.linear_elastic = l.LinearElastic(
            idx=1,
            law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW,
            stiffness=2000
        )
        
        self.linear_elastic_generic = l.LinearElasticGeneric(
            idx=2,
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=[[1.0, 0.3, 0.0], [0.3, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )
        
        self.named_const_law = l.NamedConstitutiveLaw("example_named_law")
        self.named_const_law2 = l.NamedConstitutiveLaw(["example_named_law", 1000.0])

        # Initialize DeformableAxial with required fields
        self.deformable_axial = l.DeformableAxial(
            idx=1,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.linear_elastic
        )

    def test_initialization(self):
        # Test that the DeformableAxial initializes correctly with provided values
        self.assertEqual(self.deformable_axial.node_1_label, self.node_1_label)
        self.assertEqual(self.deformable_axial.node_2_label, self.node_2_label)
        self.assertIsNone(self.deformable_axial.position_1)
        self.assertIsNone(self.deformable_axial.orientation_mat_1)
        self.assertIsNone(self.deformable_axial.position_2)
        self.assertIsNone(self.deformable_axial.orientation_mat_2)
        self.assertEqual(self.deformable_axial.const_law, self.linear_elastic)
        self.assertEqual(self.deformable_axial.idx, 1)
        self.assertEqual(self.deformable_axial.output, 'yes')

    def test_str_representation_without_optional(self):
        # Test the string output when optional positions and orientations are not provided
        expected_str = (
            f'{self.deformable_axial.element_header()}, deformable axial,\n\t{self.node_1_label},'
            f'\n\t{self.node_2_label},'
            f'\n\t{self.linear_elastic}'
            f'{self.deformable_axial.element_footer()}'
        )
        self.assertEqual(str(self.deformable_axial), expected_str)

    def test_str_representation_with_optional(self):
        # Initialize with optional positions and orientations
        deformable_axial_with_optional = l.DeformableAxial(
            idx=2,
            output='no',
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2,
            const_law=self.linear_elastic_generic
        )

        expected_str = (
            f'{deformable_axial_with_optional.element_header()}, deformable axial,\n\t{self.node_1_label},'
            f'\n\t\tposition, {self.position_1},'
            f'\n\t\torientation, {self.orientation_mat_1},'
            f'\n\t{self.node_2_label},'
            f'\n\t\tposition, {self.position_2},'
            f'\n\t\torientation, {self.orientation_mat_2},'
            f'\n\t{self.linear_elastic_generic}'
            f'{deformable_axial_with_optional.element_footer()}'
        )
        self.assertEqual(str(deformable_axial_with_optional), expected_str)

    def test_optional_none_handling(self):
        # Test to check that None is handled correctly for optional positions and orientations
        deformable_axial_without_optional = l.DeformableAxial(
            idx=3,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.named_const_law2
        )

        # Assert that optional parameters remain None
        self.assertIsNone(deformable_axial_without_optional.position_1)
        self.assertIsNone(deformable_axial_without_optional.orientation_mat_1)
        self.assertIsNone(deformable_axial_without_optional.position_2)
        self.assertIsNone(deformable_axial_without_optional.orientation_mat_2)

    def test_invalid_node_labels(self):
        # Test that passing invalid node labels raises the appropriate error
        with self.assertRaises(Exception):
            l.DeformableAxial(
                idx=4,
                node_1_label="invalid_label",  # Invalid type for node_1_label
                node_2_label=self.node_2_label,
                const_law=self.linear_elastic
            )

        with self.assertRaises(Exception):
            l.DeformableAxial(
                idx=5,
                node_1_label=self.node_1_label,
                node_2_label="invalid_label",  # Invalid type for node_2_label
                const_law=self.linear_elastic
            )

    def test_invalid_const_law(self):
        # Test that passing an invalid const_law raises the appropriate error
        with self.assertRaises(Exception):
            l.DeformableAxial(
                idx=6,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                const_law="invalid_const_law"  # Invalid type for const_law, users have to use NamedConstitutiveLaw for custom Const Laws
            )

    def test_named_constitutive_law(self):
        # Test with NamedConstitutiveLaw
        deformable_axial_with_named_law = l.DeformableAxial(
            idx=7,
            output='yes',
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.named_const_law
        )

        expected_str = (
            f'{deformable_axial_with_named_law.element_header()}, deformable axial,\n\t{self.node_1_label},'
            f'\n\t{self.node_2_label},'
            f'\n\t{self.named_const_law}'
            f'{deformable_axial_with_named_law.element_footer()}'
        )
        self.assertEqual(str(deformable_axial_with_named_law), expected_str)

    def test_named_constitutive_law2(self):
        # Test with the second NamedConstitutiveLaw instance
        deformable_axial_with_named_law2 = l.DeformableAxial(
            idx=8,
            output='no',
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.named_const_law2
        )

        expected_str = (
            f'{deformable_axial_with_named_law2.element_header()}, deformable axial,\n\t{self.node_1_label},'
            f'\n\t{self.node_2_label},'
            f'\n\t{self.named_const_law2}'
            f'{deformable_axial_with_named_law2.element_footer()}'
        )
        self.assertEqual(str(deformable_axial_with_named_law2), expected_str)

    def test_warning_for_named_constitutive_law(self):
        # Test if a warning is issued when using a string for constitutive law
        with self.assertWarns(Warning):
            l.DeformableAxial(
                idx=9,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                const_law=l.NamedConstitutiveLaw("Some const law")
            )

class TestDeformableHinge2(unittest.TestCase):

    def setUp(self):
        self.node_1_label = 1
        self.node_2_label = 2

        # Optional values for testing with positions and orientations
        self.position_1 = l.Position2(
            relative_position=[1.0, 0.0, 0.0],
            reference=''
        )
        self.orientation_mat_1 = l.Position2(
            relative_position=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )
        self.position_2 = l.Position2(
            relative_position=[0.0, 1.0, 0.0],
            reference='global'
        )
        self.orientation_mat_2 = l.Position2(
            relative_position=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            reference='node'
        )

        # Example constitutive laws
        self.linear_elastic = l.LinearElastic(
            idx=1,
            law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW,
            stiffness=2000
        )
        
        self.named_const_law = l.NamedConstitutiveLaw("example_named_law")

        # Initialize DeformableHinge2 with required fields
        self.deformable_hinge2 = l.DeformableHinge2(
            idx=1,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.linear_elastic
        )

    def test_initialization(self):
        # Test that the DeformableHinge2 initializes correctly with provided values
        self.assertEqual(self.deformable_hinge2.node_1_label, self.node_1_label)
        self.assertEqual(self.deformable_hinge2.node_2_label, self.node_2_label)
        self.assertIsNone(self.deformable_hinge2.position_1)
        self.assertIsNone(self.deformable_hinge2.orientation_mat_1)
        self.assertIsNone(self.deformable_hinge2.position_2)
        self.assertIsNone(self.deformable_hinge2.orientation_mat_2)
        self.assertEqual(self.deformable_hinge2.const_law, self.linear_elastic)
        self.assertEqual(self.deformable_hinge2.idx, 1)

    def test_str_representation_without_optional(self):
        # Test the string output when optional positions and orientations are not provided
        expected_str = (
            f'{self.deformable_hinge2.element_header()}, deformable hinge'
            f',\n\t{self.node_1_label},'
            f'\n\t{self.node_2_label},'
            f'\n\t{self.linear_elastic}'
            f'{self.deformable_hinge2.element_footer()}'
        )
        self.assertEqual(str(self.deformable_hinge2), expected_str)

    def test_str_representation_with_optional(self):
        # Initialize with optional positions and orientations
        deformable_hinge2_with_optional = l.DeformableHinge2(
            idx=2,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2,
            const_law=self.named_const_law
        )

        expected_str = (
            f'{deformable_hinge2_with_optional.element_header()}, deformable hinge'
            f',\n\t{self.node_1_label},'
            f'\n\t\tposition, {self.position_1},'
            f'\n\t\torientation, {self.orientation_mat_1},'
            f'\n\t{self.node_2_label},'
            f'\n\t\tposition, {self.position_2},'
            f'\n\t\torientation, {self.orientation_mat_2},'
            f'\n\t{self.named_const_law}'
            f'{deformable_hinge2_with_optional.element_footer()}'
        )
        self.assertEqual(str(deformable_hinge2_with_optional), expected_str)

    def test_optional_none_handling(self):
        # Test to check that None is handled correctly for optional positions and orientations
        deformable_hinge2_without_optional = l.DeformableHinge2(
            idx=3,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            const_law=self.named_const_law
        )

        # Assert that optional parameters remain None
        self.assertIsNone(deformable_hinge2_without_optional.position_1)
        self.assertIsNone(deformable_hinge2_without_optional.orientation_mat_1)
        self.assertIsNone(deformable_hinge2_without_optional.position_2)
        self.assertIsNone(deformable_hinge2_without_optional.orientation_mat_2)

    def test_invalid_node_labels(self):
        # Test that passing invalid node labels raises the appropriate error
        with self.assertRaises(Exception):
            l.DeformableHinge2(
                idx=4,
                node_1_label="invalid_label",  # Invalid type for node_1_label
                node_2_label=self.node_2_label,
                const_law=self.linear_elastic
            )

        with self.assertRaises(Exception):
            l.DeformableHinge2(
                idx=5,
                node_1_label=self.node_1_label,
                node_2_label="invalid_label",  # Invalid type for node_2_label
                const_law=self.linear_elastic
            )

    def test_invalid_const_law(self):
        # Test that passing an invalid const_law raises the appropriate error
        with self.assertRaises(Exception):
            l.DeformableHinge2(
                idx=6,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                const_law="invalid_const_law"  # Invalid type for const_law
            )

class TestNamedConstitutiveLaw(unittest.TestCase):

    def test_string_input(self):
        with self.assertWarns(UserWarning) as cm:
            law = l.NamedConstitutiveLaw("linear elastic")
        self.assertEqual(str(law), "linear elastic")
        self.assertIn("Using a string for constitutive laws is not recommended and may be removed in the future.", str(cm.warning))

    def test_list_input(self):
        with self.assertWarns(UserWarning) as cm:
            law = l.NamedConstitutiveLaw(["linear elastic", "viscoelastic"])
        self.assertEqual(str(law), "linear elastic, viscoelastic")
        self.assertIn("Using a string for constitutive laws is not recommended and may be removed in the future.", str(cm.warning))

class TestDistance(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.distance_drive = l.ConstDriveCaller(const_value=5.0)

    def test_distance_creation_valid(self):
        # Test creating a Distance instance with valid data
        distance_joint = l.Distance(
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            distance=self.distance_drive,
            idx=10,
            output='yes'
        )
        self.assertIsInstance(distance_joint, l.Distance)
        self.assertEqual(distance_joint.node_1_label, self.node_1_label)
        self.assertEqual(str(distance_joint.distance), 'const, 5.0')

    def test_distance_creation_with_from_nodes(self):
        # Test creating a Distance instance with 'from nodes' as distance
        distance_joint = l.Distance(
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            distance='from nodes',
            idx=10
        )
        self.assertIsInstance(distance_joint, l.Distance)
        self.assertEqual(distance_joint.distance, 'from nodes')

    def test_distance_creation_missing_positions(self):
        # Test creating a Distance instance without positions
        distance_joint = l.Distance(
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            distance=self.distance_drive,
            idx=10
        )
        self.assertIsInstance(distance_joint, l.Distance)
        self.assertIsNone(distance_joint.position_1)
        self.assertIsNone(distance_joint.position_2)

    def test_distance_creation_invalid_distance_string(self):
        # Test creating a Distance instance with an invalid distance string
        with self.assertRaises(ValueError):
            l.Distance(
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                distance='invalid_string'
            )

    def test_distance_str_method(self):
        # Test the __str__ method of Distance
        distance_joint = l.Distance(
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            distance=self.distance_drive,
            idx=10,
            output='yes'
        )
        expected_str = (
            f'{distance_joint.element_header()}, distance'
            f',\n\t{self.node_1_label}, position, {self.position_1}'
            f',\n\t{self.node_2_label}, position, {self.position_2}'
            f',\n\t{self.distance_drive}'
            f'{distance_joint.element_footer()}'
        )
        self.maxDiff=None
        self.assertEqual(str(distance_joint), expected_str)

    def test_distance_output_option(self):
        # Test setting the output option to 'no'
        distance_joint = l.Distance(
            idx=10,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            distance=self.distance_drive,
            output='no'
        )
        self.assertEqual(distance_joint.output, 'no')
        self.assertIn(',\n\toutput, no', str(distance_joint))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_distance_missing_required_field(self):
        # Test creating a Distance instance missing a required field (distance)
        with self.assertRaises(Exception):
            l.Distance(
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label
                # Missing distance
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_distance_invalid_distance_type(self):
        # Test passing an invalid type for distance
        with self.assertRaises(Exception):
            l.Distance(
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                distance=123  # Invalid type, should be DriveCaller2 or 'from nodes'
            )
    # TODO: First check if the <MBVar> class definition has any errors 
    # def test_distance_with_mbvar_nodes(self):
    # # Test creating a Distance instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     distance_joint = l.Distance(
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2,
    #         distance=self.distance_drive
    #     )
    #     self.assertEqual(distance_joint.node_1_label, node_var_1)
    #     self.assertEqual(distance_joint.node_2_label, node_var_2)

class TestDriveDisplacement(unittest.TestCase):
    # TODO: Implement 'TplDriveCaller' first
    pass

class TestDriveDisplacementPin(unittest.TestCase):
    # TODO: Implement 'TplDriveCaller' first
    pass

class TestDriveHinge(unittest.TestCase):
    # TODO: Implement 'TplDriveCaller' first
    pass

class TestGimbalRotation(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_1_label = 1
        self.node_2_label = 2
        self.relative_orientation_mat_1 = l.Position2(relative_position=[0.0, 0.0, 1.0], reference='global')
        self.relative_orientation_mat_2 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='global')
        self.orientation_description = "euler123"

    def test_gimbal_rotation_creation_valid(self):
        # Test creating a GimbalRotation instance with valid data
        gimbal_rotation = l.GimbalRotation(
            node_1_label=self.node_1_label,
            relative_orientation_mat_1=self.relative_orientation_mat_1,
            node_2_label=self.node_2_label,
            relative_orientation_mat_2=self.relative_orientation_mat_2,
            orientation_description=self.orientation_description,
            idx=10,
            output='yes'
        )
        self.assertIsInstance(gimbal_rotation, l.GimbalRotation)
        self.assertEqual(gimbal_rotation.node_1_label, self.node_1_label)
        self.assertEqual(gimbal_rotation.orientation_description, self.orientation_description)

    def test_gimbal_rotation_creation_without_optional_fields(self):
        # Test creating a GimbalRotation instance without optional fields
        gimbal_rotation = l.GimbalRotation(
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            idx=5
        )
        self.assertIsInstance(gimbal_rotation, l.GimbalRotation)
        self.assertEqual(gimbal_rotation.node_1_label, self.node_1_label)
        self.assertIsNone(gimbal_rotation.relative_orientation_mat_1)
        self.assertIsNone(gimbal_rotation.orientation_description)

    def test_gimbal_rotation_invalid_orientation_description(self):
        # Test creating a GimbalRotation instance with invalid orientation_description
        with self.assertRaises(ValueError) as context:
            l.GimbalRotation(
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                orientation_description="invalid_description"
            )
        self.assertIn("Invalid orientation description", str(context.exception))

    def test_gimbal_rotation_str_method(self):
        # Test the __str__ method of GimbalRotation
        gimbal_rotation = l.GimbalRotation(
            node_1_label=self.node_1_label,
            relative_orientation_mat_1=self.relative_orientation_mat_1,
            node_2_label=self.node_2_label,
            relative_orientation_mat_2=self.relative_orientation_mat_2,
            orientation_description=self.orientation_description,
            idx=10
        )
        expected_str = (
            f'{gimbal_rotation.element_header()}, gimbal rotation'
            f',\n\t{self.node_1_label}'
            f', orientation, {self.relative_orientation_mat_1}'
            f',\n\t{self.node_2_label}'
            f', orientation, {self.relative_orientation_mat_2}'
            f',\n\torientation description, {self.orientation_description}'
            f'{gimbal_rotation.element_footer()}'
        )
        self.assertEqual(str(gimbal_rotation), expected_str)

    def test_gimbal_rotation_output_option(self):
        # Test setting the output option to 'no'
        gimbal_rotation = l.GimbalRotation(
            idx=10,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(gimbal_rotation.output, 'no')
        self.assertIn(',\n\toutput, no', str(gimbal_rotation))

    # TODO: First check if there are any errors in the <MBVar> class
    # def test_gimbal_rotation_with_mbvar_nodes(self):
    #     # Test creating a GimbalRotation instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     gimbal_rotation = l.GimbalRotation(
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2
    #     )
    #     self.assertEqual(gimbal_rotation.node_1_label, node_var_1)
    #     self.assertEqual(gimbal_rotation.node_2_label, node_var_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_gimbal_rotation_missing_required_field(self):
        # Test creating a GimbalRotation instance missing a required field (node_2_label)
        with self.assertRaises(Exception):
            l.GimbalRotation(
                node_1_label=self.node_1_label
                # Missing node_2_label
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_gimbal_rotation_invalid_relative_orientation_mat(self):
        # Test passing an invalid type for relative_orientation_mat_1
        with self.assertRaises(Exception):
            l.GimbalRotation(
                node_1_label=self.node_1_label,
                relative_orientation_mat_1=123,  # Invalid type, should be Position2
                node_2_label=self.node_2_label
            )

    def test_gimbal_rotation_orientation_description_none(self):
        # Test creating a GimbalRotation instance with orientation_description as None
        gimbal_rotation = l.GimbalRotation(
            idx=10,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            orientation_description=None
        )
        self.assertIsNone(gimbal_rotation.orientation_description)
        self.assertNotIn('orientation description', str(gimbal_rotation))

class TestImposedDisplacement(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.direction = [1.0, 0.0, 0.0]
        self.relative_position_drive = l.ConstDriveCaller(const_value=5.0)

    def test_imposed_displacement_creation_valid(self):
        # Test creating an ImposedDisplacement instance with valid data
        imposed_displacement = l.ImposedDisplacement(
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            direction=self.direction,
            relative_position=self.relative_position_drive,
            idx=10,
            output='yes'
        )
        self.assertIsInstance(imposed_displacement, l.ImposedDisplacement)
        self.assertEqual(imposed_displacement.node_1_label, self.node_1_label)
        self.assertEqual(imposed_displacement.direction, self.direction)
        self.assertEqual(str(imposed_displacement.relative_position), 'const, 5.0')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_missing_required_field(self):
        # Test creating an ImposedDisplacement instance missing a required field (direction)
        with self.assertRaises(Exception):
            l.ImposedDisplacement(
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                relative_position=self.relative_position_drive
                # Missing direction
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_invalid_direction(self):
        # Test passing an invalid type for direction
        with self.assertRaises(Exception):
            l.ImposedDisplacement(
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                direction=[1.0, 0.0],  # Invalid length, should be 3 elements
                relative_position=self.relative_position_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_invalid_relative_position(self):
        # Test passing an invalid type for relative_position
        with self.assertRaises(Exception):
            l.ImposedDisplacement(
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                direction=self.direction,
                relative_position=123  # Invalid type, should be DriveCaller
            )

    def test_imposed_displacement_str_method(self):
        # Test the __str__ method of ImposedDisplacement
        imposed_displacement = l.ImposedDisplacement(
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            direction=self.direction,
            relative_position=self.relative_position_drive,
            idx=10
        )
        expected_str = (
            f'{imposed_displacement.element_header()}, imposed displacement'
            f',\n\t{self.node_1_label}, {self.position_1}'
            f',\n\t{self.node_2_label}, {self.position_2}'
            f',\n\t{self.direction}'
            f',\n\t{self.relative_position_drive}'
            f'{imposed_displacement.element_footer()}'
        )
        self.maxDiff=None
        self.assertEqual(str(imposed_displacement), expected_str)

    def test_imposed_displacement_output_option(self):
        # Test setting the output option to 'no'
        imposed_displacement = l.ImposedDisplacement(
            idx=10,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            direction=self.direction,
            relative_position=self.relative_position_drive,
            output='no'
        )
        self.assertEqual(imposed_displacement.output, 'no')
        self.assertIn(',\n\toutput, no', str(imposed_displacement))

    # TODO: Check if <MBVar> class has any errors first
    # def test_imposed_displacement_with_mbvar_nodes(self):
    #     # Test creating an ImposedDisplacement instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     imposed_displacement = l.ImposedDisplacement(
    #         node_1_label=node_var_1,
    #         position_1=self.position_1,
    #         node_2_label=node_var_2,
    #         position_2=self.position_2,
    #         direction=self.direction,
    #         relative_position=self.relative_position_drive
    #     )
    #     self.assertEqual(imposed_displacement.node_1_label, node_var_1)
    #     self.assertEqual(imposed_displacement.node_2_label, node_var_2)

    # def test_imposed_displacement_direction_with_mbvars(self):
    #     # Test passing MBVar instances in the direction vector
    #     direction = [l.MBVar('dx', 'real', 1.0), l.MBVar('dy', 'real', 0.0), l.MBVar('dz', 'real', 0.0)]
    #     imposed_displacement = l.ImposedDisplacement(
    #         node_1_label=self.node_1_label,
    #         position_1=self.position_1,
    #         node_2_label=self.node_2_label,
    #         position_2=self.position_2,
    #         direction=direction,
    #         relative_position=self.relative_position_drive
    #     )
    #     self.assertEqual(imposed_displacement.direction, direction)

class TestImposedDisplacementPin(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_label = 1
        self.node_offset = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.offset = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.direction = [1.0, 0.0, 0.0]
        self.position_drive = l.ConstDriveCaller(const_value=5.0)
        self.position_drive_with_idx = l.ConstDriveCaller(const_value=5.0, idx=10)

    def test_imposed_displacement_pin_creation_valid(self):
        # Test creating an ImposedDisplacementPin instance with valid data
        imposed_displacement_pin = l.ImposedDisplacementPin(
            node_label=self.node_label,
            node_offset=self.node_offset,
            offset=self.offset,
            direction=self.direction,
            position=self.position_drive,
            idx=20,
            output='yes'
        )
        self.assertIsInstance(imposed_displacement_pin, l.ImposedDisplacementPin)
        self.assertEqual(imposed_displacement_pin.node_label, self.node_label)
        self.assertEqual(imposed_displacement_pin.direction, self.direction)
        self.assertEqual(str(imposed_displacement_pin.position), 'const, 5.0')

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_pin_missing_required_field(self):
        # Test creating an ImposedDisplacementPin instance missing a required field (direction)
        with self.assertRaises(Exception):
            l.ImposedDisplacementPin(
                idx=20,
                node_label=self.node_label,
                node_offset=self.node_offset,
                offset=self.offset,
                position=self.position_drive
                # Missing direction
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_pin_invalid_direction(self):
        # Test passing an invalid type for direction
        with self.assertRaises(Exception):
            l.ImposedDisplacementPin(
                node_label=self.node_label,
                node_offset=self.node_offset,
                offset=self.offset,
                direction=[1.0, 0.0],  # Invalid length, should be 3 elements
                position=self.position_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_pin_invalid_position(self):
        # Test passing an invalid type for position
        with self.assertRaises(Exception):
            l.ImposedDisplacementPin(
                idx=20,
                node_label=self.node_label,
                node_offset=self.node_offset,
                offset=self.offset,
                direction=self.direction,
                position=123  # Invalid type, should be DriveCaller
            )

    def test_imposed_displacement_pin_str_method(self):
        # Test the __str__ method of ImposedDisplacementPin without position idx
        imposed_displacement_pin = l.ImposedDisplacementPin(
            node_label=self.node_label,
            node_offset=self.node_offset,
            offset=self.offset,
            direction=self.direction,
            position=self.position_drive,
            idx=20
        )
        expected_str = (
            f'{imposed_displacement_pin.element_header()}, imposed displacement pin'
            f',\n\t{self.node_label}, {self.node_offset}'
            f',\n\t{self.offset}'
            f',\n\t{self.direction}'
            f',\n\t{self.position_drive}'
            f'{imposed_displacement_pin.element_footer()}'
        )
        self.assertEqual(str(imposed_displacement_pin), expected_str)

    def test_imposed_displacement_pin_str_method_with_position_idx(self):
        # Test the __str__ method of ImposedDisplacementPin with position idx
        imposed_displacement_pin = l.ImposedDisplacementPin(
            node_label=self.node_label,
            node_offset=self.node_offset,
            offset=self.offset,
            direction=self.direction,
            position=self.position_drive_with_idx,
            idx=20
        )
        expected_str = (
            f'{imposed_displacement_pin.element_header()}, imposed displacement pin'
            f',\n\t{self.node_label}, {self.node_offset}'
            f',\n\t{self.offset}'
            f',\n\t{self.direction}'
            f',\n\treference, {self.position_drive_with_idx.idx}'
            f'{imposed_displacement_pin.element_footer()}'
        )
        self.assertEqual(str(imposed_displacement_pin), expected_str)

    def test_imposed_displacement_pin_output_option(self):
        # Test setting the output option to 'no'
        imposed_displacement_pin = l.ImposedDisplacementPin(
            idx=20,
            node_label=self.node_label,
            node_offset=self.node_offset,
            offset=self.offset,
            direction=self.direction,
            position=self.position_drive,
            output='no'
        )
        self.assertEqual(imposed_displacement_pin.output, 'no')
        self.assertIn(',\n\toutput, no', str(imposed_displacement_pin))

    # TODO: First check if MBVar class has any errors
    # def test_imposed_displacement_pin_with_mbvar_node_label(self):
    #     # Test creating an ImposedDisplacementPin instance with MBVar as node_label
    #     node_var = MBVar(name='node_var', var_type='integer', expression=100)
    #     imposed_displacement_pin = ImposedDisplacementPin(
    #         node_label=node_var,
    #         node_offset=self.node_offset,
    #         offset=self.offset,
    #         direction=self.direction,
    #         position=self.position_drive
    #     )
    #     self.assertEqual(imposed_displacement_pin.node_label, node_var)

    def test_imposed_displacement_pin_invalid_direction_values(self):
        # Test passing a direction vector that is not a unit vector
        with self.assertRaises(ValueError):
            l.ImposedDisplacementPin(
                node_label=self.node_label,
                node_offset=self.node_offset,
                offset=self.offset,
                direction=[2.0, 0.0, 0.0],  # Not a unit vector
                position=self.position_drive
            )

    # TODO: First check if MBVar class has any errors
    # def test_imposed_displacement_pin_direction_with_mbvars(self):
    #     # Test passing MBVar instances in the direction vector
    #     direction = [MBVar('dx', 'real', 1.0), MBVar('dy', 'real', 0.0), MBVar('dz', 'real', 0.0)]
    #     imposed_displacement_pin = ImposedDisplacementPin(
    #         node_label=self.node_label,
    #         node_offset=self.node_offset,
    #         offset=self.offset,
    #         direction=direction,
    #         position=self.position_drive
    #     )
    #     self.assertEqual(imposed_displacement_pin.direction, direction)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_pin_missing_node_offset(self):
        # Test creating an ImposedDisplacementPin instance missing node_offset
        with self.assertRaises(Exception):
            l.ImposedDisplacementPin(
                idx=20,
                node_label=self.node_label,
                offset=self.offset,
                direction=self.direction,
                position=self.position_drive
                # Missing node_offset
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_imposed_displacement_pin_missing_offset(self):
        # Test creating an ImposedDisplacementPin instance missing offset
        with self.assertRaises(Exception):
            l.ImposedDisplacementPin(
                idx=20,
                node_label=self.node_label,
                node_offset=self.node_offset,
                direction=self.direction,
                position=self.position_drive
                # Missing offset
            )

class TestInLine(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_1_label = 1
        self.node_2_label = 2
        self.position = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.orientation = l.Position2(relative_position=[0.0, 0.0, 1.0], reference='global')
        self.offset = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.idx = 10

    def test_inline_creation_valid(self):
        # Test creating an InLine instance with all valid data
        inline_joint = l.InLine(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position=self.position,
            orientation=self.orientation,
            node_2_label=self.node_2_label,
            offset=self.offset
        )
        self.assertIsInstance(inline_joint, l.InLine)
        self.assertEqual(inline_joint.node_1_label, self.node_1_label)
        self.assertEqual(inline_joint.position, self.position)
        self.assertEqual(inline_joint.orientation, self.orientation)
        self.assertEqual(inline_joint.node_2_label, self.node_2_label)
        self.assertEqual(inline_joint.offset, self.offset)

    def test_inline_creation_without_optional_fields(self):
        # Test creating an InLine instance without optional fields
        inline_joint = l.InLine(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )
        self.assertIsInstance(inline_joint, l.InLine)
        self.assertEqual(inline_joint.node_1_label, self.node_1_label)
        self.assertIsNone(inline_joint.position)
        self.assertIsNone(inline_joint.orientation)
        self.assertEqual(inline_joint.node_2_label, self.node_2_label)
        self.assertIsNone(inline_joint.offset)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inline_invalid_node_label_type(self):
        # Test creating an InLine instance with invalid node_1_label type
        with self.assertRaises(Exception):
            l.InLine(
                idx=self.idx,
                node_1_label="invalid_node_label",  # Should be int or MBVar
                node_2_label=self.node_2_label
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inline_missing_required_fields(self):
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.InLine(
                idx=self.idx,
                node_2_label=self.node_2_label
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.InLine(
                idx=self.idx,
                node_1_label=self.node_1_label
            )

    def test_inline_str_method(self):
        # Test the __str__ method of InLine
        inline_joint = l.InLine(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position=self.position,
            orientation=self.orientation,
            node_2_label=self.node_2_label,
            offset=self.offset
        )
        expected_str = (
            f'{inline_joint.element_header()}, in line'
            f',\n\t{self.node_1_label}'
            f', position, {self.position}'
            f'\n\t, orientation, {self.orientation}'
            f',\n\t{self.node_2_label}'
            f', offset, {self.offset}'
            f'{inline_joint.element_footer()}'
        )
        self.assertEqual(str(inline_joint), expected_str)

    def test_inline_output_option(self):
        # Test setting the output option to 'no'
        inline_joint = l.InLine(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(inline_joint.output, 'no')
        self.assertIn(',\n\toutput, no', str(inline_joint))

    # TODO: Ensure MBVar class works correctly before running this test
    # def test_inline_with_mbvar_node_labels(self):
    #     # Test creating an InLine instance with MBVar as node labels
    #     # Ensure MBVar class works correctly before running this test
    #     try:
    #         node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #         node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #         inline_joint = l.InLine(
    #             idx=self.idx,
    #             node_1_label=node_var_1,
    #             node_2_label=node_var_2,
    #             position=self.position,
    #             orientation=self.orientation,
    #             offset=self.offset
    #         )
    #         self.assertEqual(inline_joint.node_1_label, node_var_1)
    #         self.assertEqual(inline_joint.node_2_label, node_var_2)
    #     except Exception as e:
    #         self.skipTest(f"MBVar class has errors: {e}")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inline_invalid_position_type(self):
        # Test creating an InLine instance with invalid position type
        with self.assertRaises(Exception):
            l.InLine(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position="invalid_position",  # Should be Position2 or None
                node_2_label=self.node_2_label
            )

    def test_inline_missing_idx(self):
        # Test creating an InLine instance without idx
        with self.assertRaises(Exception):
            l.InLine(
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label
            )

class TestInPlane(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_1_label = 1
        self.node_2_label = 2
        self.position = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.offset = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.relative_direction_unit = [1.0, 0.0, 0.0]  # Unit vector
        self.relative_direction_non_unit = [2.0, 0.0, 0.0]  # Non-unit vector
        self.idx = 10

    def test_inplane_creation_valid(self):
        # Test creating an InPlane instance with all valid data
        inplane_joint = l.InPlane(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position=self.position,
            relative_direction=self.relative_direction_unit,
            node_2_label=self.node_2_label,
            offset=self.offset
        )
        self.assertIsInstance(inplane_joint, l.InPlane)
        self.assertEqual(inplane_joint.node_1_label, self.node_1_label)
        self.assertEqual(inplane_joint.position, self.position)
        self.assertEqual(inplane_joint.relative_direction, self.relative_direction_unit)
        self.assertEqual(inplane_joint.node_2_label, self.node_2_label)
        self.assertEqual(inplane_joint.offset, self.offset)

    def test_inplane_creation_without_optional_fields(self):
        # Test creating an InPlane instance without optional fields
        inplane_joint = l.InPlane(
            idx=self.idx,
            node_1_label=self.node_1_label,
            relative_direction=self.relative_direction_unit,
            node_2_label=self.node_2_label
        )
        self.assertIsInstance(inplane_joint, l.InPlane)
        self.assertEqual(inplane_joint.node_1_label, self.node_1_label)
        self.assertIsNone(inplane_joint.position)
        self.assertEqual(inplane_joint.relative_direction, self.relative_direction_unit)
        self.assertEqual(inplane_joint.node_2_label, self.node_2_label)
        self.assertIsNone(inplane_joint.offset)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_invalid_relative_direction_unit_vector(self):
        # Test creating an InPlane instance with non-unit relative_direction
        with self.assertRaises(Exception) as context:
            l.InPlane(
                idx=self.idx,
                node_1_label=self.node_1_label,
                relative_direction=self.relative_direction_non_unit,
                node_2_label=self.node_2_label
            )
        self.assertIn("relative_direction must be a unit vector", str(context.exception))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_invalid_relative_direction_length(self):
        # Test creating an InPlane instance with invalid relative_direction length
        with self.assertRaises(Exception) as context:
            l.InPlane(
                idx=self.idx,
                node_1_label=self.node_1_label,
                relative_direction=[1.0, 0.0],  # Invalid length
                node_2_label=self.node_2_label
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_invalid_node_label_type(self):
        # Test creating an InPlane instance with invalid node_1_label type
        with self.assertRaises(Exception):
            l.InPlane(
                idx=self.idx,
                node_1_label="invalid_node_label",  # Should be int or MBVar
                relative_direction=self.relative_direction_unit,
                node_2_label=self.node_2_label
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_missing_required_fields(self):
        # Missing relative_direction
        with self.assertRaises(Exception) as context:
            l.InPlane(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label
            )

    def test_inplane_str_method(self):
        # Test the __str__ method of InPlane
        inplane_joint = l.InPlane(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position=self.position,
            relative_direction=self.relative_direction_unit,
            node_2_label=self.node_2_label,
            offset=self.offset
        )
        expected_str = (
            f'{inplane_joint.element_header()}, in plane'
            f',\n\t{self.node_1_label}'
            f', position, {self.position}'
            f',\n\t{self.relative_direction_unit}'
            f',\n\t{self.node_2_label}'
            f', offset, {self.offset}'
            f'{inplane_joint.element_footer()}'
        )
        self.assertEqual(str(inplane_joint), expected_str)

    def test_inplane_output_option(self):
        # Test setting the output option to 'no'
        inplane_joint = l.InPlane(
            idx=self.idx,
            node_1_label=self.node_1_label,
            relative_direction=self.relative_direction_unit,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(inplane_joint.output, 'no')
        self.assertIn(',\n\toutput, no', str(inplane_joint))

    # # TODO: Ensure MBVar class works correctly before running this test
    # def test_inplane_with_mbvar_node_labels(self):
    #     # Test creating an InPlane instance with MBVar as node labels
    #     try:
    #         node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #         node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #         inplane_joint = l.InPlane(
    #             idx=self.idx,
    #             node_1_label=node_var_1,
    #             relative_direction=self.relative_direction_unit,
    #             node_2_label=node_var_2,
    #             position=self.position,
    #             offset=self.offset
    #         )
    #         self.assertEqual(inplane_joint.node_1_label, node_var_1)
    #         self.assertEqual(inplane_joint.node_2_label, node_var_2)
    #     except Exception as e:
    #         self.skipTest(f"MBVar class has errors: {e}")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_invalid_position_type(self):
        # Test creating an InPlane instance with invalid position type
        with self.assertRaises(Exception):
            l.InPlane(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position="invalid_position",  # Should be Position2 or None
                relative_direction=self.relative_direction_unit,
                node_2_label=self.node_2_label
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_inplane_missing_idx(self):
        # Test creating an InPlane instance without idx
        with self.assertRaises(Exception):
            l.InPlane(
                node_1_label=self.node_1_label,
                relative_direction=self.relative_direction_unit,
                node_2_label=self.node_2_label
            )

class TestLinearAcceleration(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_label = 1
        self.relative_direction_unit = [1.0, 0.0, 0.0]  # Unit vector
        self.relative_direction_non_unit = [2.0, 0.0, 0.0]  # Non-unit vector
        self.acceleration_drive = l.ConstDriveCaller(const_value=5.0)
        self.acceleration_drive_with_idx = l.ConstDriveCaller(const_value=5.0, idx=10)
        self.idx = 10

    def test_linear_acceleration_creation_valid(self):
        # Test creating a LinearAcceleration instance with all valid data
        linear_acceleration_joint = l.LinearAcceleration(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            acceleration=self.acceleration_drive
        )
        self.assertIsInstance(linear_acceleration_joint, l.LinearAcceleration)
        self.assertEqual(linear_acceleration_joint.node_label, self.node_label)
        self.assertEqual(linear_acceleration_joint.relative_direction, self.relative_direction_unit)
        self.assertEqual(linear_acceleration_joint.acceleration, self.acceleration_drive)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_acceleration_invalid_relative_direction_unit_vector(self):
        # Test creating a LinearAcceleration instance with non-unit relative_direction
        with self.assertRaises(Exception) as context:
            l.LinearAcceleration(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_non_unit,
                acceleration=self.acceleration_drive
            )
        self.assertIn("relative_direction must be a unit vector", str(context.exception))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_acceleration_invalid_relative_direction_length(self):
        # Test creating a LinearAcceleration instance with invalid relative_direction length
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=[1.0, 0.0],  # Invalid length
                acceleration=self.acceleration_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_acceleration_missing_required_fields(self):
        # Missing relative_direction
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                idx=self.idx,
                node_label=self.node_label,
                acceleration=self.acceleration_drive
            )
        # Missing acceleration
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit
            )
        # Missing node_label
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                idx=self.idx,
                relative_direction=self.relative_direction_unit,
                acceleration=self.acceleration_drive
            )

    def test_linear_acceleration_str_method(self):
        # Test the __str__ method of LinearAcceleration
        linear_acceleration_joint = l.LinearAcceleration(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            acceleration=self.acceleration_drive
        )
        expected_str = (
            f'{linear_acceleration_joint.element_header()}, linear acceleration'
            f',\n\t{self.node_label}'
            f',\n\t {self.relative_direction_unit}'
            f',\n\t{self.acceleration_drive}'
            f'{linear_acceleration_joint.element_footer()}'
        )
        self.assertEqual(str(linear_acceleration_joint), expected_str)

    def test_linear_acceleration_str_method_with_acceleration_idx(self):
        # Test the __str__ method when acceleration.idx is provided and non-negative
        linear_acceleration_joint = l.LinearAcceleration(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            acceleration=self.acceleration_drive_with_idx
        )
        expected_str = (
            f'{linear_acceleration_joint.element_header()}, linear acceleration'
            f',\n\t{self.node_label}'
            f',\n\t {self.relative_direction_unit}'
            f',\n\treference, {self.acceleration_drive_with_idx.idx}'
            f'{linear_acceleration_joint.element_footer()}'
        )
        self.assertEqual(str(linear_acceleration_joint), expected_str)

    def test_linear_acceleration_output_option(self):
        # Test setting the output option to 'no'
        linear_acceleration_joint = l.LinearAcceleration(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            acceleration=self.acceleration_drive,
            output='no'
        )
        self.assertEqual(linear_acceleration_joint.output, 'no')
        self.assertIn(',\n\toutput, no', str(linear_acceleration_joint))

    # TODO: First check if MBVar class has any errors
    # def test_linear_acceleration_with_mbvar_node_label(self):
    #     # Test creating a LinearAcceleration instance with MBVar as node_label
    #     try:
    #         node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #         linear_acceleration_joint = l.LinearAcceleration(
    #             idx=self.idx,
    #             node_label=node_var,
    #             relative_direction=self.relative_direction_unit,
    #             acceleration=self.acceleration_drive
    #         )
    #         self.assertEqual(linear_acceleration_joint.node_label, node_var)
    #     except Exception as e:
    #         self.skipTest(f"MBVar class has errors: {e}")

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_acceleration_invalid_acceleration_type(self):
        # Test creating a LinearAcceleration instance with invalid acceleration type
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit,
                acceleration=123  # Invalid type, should be DriveCaller or DriveCaller2
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_acceleration_missing_idx(self):
        # Test creating a LinearAcceleration instance without idx
        with self.assertRaises(Exception):
            l.LinearAcceleration(
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit,
                acceleration=self.acceleration_drive
            )

class TestLinearVelocity(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.node_label = 1
        self.relative_direction_unit = [1.0, 0.0, 0.0]  # Unit vector
        self.relative_direction_non_unit = [2.0, 0.0, 0.0]  # Non-unit vector
        self.velocity_drive = l.ConstDriveCaller(const_value=5.0)
        self.velocity_drive_with_idx = l.ConstDriveCaller(const_value=5.0, idx=10)
        self.idx = 10

    def test_linear_velocity_creation_valid(self):
        # Test creating a LinearVelocity instance with all valid data
        linear_velocity_joint = l.LinearVelocity(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            velocity=self.velocity_drive
        )
        self.assertIsInstance(linear_velocity_joint, l.LinearVelocity)
        self.assertEqual(linear_velocity_joint.node_label, self.node_label)
        self.assertEqual(linear_velocity_joint.relative_direction, self.relative_direction_unit)
        self.assertEqual(linear_velocity_joint.velocity, self.velocity_drive)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_velocity_invalid_relative_direction_unit_vector(self):
        # Test creating a LinearVelocity instance with non-unit relative_direction
        with self.assertRaises(Exception) as context:
            l.LinearVelocity(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_non_unit,
                velocity=self.velocity_drive
            )
        self.assertIn("relative_direction must be a unit vector", str(context.exception))

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_velocity_invalid_relative_direction_length(self):
        # Test creating a LinearVelocity instance with invalid relative_direction length
        with self.assertRaises(Exception):
            l.LinearVelocity(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=[1.0, 0.0],  # Invalid length
                velocity=self.velocity_drive
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_velocity_missing_required_fields(self):
        # Missing relative_direction
        with self.assertRaises(Exception):
            l.LinearVelocity(
                idx=self.idx,
                node_label=self.node_label,
                velocity=self.velocity_drive
            )
        # Missing velocity
        with self.assertRaises(Exception):
            l.LinearVelocity(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit
            )
        # Missing node_label
        with self.assertRaises(Exception):
            l.LinearVelocity(
                idx=self.idx,
                relative_direction=self.relative_direction_unit,
                velocity=self.velocity_drive
            )

    def test_linear_velocity_str_method(self):
        # Test the __str__ method of LinearVelocity
        linear_velocity_joint = l.LinearVelocity(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            velocity=self.velocity_drive
        )
        expected_str = (
            f'{linear_velocity_joint.element_header()}, linear velocity'
            f',\n\t{self.node_label}'
            f',\n\t {self.relative_direction_unit}'
            f',\n\t{self.velocity_drive}'
            f'{linear_velocity_joint.element_footer()}'
        )
        self.assertEqual(str(linear_velocity_joint), expected_str)

    def test_linear_velocity_str_method_with_velocity_idx(self):
        # Test the __str__ method when velocity.idx is provided and non-negative
        linear_velocity_joint = l.LinearVelocity(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            velocity=self.velocity_drive_with_idx
        )
        expected_str = (
            f'{linear_velocity_joint.element_header()}, linear velocity'
            f',\n\t{self.node_label}'
            f',\n\t {self.relative_direction_unit}'
            f',\n\treference, {self.velocity_drive_with_idx.idx}'
            f'{linear_velocity_joint.element_footer()}'
        )
        self.assertEqual(str(linear_velocity_joint), expected_str)

    def test_linear_velocity_output_option(self):
        # Test setting the output option to 'no'
        linear_velocity_joint = l.LinearVelocity(
            idx=self.idx,
            node_label=self.node_label,
            relative_direction=self.relative_direction_unit,
            velocity=self.velocity_drive,
            output='no'
        )
        self.assertEqual(linear_velocity_joint.output, 'no')
        self.assertIn(',\n\toutput, no', str(linear_velocity_joint))

    # # TODO: First check if MBVar class has any errors
    # def test_linear_velocity_with_mbvar_node_label(self):
    #     # Test creating a LinearVelocity instance with MBVar as node_label
    #     node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #     linear_velocity_joint = l.LinearVelocity(
    #         idx=self.idx,
    #         node_label=node_var,
    #         relative_direction=self.relative_direction_unit,
    #         velocity=self.velocity_drive
    #     )
    #     self.assertEqual(linear_velocity_joint.node_label, node_var)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_velocity_invalid_velocity_type(self):
        # Test creating a LinearVelocity instance with invalid velocity type
        with self.assertRaises(Exception):
            l.LinearVelocity(
                idx=self.idx,
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit,
                velocity=123  # Invalid type, should be DriveCaller or DriveCaller2
            )

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_linear_velocity_missing_idx(self):
        # Test creating a LinearVelocity instance without idx
        with self.assertRaises(Exception):
            l.LinearVelocity(
                node_label=self.node_label,
                relative_direction=self.relative_direction_unit,
                velocity=self.velocity_drive
            )

class TestPlaneDisplacement(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 10
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.orientation_mat_1 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='node')
        self.orientation_mat_2 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='node')

    def test_plane_displacement_creation_valid(self):
        # Test creating a PlaneDisplacement instance with all valid data
        plane_displacement = l.PlaneDisplacement(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        self.assertIsInstance(plane_displacement, l.PlaneDisplacement)
        self.assertEqual(plane_displacement.node_1_label, self.node_1_label)
        self.assertEqual(plane_displacement.position_1, self.position_1)
        self.assertEqual(plane_displacement.orientation_mat_1, self.orientation_mat_1)
        self.assertEqual(plane_displacement.node_2_label, self.node_2_label)
        self.assertEqual(plane_displacement.position_2, self.position_2)
        self.assertEqual(plane_displacement.orientation_mat_2, self.orientation_mat_2)

    def test_plane_displacement_creation_without_optional_fields(self):
        # Test creating a PlaneDisplacement instance without optional orientation matrices
        plane_displacement = l.PlaneDisplacement(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2
        )
        self.assertIsInstance(plane_displacement, l.PlaneDisplacement)
        self.assertIsNone(plane_displacement.orientation_mat_1)
        self.assertIsNone(plane_displacement.orientation_mat_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_plane_displacement_missing_required_fields(self):
        # Missing position_1
        with self.assertRaises(Exception):
            l.PlaneDisplacement(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                position_2=self.position_2
            )
        # Missing position_2
        with self.assertRaises(Exception):
            l.PlaneDisplacement(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label
            )

    def test_plane_displacement_str_method(self):
        # Test the __str__ method of PlaneDisplacement
        plane_displacement = l.PlaneDisplacement(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        expected_str = (
            f'{plane_displacement.element_header()}, plane displacement'
            f',\n\t{self.node_1_label}, position, {self.position_1}'
            f',\n\torientation, {self.orientation_mat_1}'
            f',\n\t{self.node_2_label}, position, {self.position_2}'
            f',\n\torientation, {self.orientation_mat_2}'
            f'{plane_displacement.element_footer()}'
        )
        self.assertEqual(str(plane_displacement), expected_str)

    def test_plane_displacement_output_option(self):
        # Test setting the output option to 'no'
        plane_displacement = l.PlaneDisplacement(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            output='no'
        )
        self.assertEqual(plane_displacement.output, 'no')
        self.assertIn(',\n\toutput, no', str(plane_displacement))

    # #TODO: First check if the MBVar class has any errors
    # def test_plane_displacement_with_mbvar_node_labels(self):
    #     # Test creating a PlaneDisplacement instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     plane_displacement = l.PlaneDisplacement(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         position_1=self.position_1,
    #         node_2_label=node_var_2,
    #         position_2=self.position_2
    #     )
    #     self.assertEqual(plane_displacement.node_1_label, node_var_1)
    #     self.assertEqual(plane_displacement.node_2_label, node_var_2)

class TestPlaneDisplacementPin(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 20
        self.node_label = 1
        self.relative_offset = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.absolute_pin_position = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.relative_orientation_mat = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='node')
        self.absolute_pin_orientation_mat = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='global')

    def test_plane_displacement_pin_creation_valid(self):
        # Test creating a PlaneDisplacementPin instance with all valid data
        plane_displacement_pin = l.PlaneDisplacementPin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            relative_orientation_mat=self.relative_orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_pin_orientation_mat=self.absolute_pin_orientation_mat
        )
        self.assertIsInstance(plane_displacement_pin, l.PlaneDisplacementPin)
        self.assertEqual(plane_displacement_pin.node_label, self.node_label)
        self.assertEqual(plane_displacement_pin.relative_offset, self.relative_offset)
        self.assertEqual(plane_displacement_pin.relative_orientation_mat, self.relative_orientation_mat)
        self.assertEqual(plane_displacement_pin.absolute_pin_position, self.absolute_pin_position)
        self.assertEqual(plane_displacement_pin.absolute_pin_orientation_mat, self.absolute_pin_orientation_mat)

    def test_plane_displacement_pin_creation_without_optional_fields(self):
        # Test creating a PlaneDisplacementPin instance without optional orientation matrices
        plane_displacement_pin = l.PlaneDisplacementPin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            absolute_pin_position=self.absolute_pin_position
        )
        self.assertIsInstance(plane_displacement_pin, l.PlaneDisplacementPin)
        self.assertIsNone(plane_displacement_pin.relative_orientation_mat)
        self.assertIsNone(plane_displacement_pin.absolute_pin_orientation_mat)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_plane_displacement_pin_missing_required_fields(self):
        # Missing relative_offset
        with self.assertRaises(Exception):
            l.PlaneDisplacementPin(
                idx=self.idx,
                node_label=self.node_label,
                absolute_pin_position=self.absolute_pin_position
            )
        # Missing absolute_pin_position
        with self.assertRaises(Exception):
            l.PlaneDisplacementPin(
                idx=self.idx,
                node_label=self.node_label,
                relative_offset=self.relative_offset
            )

    def test_plane_displacement_pin_str_method(self):
        # Test the __str__ method of PlaneDisplacementPin
        plane_displacement_pin = l.PlaneDisplacementPin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            relative_orientation_mat=self.relative_orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_pin_orientation_mat=self.absolute_pin_orientation_mat
        )
        expected_str = (
            f'{plane_displacement_pin.element_header()}, plane displacement pin'
            f',\n\t{self.node_label}'
            f',\n\t\tposition, {self.relative_offset}'
            f',\n\t\torientation, {self.relative_orientation_mat}'
            f',\n\tposition, {self.absolute_pin_position}'
            f',\n\torientation, {self.absolute_pin_orientation_mat}'
            f'{plane_displacement_pin.element_footer()}'
        )
        self.assertEqual(str(plane_displacement_pin), expected_str)

    def test_plane_displacement_pin_output_option(self):
        # Test setting the output option to 'no'
        plane_displacement_pin = l.PlaneDisplacementPin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            absolute_pin_position=self.absolute_pin_position,
            output='no'
        )
        self.assertEqual(plane_displacement_pin.output, 'no')
        self.assertIn(',\n\toutput, no', str(plane_displacement_pin))

    # TODO: First check if MBVar class has any errors
    # def test_plane_displacement_pin_with_mbvar_node_label(self):
    #     # Test creating a PlaneDisplacementPin instance with MBVar as node_label
    #     node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #     plane_displacement_pin = l.PlaneDisplacementPin(
    #         idx=self.idx,
    #         node_label=node_var,
    #         relative_offset=self.relative_offset,
    #         absolute_pin_position=self.absolute_pin_position
    #     )
    #     self.assertEqual(plane_displacement_pin.node_label, node_var)

class TestPrismatic(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 30
        self.node_1_label = 1
        self.node_2_label = 2
        self.relative_orientation_mat_1 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='other node')
        self.relative_orientation_mat_2 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='other node')

    def test_prismatic_creation_valid(self):
        # Test creating a Prismatic instance with all valid data
        prismatic = l.Prismatic(
            idx=self.idx,
            node_1_label=self.node_1_label,
            relative_orientation_mat_1=self.relative_orientation_mat_1,
            node_2_label=self.node_2_label,
            relative_orientation_mat_2=self.relative_orientation_mat_2
        )
        self.assertIsInstance(prismatic, l.Prismatic)
        self.assertEqual(prismatic.node_1_label, self.node_1_label)
        self.assertEqual(prismatic.relative_orientation_mat_1, self.relative_orientation_mat_1)
        self.assertEqual(prismatic.node_2_label, self.node_2_label)
        self.assertEqual(prismatic.relative_orientation_mat_2, self.relative_orientation_mat_2)

    def test_prismatic_creation_without_optional_fields(self):
        # Test creating a Prismatic instance without optional orientation matrices
        prismatic = l.Prismatic(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )
        self.assertIsInstance(prismatic, l.Prismatic)
        self.assertIsNone(prismatic.relative_orientation_mat_1)
        self.assertIsNone(prismatic.relative_orientation_mat_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_prismatic_missing_required_fields(self):
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.Prismatic(
                idx=self.idx,
                node_2_label=self.node_2_label
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.Prismatic(
                idx=self.idx,
                node_1_label=self.node_1_label
            )

    def test_prismatic_str_method(self):
        # Test the __str__ method of Prismatic
        prismatic = l.Prismatic(
            idx=self.idx,
            node_1_label=self.node_1_label,
            relative_orientation_mat_1=self.relative_orientation_mat_1,
            node_2_label=self.node_2_label,
            relative_orientation_mat_2=self.relative_orientation_mat_2
        )
        expected_str = (
            f'{prismatic.element_header()}, prismatic'
            f',\n\t{self.node_1_label}, orientation, {self.relative_orientation_mat_1}'
            f',\n\t{self.node_2_label}, orientation, {self.relative_orientation_mat_2}'
            f'{prismatic.element_footer()}'
        )
        self.assertEqual(str(prismatic), expected_str)

    def test_prismatic_output_option(self):
        # Test setting the output option to 'no'
        prismatic = l.Prismatic(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(prismatic.output, 'no')
        self.assertIn(',\n\toutput, no', str(prismatic))

    # # TODO: Test if MBVar class has any errors 
    # def test_prismatic_with_mbvar_node_labels(self):
    #     # Test creating a Prismatic instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     prismatic = l.Prismatic(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2
    #     )
    #     self.assertEqual(prismatic.node_1_label, node_var_1)
    #     self.assertEqual(prismatic.node_2_label, node_var_2)

class TestRevoluteHinge(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 10
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='node')
        self.orientation_mat_1 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='other node')
        self.orientation_mat_2 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='other node')
        self.initial_theta = 0.0
        self.friction = 0.5
        self.preload = 100.0
        self.friction_model = 'coulomb'
        self.shape_function = 'linear'

    def test_revolute_hinge_creation_valid(self):
        # Test creating a RevoluteHinge instance with all valid data, including friction parameters
        revolute_hinge = l.RevoluteHinge(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2,
            initial_theta=self.initial_theta,
            friction=self.friction,
            preload=self.preload,
            friction_model=self.friction_model,
            shape_function=self.shape_function
        )
        self.assertIsInstance(revolute_hinge, l.RevoluteHinge)
        self.assertEqual(revolute_hinge.node_1_label, self.node_1_label)
        self.assertEqual(revolute_hinge.position_1, self.position_1)
        self.assertEqual(revolute_hinge.orientation_mat_1, self.orientation_mat_1)
        self.assertEqual(revolute_hinge.node_2_label, self.node_2_label)
        self.assertEqual(revolute_hinge.position_2, self.position_2)
        self.assertEqual(revolute_hinge.orientation_mat_2, self.orientation_mat_2)
        self.assertEqual(revolute_hinge.initial_theta, self.initial_theta)
        self.assertEqual(revolute_hinge.friction, self.friction)
        self.assertEqual(revolute_hinge.preload, self.preload)
        self.assertEqual(revolute_hinge.friction_model, self.friction_model)
        self.assertEqual(revolute_hinge.shape_function, self.shape_function)

    def test_revolute_hinge_creation_without_optional_fields(self):
        # Test creating a RevoluteHinge instance without optional fields
        revolute_hinge = l.RevoluteHinge(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2
        )
        self.assertIsInstance(revolute_hinge, l.RevoluteHinge)
        self.assertIsNone(revolute_hinge.orientation_mat_1)
        self.assertIsNone(revolute_hinge.orientation_mat_2)
        self.assertIsNone(revolute_hinge.initial_theta)
        self.assertIsNone(revolute_hinge.friction)
        self.assertIsNone(revolute_hinge.preload)
        self.assertIsNone(revolute_hinge.friction_model)
        self.assertIsNone(revolute_hinge.shape_function)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_revolute_hinge_missing_required_fields(self):
        # Missing position_1
        with self.assertRaises(Exception):
            l.RevoluteHinge(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                position_2=self.position_2
            )
        # Missing position_2
        with self.assertRaises(Exception):
            l.RevoluteHinge(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label
            )

    def test_revolute_hinge_str_method(self):
        # Test the __str__ method of RevoluteHinge
        revolute_hinge = l.RevoluteHinge(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2,
            initial_theta=self.initial_theta,
            friction=self.friction,
            preload=self.preload,
            friction_model=self.friction_model,
            shape_function=self.shape_function
        )
        expected_str = (
            f'{revolute_hinge.element_header()}, revolute hinge'
            f',\n\t{self.node_1_label}'
            f',\n\t\tposition, {self.position_1}'
            f',\n\t\torientation, {self.orientation_mat_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\tposition, {self.position_2}'
            f',\n\t\torientation, {self.orientation_mat_2}'
            f',\n\tinitial theta, {self.initial_theta}'
            f',\n\tfriction, {self.friction}'
            f',\n\t\tpreload, {self.preload}'
            f',\n\t\t{self.friction_model}'
            f',\n\t\t{self.shape_function}'
            f'{revolute_hinge.element_footer()}'
        )
        self.assertEqual(str(revolute_hinge), expected_str)

    def test_revolute_hinge_output_option(self):
        # Test setting the output option to 'no'
        revolute_hinge = l.RevoluteHinge(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            output='no'
        )
        self.assertEqual(revolute_hinge.output, 'no')
        self.assertIn(',\n\toutput, no', str(revolute_hinge))

    # TODO: Check if MBVar class has any errors
    # def test_revolute_hinge_with_mbvar_node_labels(self):
    #     # Test creating a RevoluteHinge instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     revolute_hinge = l.RevoluteHinge(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         position_1=self.position_1,
    #         node_2_label=node_var_2,
    #         position_2=self.position_2
    #     )
    #     self.assertEqual(revolute_hinge.node_1_label, node_var_1)
    #     self.assertEqual(revolute_hinge.node_2_label, node_var_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_revolute_hinge_friction_validation(self):
        # Test that friction parameters are validated correctly
        # Case when friction is specified without friction_model or shape_function
        with self.assertRaises(Exception):
            l.RevoluteHinge(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                friction=self.friction
            )
        # Case when friction_model and shape_function are specified without friction
        with self.assertRaises(Exception):
            l.RevoluteHinge(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                friction_model=self.friction_model,
                shape_function=self.shape_function
            )
        # Valid case when friction is not specified and friction-related parameters are not specified
        try:
            revolute_hinge = l.RevoluteHinge(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2
            )
            self.assertIsNone(revolute_hinge.friction)
            self.assertIsNone(revolute_hinge.preload)
            self.assertIsNone(revolute_hinge.friction_model)
            self.assertIsNone(revolute_hinge.shape_function)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")

    def test_revolute_hinge_initial_theta(self):
        # Test setting initial_theta
        revolute_hinge = l.RevoluteHinge(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            initial_theta=self.initial_theta
        )
        self.assertEqual(revolute_hinge.initial_theta, self.initial_theta)
        self.assertIn(f',\n\tinitial theta, {self.initial_theta}', str(revolute_hinge))

class TestRevolutePin(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 20
        self.node_label = 1
        self.relative_offset = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.relative_orientation_mat = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='other node')
        self.absolute_pin_position = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='global')
        self.absolute_pin_orientation_mat = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='global')
        self.initial_theta = 0.0

    def test_revolute_pin_creation_valid(self):
        # Test creating a RevolutePin instance with all valid data
        revolute_pin = l.RevolutePin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            relative_orientation_mat=self.relative_orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_pin_orientation_mat=self.absolute_pin_orientation_mat,
            initial_theta=self.initial_theta
        )
        self.assertIsInstance(revolute_pin, l.RevolutePin)
        self.assertEqual(revolute_pin.node_label, self.node_label)
        self.assertEqual(revolute_pin.relative_offset, self.relative_offset)
        self.assertEqual(revolute_pin.relative_orientation_mat, self.relative_orientation_mat)
        self.assertEqual(revolute_pin.absolute_pin_position, self.absolute_pin_position)
        self.assertEqual(revolute_pin.absolute_pin_orientation_mat, self.absolute_pin_orientation_mat)
        self.assertEqual(revolute_pin.initial_theta, self.initial_theta)

    def test_revolute_pin_creation_without_optional_fields(self):
        # Test creating a RevolutePin instance without optional fields
        revolute_pin = l.RevolutePin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            absolute_pin_position=self.absolute_pin_position
        )
        self.assertIsInstance(revolute_pin, l.RevolutePin)
        self.assertIsNone(revolute_pin.relative_orientation_mat)
        self.assertIsNone(revolute_pin.absolute_pin_orientation_mat)
        self.assertIsNone(revolute_pin.initial_theta)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_revolute_pin_missing_required_fields(self):
        # Missing relative_offset
        with self.assertRaises(Exception):
            l.RevolutePin(
                idx=self.idx,
                node_label=self.node_label,
                absolute_pin_position=self.absolute_pin_position
            )
        # Missing absolute_pin_position
        with self.assertRaises(Exception):
            l.RevolutePin(
                idx=self.idx,
                node_label=self.node_label,
                relative_offset=self.relative_offset
            )
        # Missing node_label
        with self.assertRaises(Exception):
            l.RevolutePin(
                idx=self.idx,
                relative_offset=self.relative_offset,
                absolute_pin_position=self.absolute_pin_position
            )

    def test_revolute_pin_str_method(self):
        # Test the __str__ method of RevolutePin
        revolute_pin = l.RevolutePin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            relative_orientation_mat=self.relative_orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_pin_orientation_mat=self.absolute_pin_orientation_mat,
            initial_theta=self.initial_theta
        )
        expected_str = (
            f'{revolute_pin.element_header()}, revolute pin'
            f',\n\t{self.node_label}'
            f',\n\t\tposition, {self.relative_offset}'
            f',\n\t\torientation, {self.relative_orientation_mat}'
            f',\n\tposition, {self.absolute_pin_position}'
            f',\n\torientation, {self.absolute_pin_orientation_mat}'
            f',\n\tinitial theta, {self.initial_theta}'
            f'{revolute_pin.element_footer()}'
        )
        self.assertEqual(str(revolute_pin), expected_str)

    # # TODO: Check if MBVar class has errors
    # def test_revolute_pin_with_mbvar_node_label(self):
    #     # Test creating a RevolutePin instance with MBVar as node label
    #     node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #     revolute_pin = l.RevolutePin(
    #         idx=self.idx,
    #         node_label=node_var,
    #         relative_offset=self.relative_offset,
    #         absolute_pin_position=self.absolute_pin_position
    #     )
    #     self.assertEqual(revolute_pin.node_label, node_var)

    def test_revolute_pin_initial_theta(self):
        # Test setting initial_theta
        revolute_pin = l.RevolutePin(
            idx=self.idx,
            node_label=self.node_label,
            relative_offset=self.relative_offset,
            absolute_pin_position=self.absolute_pin_position,
            initial_theta=self.initial_theta
        )
        self.assertEqual(revolute_pin.initial_theta, self.initial_theta)
        self.assertIn(f',\n\tinitial theta, {self.initial_theta}', str(revolute_pin))

class TestRevoluteRotation(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 40
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='node')
        self.orientation_mat_1 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='global')
        self.orientation_mat_2 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='global')

    def test_revolute_rotation_creation_valid(self):
        # Test creating a RevoluteRotation instance with all valid data
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        self.assertIsInstance(revolute_rotation, l.RevoluteRotation)
        self.assertEqual(revolute_rotation.node_1_label, self.node_1_label)
        self.assertEqual(revolute_rotation.position_1, self.position_1)
        self.assertEqual(revolute_rotation.orientation_mat_1, self.orientation_mat_1)
        self.assertEqual(revolute_rotation.node_2_label, self.node_2_label)
        self.assertEqual(revolute_rotation.position_2, self.position_2)
        self.assertEqual(revolute_rotation.orientation_mat_2, self.orientation_mat_2)

    def test_revolute_rotation_creation_without_optional_fields(self):
        # Test creating a RevoluteRotation instance without optional fields
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )
        self.assertIsInstance(revolute_rotation, l.RevoluteRotation)
        self.assertIsNone(revolute_rotation.position_1)
        self.assertIsNone(revolute_rotation.orientation_mat_1)
        self.assertIsNone(revolute_rotation.position_2)
        self.assertIsNone(revolute_rotation.orientation_mat_2)

    @unittest.skipIf(pydantic is None, "Depends on Pydantic library")
    def test_revolute_rotation_missing_required_fields(self):
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.RevoluteRotation(
                idx=self.idx,
                node_2_label=self.node_2_label
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.RevoluteRotation(
                idx=self.idx,
                node_1_label=self.node_1_label
            )

    def test_revolute_rotation_str_method(self):
        # Test the __str__ method of RevoluteRotation
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        expected_str = (
            f'{revolute_rotation.element_header()}, revolute rotation'
            f',\n\t{self.node_1_label}'
            f',\n\t\tposition, {self.position_1}'
            f',\n\t\torientation, {self.orientation_mat_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\tposition, {self.position_2}'
            f',\n\t\torientation, {self.orientation_mat_2}'
            f'{revolute_rotation.element_footer()}'
        )
        self.assertEqual(str(revolute_rotation), expected_str)

    # # TODO: Check if MBVar class has errors
    # def test_revolute_rotation_with_mbvar_node_labels(self):
    #     # Test creating a RevoluteRotation instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     revolute_rotation = l.RevoluteRotation(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2
    #     )
    #     self.assertEqual(revolute_rotation.node_1_label, node_var_1)
    #     self.assertEqual(revolute_rotation.node_2_label, node_var_2)

    def test_revolute_rotation_optional_positions(self):
        # Test setting only position_1 and position_2
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2
        )
        self.assertEqual(revolute_rotation.position_1, self.position_1)
        self.assertIsNone(revolute_rotation.orientation_mat_1)
        self.assertEqual(revolute_rotation.position_2, self.position_2)
        self.assertIsNone(revolute_rotation.orientation_mat_2)
        # Check string representation
        expected_str = (
            f'{revolute_rotation.element_header()}, revolute rotation'
            f',\n\t{self.node_1_label}'
            f',\n\t\tposition, {self.position_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\tposition, {self.position_2}'
            f'{revolute_rotation.element_footer()}'
        )
        self.assertEqual(str(revolute_rotation), expected_str)

    def test_revolute_rotation_optional_orientations(self):
        # Test setting only orientation_mat_1 and orientation_mat_2
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            orientation_mat_2=self.orientation_mat_2
        )
        self.assertIsNone(revolute_rotation.position_1)
        self.assertEqual(revolute_rotation.orientation_mat_1, self.orientation_mat_1)
        self.assertIsNone(revolute_rotation.position_2)
        self.assertEqual(revolute_rotation.orientation_mat_2, self.orientation_mat_2)
        # Check string representation
        expected_str = (
            f'{revolute_rotation.element_header()}, revolute rotation'
            f',\n\t{self.node_1_label}'
            f',\n\t\torientation, {self.orientation_mat_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\torientation, {self.orientation_mat_2}'
            f'{revolute_rotation.element_footer()}'
        )
        self.assertEqual(str(revolute_rotation), expected_str)

    def test_revolute_rotation_output_option(self):
        # Test setting the output option to 'no'
        revolute_rotation = l.RevoluteRotation(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(revolute_rotation.output, 'no')
        self.assertIn(',\n\toutput, no', str(revolute_rotation))

class TestRod2(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 50
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='node')
        self.rod_length = 10.0
        self.const_law_valid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW,
            stiffness=1e6
        )
        self.const_law_invalid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=1e6
        )

    def test_rod2_creation_valid(self):
        # Test creating a Rod2 instance with all valid data
        rod2 = l.Rod2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            rod_length=self.rod_length,
            const_law=self.const_law_valid
        )
        self.assertIsInstance(rod2, l.Rod2)
        self.assertEqual(rod2.node_1_label, self.node_1_label)
        self.assertEqual(rod2.position_1, self.position_1)
        self.assertEqual(rod2.node_2_label, self.node_2_label)
        self.assertEqual(rod2.position_2, self.position_2)
        self.assertEqual(rod2.rod_length, self.rod_length)
        self.assertEqual(rod2.const_law, self.const_law_valid)

    def test_rod2_creation_without_optional_fields(self):
        # Test creating a Rod2 instance without optional position fields
        rod2 = l.Rod2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            rod_length='from nodes',
            const_law=self.const_law_valid
        )
        self.assertIsInstance(rod2, l.Rod2)
        self.assertIsNone(rod2.position_1)
        self.assertIsNone(rod2.position_2)
        self.assertEqual(rod2.rod_length, 'from nodes')

    @unittest.skipIf(pydantic is None, "Depends on Pydantic library")
    def test_rod2_missing_required_fields(self):
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.Rod2(
                idx=self.idx,
                node_2_label=self.node_2_label,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing rod_length
        with self.assertRaises(Exception):
            l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                const_law=self.const_law_valid
            )
        # Missing const_law
        with self.assertRaises(Exception):
            l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                rod_length=self.rod_length
            )

    def test_rod2_str_method(self):
        # Test the __str__ method of Rod2
        rod2 = l.Rod2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            rod_length=self.rod_length,
            const_law=self.const_law_valid
        )
        expected_str = (
            f'{rod2.element_header()}, rod'
            f',\n\t{self.node_1_label}'
            f',\n\t\tposition, {self.position_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\tposition, {self.position_2}'
            f',\n\t{self.rod_length}'
            f',\n\t{self.const_law_valid}'
            f'{rod2.element_footer()}'
        )
        self.maxDiff=None
        self.assertEqual(str(rod2), expected_str)

    def test_rod2_output_option(self):
        # Test setting the output option to 'no'
        rod2 = l.Rod2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            rod_length='from nodes',
            const_law=self.const_law_valid,
            output='no'
        )
        self.assertEqual(rod2.output, 'no')
        self.assertIn(',\n\toutput, no', str(rod2))

    # # TODO: Check if MBVar class has errors
    # def test_rod2_with_mbvar_node_labels(self):
    #     # Test creating a Rod2 instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     rod2 = l.Rod2(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2,
    #         rod_length=self.rod_length,
    #         const_law=self.const_law_valid
    #     )
    #     self.assertEqual(rod2.node_1_label, node_var_1)
    #     self.assertEqual(rod2.node_2_label, node_var_2)

    def test_rod2_rod_length_validation(self):
        # Test the rod_length field validator
        # Valid cases
        try:
            rod2_float_length = l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                rod_length=15.0,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod2_float_length.rod_length, 15.0)
            rod2_str_length = l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                rod_length='from nodes',
                const_law=self.const_law_valid
            )
            self.assertEqual(rod2_str_length.rod_length, 'from nodes')
            # TODO: Check if MBVar class has errors
            # rod2_mbvar_length = l.Rod2(
            #     idx=self.idx,
            #     node_1_label=self.node_1_label,
            #     node_2_label=self.node_2_label,
            #     rod_length=l.MBVar(name='rod_length_var', var_type='real', expression=20.0),
            #     const_law=self.const_law_valid
            # )
            # self.assertIsInstance(rod2_mbvar_length.rod_length, l.MBVar)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.Rod2(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    node_2_label=self.node_2_label,
                    rod_length='invalid string',
                    const_law=self.const_law_valid
                )

    def test_rod2_const_law_validation(self):
        # Test the const_law field validator
        # Valid case
        try:
            rod2 = l.Rod2(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod2.const_law, self.const_law_valid)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            # Invalid case: const_law is not a ConstitutiveLaw instance
            with self.assertRaises(Exception):
                l.Rod2(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    node_2_label=self.node_2_label,
                    rod_length=self.rod_length,
                    const_law='invalid_const_law'
                )
            # Invalid case: const_law with wrong law_type
            with self.assertRaises(Exception):
                l.Rod2(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    node_2_label=self.node_2_label,
                    rod_length=self.rod_length,
                    const_law=self.const_law_invalid
                )

class TestRodWithOffset(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 60
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='node')
        self.rod_length = 10.0
        self.const_law_valid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW,
            stiffness=1e6
        )
        self.const_law_invalid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=1e6
        )

    def test_rod_with_offset_creation_valid(self):
        # Test creating a RodWithOffset instance with all valid data
        rod_with_offset = l.RodWithOffset(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            rod_length=self.rod_length,
            const_law=self.const_law_valid
        )
        self.assertIsInstance(rod_with_offset, l.RodWithOffset)
        self.assertEqual(rod_with_offset.node_1_label, self.node_1_label)
        self.assertEqual(rod_with_offset.position_1, self.position_1)
        self.assertEqual(rod_with_offset.node_2_label, self.node_2_label)
        self.assertEqual(rod_with_offset.position_2, self.position_2)
        self.assertEqual(rod_with_offset.rod_length, self.rod_length)
        self.assertEqual(rod_with_offset.const_law, self.const_law_valid)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_rod_with_offset_missing_required_fields(self):
        # Missing position_1
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing position_2
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing rod_length
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                const_law=self.const_law_valid
            )
        # Missing const_law
        with self.assertRaises(Exception):
            l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length=self.rod_length
            )

    def test_rod_with_offset_str_method(self):
        # Test the __str__ method of RodWithOffset
        rod_with_offset = l.RodWithOffset(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            rod_length=self.rod_length,
            const_law=self.const_law_valid
        )
        expected_str = (
            f'{rod_with_offset.element_header()}, rod with offset'
            f',\n\t{self.node_1_label}'
            f',\n\t\t{self.position_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\t{self.position_2}'
            f',\n\t{self.rod_length}'
            f',\n\t{self.const_law_valid}'
            f'{rod_with_offset.element_footer()}'
        )
        self.assertEqual(str(rod_with_offset), expected_str)

    # # TODO: Check if MBVar class has errors
    # def test_rod_with_offset_with_mbvar_node_labels(self):
    #     # Test creating a RodWithOffset instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     rod_with_offset = l.RodWithOffset(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         position_1=self.position_1,
    #         node_2_label=node_var_2,
    #         position_2=self.position_2,
    #         rod_length=self.rod_length,
    #         const_law=self.const_law_valid
    #     )
    #     self.assertEqual(rod_with_offset.node_1_label, node_var_1)
    #     self.assertEqual(rod_with_offset.node_2_label, node_var_2)
    
    def test_rod_with_offset_rod_length_validation(self):
        # Test the rod_length field validator
        # Valid cases
        try:
            rod_float_length = l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length=15.0,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_float_length.rod_length, 15.0)
            rod_str_length = l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length='from nodes',
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_str_length.rod_length, 'from nodes')
            # TODO: Check if MBVar class has errors
            # rod_mbvar_length = l.RodWithOffset(
            #     idx=self.idx,
            #     node_1_label=self.node_1_label,
            #     position_1=self.position_1,
            #     node_2_label=self.node_2_label,
            #     position_2=self.position_2,
            #     rod_length=l.MBVar(name='rod_length_var', var_type='real', expression=20.0),
            #     const_law=self.const_law_valid
            # )
            # self.assertIsInstance(rod_mbvar_length.rod_length, l.MBVar)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodWithOffset(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    node_2_label=self.node_2_label,
                    position_2=self.position_2,
                    rod_length='invalid string',
                    const_law=self.const_law_valid
                )

    def test_rod_with_offset_const_law_validation(self):
        # Test the const_law field validator
        # Valid case
        try:
            rod_with_offset = l.RodWithOffset(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                node_2_label=self.node_2_label,
                position_2=self.position_2,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_with_offset.const_law, self.const_law_valid)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case: const_law is not a ConstitutiveLaw instance
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodWithOffset(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    node_2_label=self.node_2_label,
                    position_2=self.position_2,
                    rod_length=self.rod_length,
                    const_law='invalid_const_law'
                )
            # Invalid case: const_law with wrong law_type
            with self.assertRaises(Exception):
                l.RodWithOffset(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    node_2_label=self.node_2_label,
                    position_2=self.position_2,
                    rod_length=self.rod_length,
                    const_law=self.const_law_invalid
                )

    def test_rod_with_offset_output_option(self):
        # Test setting the output option to 'no'
        rod_with_offset = l.RodWithOffset(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            rod_length='from nodes',
            const_law=self.const_law_valid,
            output='no'
        )
        self.assertEqual(rod_with_offset.output, 'no')
        self.assertIn(',\n\toutput, no', str(rod_with_offset))

class TestRodBezier(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 70
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.position_2 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='node')
        self.position_3 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='node')
        self.position_4 = l.Position2(relative_position=[1.0, 1.0, 0.0], reference='node')
        self.rod_length = 10.0
        self.const_law_valid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.SCALAR_ISOTROPIC_LAW,
            stiffness=1e6
        )
        self.const_law_invalid = l.LinearElastic(
            law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW,
            stiffness=1e6
        )

    def test_rod_bezier_creation_valid(self):
        # Test creating a RodBezier instance with all valid data
        rod_bezier = l.RodBezier(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            position_2=self.position_2,
            node_2_label=self.node_2_label,
            position_3=self.position_3,
            position_4=self.position_4,
            rod_length=self.rod_length,
            const_law=self.const_law_valid,
            integration_order=5,
            integration_segments=4
        )
        self.assertIsInstance(rod_bezier, l.RodBezier)
        self.assertEqual(rod_bezier.node_1_label, self.node_1_label)
        self.assertEqual(rod_bezier.position_1, self.position_1)
        self.assertEqual(rod_bezier.position_2, self.position_2)
        self.assertEqual(rod_bezier.node_2_label, self.node_2_label)
        self.assertEqual(rod_bezier.position_3, self.position_3)
        self.assertEqual(rod_bezier.position_4, self.position_4)
        self.assertEqual(rod_bezier.rod_length, self.rod_length)
        self.assertEqual(rod_bezier.const_law, self.const_law_valid)
        self.assertEqual(rod_bezier.integration_order, 5)
        self.assertEqual(rod_bezier.integration_segments, 4)

    def test_rod_bezier_creation_with_defaults(self):
        # Test creating a RodBezier instance with default integration parameters
        rod_bezier = l.RodBezier(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            position_2=self.position_2,
            node_2_label=self.node_2_label,
            position_3=self.position_3,
            position_4=self.position_4,
            rod_length='from nodes',
            const_law=self.const_law_valid
        )
        self.assertIsInstance(rod_bezier, l.RodBezier)
        self.assertEqual(rod_bezier.integration_order, 2)
        self.assertEqual(rod_bezier.integration_segments, 3)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_rod_bezier_missing_required_fields(self):
        # Missing position_1
        with self.assertRaises(Exception):
            l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
        # Missing const_law
        with self.assertRaises(Exception):
            l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=self.rod_length
            )

    def test_rod_bezier_str_method(self):
        # Test the __str__ method of RodBezier
        rod_bezier = l.RodBezier(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            position_2=self.position_2,
            node_2_label=self.node_2_label,
            position_3=self.position_3,
            position_4=self.position_4,
            rod_length=self.rod_length,
            const_law=self.const_law_valid,
            integration_order=5,
            integration_segments=4
        )
        expected_str = (
            f'{rod_bezier.element_header()}, rod bezier'
            f',\n\t{self.node_1_label}'
            f',\n\t\t{self.position_1}'
            f',\n\t\t{self.position_2}'
            f',\n\t{self.node_2_label}'
            f',\n\t\t{self.position_3}'
            f',\n\t\t{self.position_4}'
            f',\n\t{self.rod_length}'
            f',\n\tintegration order, {rod_bezier.integration_order}'
            f',\n\tintegration segments, {rod_bezier.integration_segments}'
            f',\n\t{self.const_law_valid}'
            f'{rod_bezier.element_footer()}'
        )
        self.assertEqual(str(rod_bezier), expected_str)

    # TODO: Check if MBVar class has errors
    # def test_rod_bezier_with_mbvar_node_labels(self):
    #     # Test creating a RodBezier instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     rod_bezier = l.RodBezier(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         position_1=self.position_1,
    #         position_2=self.position_2,
    #         node_2_label=node_var_2,
    #         position_3=self.position_3,
    #         position_4=self.position_4,
    #         rod_length=self.rod_length,
    #         const_law=self.const_law_valid
    #     )
    #     self.assertEqual(rod_bezier.node_1_label, node_var_1)
    #     self.assertEqual(rod_bezier.node_2_label, node_var_2)

    def test_rod_bezier_rod_length_validation(self):
        # Test the rod_length field validator
        # Valid cases
        try:
            rod_float_length = l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=15.0,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_float_length.rod_length, 15.0)
            rod_str_length = l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length='from nodes',
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_str_length.rod_length, 'from nodes')
            # TODO: Check if MBVar class has errors
            # rod_mbvar_length = l.RodBezier(
            #     idx=self.idx,
            #     node_1_label=self.node_1_label,
            #     position_1=self.position_1,
            #     position_2=self.position_2,
            #     node_2_label=self.node_2_label,
            #     position_3=self.position_3,
            #     position_4=self.position_4,
            #     rod_length=l.MBVar(name='rod_length_var', var_type='real', expression=20.0),
            #     const_law=self.const_law_valid
            # )
            # self.assertIsInstance(rod_mbvar_length.rod_length, l.MBVar)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodBezier(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    position_2=self.position_2,
                    node_2_label=self.node_2_label,
                    position_3=self.position_3,
                    position_4=self.position_4,
                    rod_length='invalid string',
                    const_law=self.const_law_valid
                )

    def test_rod_bezier_const_law_validation(self):
        # Test the const_law field validator
        # Valid case
        try:
            rod_bezier = l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=self.rod_length,
                const_law=self.const_law_valid
            )
            self.assertEqual(rod_bezier.const_law, self.const_law_valid)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case: const_law is not a ConstitutiveLaw instance
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodBezier(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    position_2=self.position_2,
                    node_2_label=self.node_2_label,
                    position_3=self.position_3,
                    position_4=self.position_4,
                    rod_length=self.rod_length,
                    const_law='invalid_const_law'
                )
            # Invalid case: const_law with wrong law_type
            with self.assertRaises(Exception):
                l.RodBezier(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    position_2=self.position_2,
                    node_2_label=self.node_2_label,
                    position_3=self.position_3,
                    position_4=self.position_4,
                    rod_length=self.rod_length,
                    const_law=self.const_law_invalid
                )

    def test_rod_bezier_integration_order_validation(self):
        # Test the integration_order field validator
        # Valid case
        try:
            rod_bezier = l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=self.rod_length,
                const_law=self.const_law_valid,
                integration_order=5
            )
            self.assertEqual(rod_bezier.integration_order, 5)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case: integration_order out of bounds
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodBezier(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    position_2=self.position_2,
                    node_2_label=self.node_2_label,
                    position_3=self.position_3,
                    position_4=self.position_4,
                    rod_length=self.rod_length,
                    const_law=self.const_law_valid,
                    integration_order=11  # Invalid value
                )

    def test_rod_bezier_integration_segments_validation(self):
        # Test the integration_segments field validator
        # Valid case
        try:
            rod_bezier = l.RodBezier(
                idx=self.idx,
                node_1_label=self.node_1_label,
                position_1=self.position_1,
                position_2=self.position_2,
                node_2_label=self.node_2_label,
                position_3=self.position_3,
                position_4=self.position_4,
                rod_length=self.rod_length,
                const_law=self.const_law_valid,
                integration_segments=5
            )
            self.assertEqual(rod_bezier.integration_segments, 5)
        except Exception as e:
            self.fail(f"Unexpected exception occurred: {e}")
        # Invalid case: integration_segments not positive
        if pydantic is None:
            self.skipTest("Pydantic not available, skipping invalid input test")
        else:
            with self.assertRaises(Exception):
                l.RodBezier(
                    idx=self.idx,
                    node_1_label=self.node_1_label,
                    position_1=self.position_1,
                    position_2=self.position_2,
                    node_2_label=self.node_2_label,
                    position_3=self.position_3,
                    position_4=self.position_4,
                    rod_length=self.rod_length,
                    const_law=self.const_law_valid,
                    integration_segments=0  # Invalid value
                )

    def test_rod_bezier_output_option(self):
        # Test setting the output option to 'no'
        rod_bezier = l.RodBezier(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            position_2=self.position_2,
            node_2_label=self.node_2_label,
            position_3=self.position_3,
            position_4=self.position_4,
            rod_length='from nodes',
            const_law=self.const_law_valid,
            output='no'
        )
        self.assertEqual(rod_bezier.output, 'no')
        self.assertIn(',\n\toutput, no', str(rod_bezier))
        
class TestSphericalHinge2(unittest.TestCase):
    def setUp(self):
        # Common variables used in tests
        self.idx = 100
        self.node_1_label = 1
        self.node_2_label = 2
        self.position_1 = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='node')
        self.orientation_mat_1 = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='global')
        self.position_2 = l.Position2(relative_position=[1.0, 1.0, 1.0], reference='node')
        self.orientation_mat_2 = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='global')

    def test_spherical_hinge_creation_valid(self):
        # Test creating a SphericalHinge instance with all valid data
        spherical_hinge = l.SphericalHinge2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        self.assertIsInstance(spherical_hinge, l.SphericalHinge2)
        self.assertEqual(spherical_hinge.node_1_label, self.node_1_label)
        self.assertEqual(spherical_hinge.position_1, self.position_1)
        self.assertEqual(spherical_hinge.orientation_mat_1, self.orientation_mat_1)
        self.assertEqual(spherical_hinge.node_2_label, self.node_2_label)
        self.assertEqual(spherical_hinge.position_2, self.position_2)
        self.assertEqual(spherical_hinge.orientation_mat_2, self.orientation_mat_2)

    def test_spherical_hinge_creation_without_optional_fields(self):
        # Test creating a SphericalHinge instance without optional fields
        spherical_hinge = l.SphericalHinge2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label
        )
        self.assertIsInstance(spherical_hinge, l.SphericalHinge2)
        self.assertIsNone(spherical_hinge.position_1)
        self.assertIsNone(spherical_hinge.orientation_mat_1)
        self.assertIsNone(spherical_hinge.position_2)
        self.assertIsNone(spherical_hinge.orientation_mat_2)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_spherical_hinge_missing_required_fields(self):
        # Missing node_1_label
        with self.assertRaises(Exception):
            l.SphericalHinge2(
                idx=self.idx,
                node_2_label=self.node_2_label
            )
        # Missing node_2_label
        with self.assertRaises(Exception):
            l.SphericalHinge2(
                idx=self.idx,
                node_1_label=self.node_1_label
            )

    def test_spherical_hinge_str_method(self):
        # Test the __str__ method of SphericalHinge
        spherical_hinge = l.SphericalHinge2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            position_1=self.position_1,
            orientation_mat_1=self.orientation_mat_1,
            node_2_label=self.node_2_label,
            position_2=self.position_2,
            orientation_mat_2=self.orientation_mat_2
        )
        expected_str = (
            f'{spherical_hinge.element_header()}, spherical hinge'
            f',\n\t{self.node_1_label}'
            f',\n\t\tposition, {self.position_1}'
            f',\n\t\torientation, {self.orientation_mat_1}'
            f',\n\t{self.node_2_label}'
            f',\n\t\tposition, {self.position_2}'
            f',\n\t\torientation, {self.orientation_mat_2}'
            f'{spherical_hinge.element_footer()}'
        )
        self.assertEqual(str(spherical_hinge), expected_str)

    # # TODO: Check if MBVar class has errors
    # def test_spherical_hinge_with_mbvar_node_labels(self):
    #     # Test creating a SphericalHinge instance with MBVar as node labels
    #     node_var_1 = l.MBVar(name='node_var_1', var_type='integer', expression=100)
    #     node_var_2 = l.MBVar(name='node_var_2', var_type='integer', expression=200)
    #     spherical_hinge = l.SphericalHinge2(
    #         idx=self.idx,
    #         node_1_label=node_var_1,
    #         node_2_label=node_var_2
    #     )
    #     self.assertEqual(spherical_hinge.node_1_label, node_var_1)
    #     self.assertEqual(spherical_hinge.node_2_label, node_var_2)

    def test_spherical_hinge_output_option(self):
        # Test setting the output option to 'no'
        spherical_hinge = l.SphericalHinge2(
            idx=self.idx,
            node_1_label=self.node_1_label,
            node_2_label=self.node_2_label,
            output='no'
        )
        self.assertEqual(spherical_hinge.output, 'no')
        self.assertIn(',\n\toutput, no', str(spherical_hinge))

class TestSphericalPin(unittest.TestCase):

    def setUp(self):
        # Common variables used in tests
        self.idx = 100
        self.node_label = 1
        self.position = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.orientation_mat = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='node')
        self.absolute_pin_position = l.Position2(relative_position=[2.0, 2.0, 2.0], reference='global')
        self.absolute_orientation_mat = l.Position2(relative_position=[0.0, 1.0, 0.0], reference='node')

    def test_spherical_pin_creation_valid(self):
        # Test creating a SphericalPin instance with all valid data
        spherical_pin = l.SphericalPin(
            idx=self.idx,
            node_label=self.node_label,
            position=self.position,
            orientation_mat=self.orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_orientation_mat=self.absolute_orientation_mat
        )
        self.assertIsInstance(spherical_pin, l.SphericalPin)
        self.assertEqual(spherical_pin.node_label, self.node_label)
        self.assertEqual(spherical_pin.position, self.position)
        self.assertEqual(spherical_pin.orientation_mat, self.orientation_mat)
        self.assertEqual(spherical_pin.absolute_pin_position, self.absolute_pin_position)
        self.assertEqual(spherical_pin.absolute_orientation_mat, self.absolute_orientation_mat)

    def test_spherical_pin_creation_without_optional_fields(self):
        # Test creating a SphericalPin instance without optional fields
        spherical_pin = l.SphericalPin(
            idx=self.idx,
            node_label=self.node_label,
            absolute_pin_position=self.absolute_pin_position
        )
        self.assertIsInstance(spherical_pin, l.SphericalPin)
        self.assertEqual(spherical_pin.node_label, self.node_label)
        self.assertIsNone(spherical_pin.position)
        self.assertIsNone(spherical_pin.orientation_mat)
        self.assertEqual(spherical_pin.absolute_pin_position, self.absolute_pin_position)
        self.assertIsNone(spherical_pin.absolute_orientation_mat)

    def test_spherical_pin_str_method(self):
        # Test the __str__ method of SphericalPin
        spherical_pin = l.SphericalPin(
            idx=self.idx,
            node_label=self.node_label,
            position=self.position,
            orientation_mat=self.orientation_mat,
            absolute_pin_position=self.absolute_pin_position,
            absolute_orientation_mat=self.absolute_orientation_mat
        )
        expected_str = (
            f'{spherical_pin.element_header()}, spherical pin'
            f',\n\t{self.node_label}'
            f',\n\t\tposition, {self.position}'
            f',\n\t\torientation, {self.orientation_mat}'
            f',\n\tposition, {self.absolute_pin_position}'
            f',\n\torientation, {self.absolute_orientation_mat}'
            f'{spherical_pin.element_footer()}'
        )
        self.assertEqual(str(spherical_pin), expected_str)

    def test_spherical_pin_with_minimal_data(self):
        # Test creating SphericalPin with only required fields
        spherical_pin = l.SphericalPin(
            idx=self.idx,
            node_label=self.node_label,
            absolute_pin_position=self.absolute_pin_position
        )
        self.assertEqual(spherical_pin.node_label, self.node_label)
        self.assertIsNone(spherical_pin.position)
        self.assertIsNone(spherical_pin.orientation_mat)
        self.assertEqual(spherical_pin.absolute_pin_position, self.absolute_pin_position)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_spherical_pin_missing_required_fields(self):
        # Test creating a SphericalPin instance without required fields
        with self.assertRaises(Exception):
            l.SphericalPin(
                idx=self.idx,
                node_label=self.node_label
            )

    def test_spherical_pin_output_option(self):
        # Test setting the output option to 'no'
        spherical_pin = l.SphericalPin(
            idx=self.idx,
            node_label=self.node_label,
            absolute_pin_position=self.absolute_pin_position,
            output='no'
        )
        self.assertEqual(spherical_pin.output, 'no')
        self.assertIn(',\n\toutput, no', str(spherical_pin))

    # TODO: Check if MBVar class has errors
    # # Test with MBVar (assuming MBVar implementation is valid)
    # def test_spherical_pin_with_mbvar_node_label(self):
    #     node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #     spherical_pin = l.SphericalPin(
    #         idx=self.idx,
    #         node_label=node_var,
    #         absolute_pin_position=self.absolute_pin_position
    #     )
    #     self.assertEqual(spherical_pin.node_label, node_var)

class TestViscousBody(unittest.TestCase):

    def setUp(self):
        # Common variables used in tests
        self.idx = 100
        self.node_label = 1
        self.position = l.Position2(relative_position=[0.0, 0.0, 0.0], reference='global')
        self.orientation_mat = l.Position2(relative_position=[1.0, 0.0, 0.0], reference='node')
        self.const_law_valid = l.LinearViscous(viscosity=5.0, law_type=l.ConstitutiveLaw.LawType.D6_ISOTROPIC_LAW)
        self.const_law_invalid = l.LinearViscous(viscosity=5.0, law_type=l.ConstitutiveLaw.LawType.D3_ISOTROPIC_LAW)
        self.named_const_law = l.NamedConstitutiveLaw("linear viscous generic")

    def test_viscous_body_creation_valid(self):
        # Test creating a ViscousBody instance with valid 6D const_law
        viscous_body = l.ViscousBody(
            idx=self.idx,
            node_label=self.node_label,
            position=self.position,
            orientation_mat=self.orientation_mat,
            const_law=self.const_law_valid
        )
        self.assertIsInstance(viscous_body, l.ViscousBody)
        self.assertEqual(viscous_body.node_label, self.node_label)
        self.assertEqual(viscous_body.const_law, self.const_law_valid)

    def test_viscous_body_with_named_const_law(self):
    # Test creating a ViscousBody instance with NamedConstitutiveLaw
        viscous_body = l.ViscousBody(
            idx=self.idx,
            node_label=self.node_label,
            const_law=self.named_const_law
        )
        self.assertIsInstance(viscous_body, l.ViscousBody)
        self.assertEqual(viscous_body.const_law, self.named_const_law)


    def test_viscous_body_creation_without_optional_fields(self):
        # Test creating a ViscousBody instance without optional fields
        viscous_body = l.ViscousBody(
            idx=self.idx,
            node_label=self.node_label,
            const_law=self.const_law_valid
        )
        self.assertIsInstance(viscous_body, l.ViscousBody)
        self.assertEqual(viscous_body.node_label, self.node_label)
        self.assertIsNone(viscous_body.position)
        self.assertIsNone(viscous_body.orientation_mat)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_viscous_body_invalid_const_law(self):
        # Test creating a ViscousBody instance with invalid 3D const_law
        with self.assertRaises(Exception) as context:
            l.ViscousBody(
                idx=self.idx,
                node_label=self.node_label,
                const_law=self.const_law_invalid
            )
        self.assertIn("const_law must be a 6D constitutive law with law_type 'D6_ISOTROPIC_LAW'", str(context.exception))

    def test_viscous_body_str_method(self):
        # Test the __str__ method of ViscousBody
        viscous_body = l.ViscousBody(
            idx=self.idx,
            node_label=self.node_label,
            position=self.position,
            orientation_mat=self.orientation_mat,
            const_law=self.const_law_valid
        )
        expected_str = (
            f'{viscous_body.element_header()}, viscous body'
            f',\n\t{self.node_label}'
            f',\n\t\tposition, {self.position}'
            f',\n\t\torientation, {self.orientation_mat}'
            f',\n\t{self.const_law_valid}'
            f'{viscous_body.element_footer()}'
        )
        self.assertEqual(str(viscous_body), expected_str)

    def test_viscous_body_output_option(self):
        # Test setting the output option to 'no'
        viscous_body = l.ViscousBody(
            idx=self.idx,
            node_label=self.node_label,
            const_law=self.const_law_valid,
            output='no'
        )
        self.assertEqual(viscous_body.output, 'no')
        self.assertIn(',\n\toutput, no', str(viscous_body))

    # # TODO: Check if MBVar class has errors
    # def test_viscous_body_with_mbvar_node_label(self):
    #     node_var = l.MBVar(name='node_var', var_type='integer', expression=100)
    #     viscous_body = l.ViscousBody(
    #         idx=self.idx,
    #         node_label=node_var,
    #         const_law=self.const_law_valid
    #     )
    #     self.assertEqual(viscous_body.node_label, node_var)

    @unittest.skipIf(pydantic is None, "depends on library, since it doesn't prevent correct models from running")
    def test_viscous_body_missing_required_field(self):
    # Test missing required 'const_law' field
        with self.assertRaises(Exception):
            viscous_body = l.ViscousBody(
                idx=self.idx,
                node_label=self.node_label
                # const_law is missing here
            )

if __name__ == '__main__':
    unittest.main()
