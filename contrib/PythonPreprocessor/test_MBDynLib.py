from io import StringIO
import unittest

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
