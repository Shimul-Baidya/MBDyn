import unittest
import sys
import os
import tempfile

# Adjust path to import from the RPM package
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import MBDynLib as l
from RPM.core.processor import ModelProcessor
from RPM.components.base import Rotorcraft
from RPM.components.airframe import Airframe
from RPM.components.rotor import Rotor
from RPM.components.blade import Blade
from RPM.components.hub import Hub
from RPM.core.datatypes import ReferenceSystem, PhysicalQuantity


class TestModelProcessor(unittest.TestCase):
    """Comprehensive tests for ModelProcessor label generation and processing."""

    def test_01_initialization_loads_labels(self):
        """Test that ModelProcessor correctly loads Labels.yml."""
        processor = ModelProcessor()
        self.assertIsNotNone(processor.label_conventions)
        self.assertIsInstance(processor.label_conventions, dict)
        # Verify it loaded actual file with expected keys
        self.assertIn('airframe_1', processor.label_conventions)
        self.assertIn('rotor_1', processor.label_conventions)
        self.assertIn('blade', processor.label_conventions)

    def test_02_root_component_fixed_label(self):
        """Test processing a root component with fixed label (AIRFRAME_1)."""
        rotorcraft = Rotorcraft(name='TestCraft')
        airframe = Airframe(
            reference_system=ReferenceSystem(name='airframe_1', base_reference='global')
        )
        rotorcraft.add_root_component(airframe)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Check label assignment
        self.assertIsNotNone(airframe.reference_system.mbdyn_label)
        self.assertEqual(airframe.reference_system.mbdyn_label.name, 'AIRFRAME_1')
        
        # Check value (from actual Labels.yml: AIRFRAME_1: 80)
        self.assertEqual(processor._value_registry['AIRFRAME_1'], 80)
        
        # Check MBVar registry
        self.assertIn('AIRFRAME_1', processor._mbvar_registry)
        self.assertIsInstance(processor._mbvar_registry['AIRFRAME_1'], l.MBVar)
    def test_03_root_component_naming_with_digit_suffix(self):
        """Test that 'rotor_1' creates 'ROTOR_1' not 'ROTOR_1_1'."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Should be 'ROTOR_1' not 'ROTOR_1_1'
        self.assertEqual(rotor.reference_system.mbdyn_label.name, 'ROTOR_1')
        # From actual Labels.yml: ROTOR_1: 10000
        self.assertEqual(processor._value_registry['ROTOR_1'], 10000)

    def test_04_duplicate_components_counter_system(self):
        """Test that creating multiple components with same name appends counter."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor1 = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        rotor2 = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        rotorcraft.add_root_component(rotor1)
        rotorcraft.add_root_component(rotor2)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # First rotor gets clean name
        self.assertEqual(rotor1.reference_system.mbdyn_label.name, 'ROTOR_1')
        # Second rotor gets counter suffix
        self.assertEqual(rotor2.reference_system.mbdyn_label.name, 'ROTOR_1_2')
        
        # Both should have same base value from 'rotor' conventions
        self.assertEqual(processor._value_registry['ROTOR_1'], 10000)
        self.assertEqual(processor._value_registry['ROTOR_1_2'], 10000)

    def test_05_sub_component_fixed_label(self):
        """Test sub-component with fixed label (hub under rotor)."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        hub = Hub(
            reference_system=ReferenceSystem(name='hub', base_reference='ROTOR_1')
        )
        rotor.add_sub_component(hub)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Hub should have parent prefix
        self.assertEqual(hub.reference_system.mbdyn_label.name, 'ROTOR_1_HUB_1')
        
        # Hub value should be parent + base: 10000 + 200 = 10200
        self.assertEqual(processor._value_registry['ROTOR_1_HUB_1'], 10200)

    def test_06_sub_component_incremental_label(self):
        """Test sub-component with incremental label (blade under rotor)."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        
        # Add 3 blades
        for i in range(3):
            blade = Blade(
                reference_system=ReferenceSystem(name='blade', base_reference='ROTOR_1')
            )
            rotor.add_sub_component(blade)
        
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Blade values should be incremental (CURR_BLADE: 1000)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1'], 11000)  # 10000 + 1000*1
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_2'], 12000)  # 10000 + 1000*2
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_3'], 13000)  # 10000 + 1000*3

    def test_07_internal_label_node(self):
        """Test creating internal node labels."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='ROTOR_1')
        )
        rotor.add_sub_component(blade)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create internal node labels
        node0 = processor.create_internal_label(blade, 'NODE', 0)
        node1 = processor.create_internal_label(blade, 'NODE', 1)
        
        # Check names
        self.assertEqual(node0.name, 'ROTOR_1_BLADE_1_NODE_0')
        self.assertEqual(node1.name, 'ROTOR_1_BLADE_1_NODE_1')
        
        # Check values (no offset for nodes)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_NODE_0'], 11000)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_NODE_1'], 11001)

    def test_08_internal_label_ref_feathering(self):
        """Test creating feathering reference labels with offset."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='ROTOR_1')
        )
        rotor.add_sub_component(blade)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create feathering ref labels
        feath0 = processor.create_internal_label(blade, 'REF_FEATHERING', 0)
        feath1 = processor.create_internal_label(blade, 'REF_FEATHERING', 1)
        
        # Check names
        self.assertEqual(feath0.name, 'ROTOR_1_BLADE_1_REF_FEATHERING_0')
        
        # Check values (40000 offset)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_REF_FEATHERING_0'], 51000)  # 11000 + 40000
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_REF_FEATHERING_1'], 51001)  # 11000 + 40000 + 1

    def test_09_internal_label_ref_neutral(self):
        """Test creating neutral reference labels with offset."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='ROTOR_1')
        )
        rotor.add_sub_component(blade)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create neutral ref labels
        neutral0 = processor.create_internal_label(blade, 'REF_NEUTRAL', 0)
        
        # Check values (80000 offset)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_REF_NEUTRAL_0'], 91000)  # 11000 + 80000

    def test_10_internal_label_ref_precone(self):
        """Test creating precone reference label."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='ROTOR_1')
        )
        rotor.add_sub_component(blade)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create precone ref label
        precone = processor.create_internal_label(blade, 'REF_PRECONE', 0)
        
        # Check values (offset = 1)
        self.assertEqual(processor._value_registry['ROTOR_1_BLADE_1_REF_PRECONE_0'], 11001)  # 11000 + 1

    def test_11_internal_label_deduplication(self):
        """Test that calling create_internal_label twice returns same MBVar."""
        rotorcraft = Rotorcraft(name='TestCraft')
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create same internal label twice
        node1 = processor.create_internal_label(rotor, 'NODE', 0)
        node2 = processor.create_internal_label(rotor, 'NODE', 0)
        
        # Should return same object
        self.assertIs(node1, node2)

    def test_12_missing_label_warning(self):
        """Test that missing label defaults to 0 with warning."""
        rotorcraft = Rotorcraft(name='TestCraft')
        unknown = Airframe(
            reference_system=ReferenceSystem(name='unknown_component', base_reference='global')
        )
        rotorcraft.add_root_component(unknown)
        
        processor = ModelProcessor()
        
        # Capture print output
        import io
        from contextlib import redirect_stdout
        
        f = io.StringIO()
        with redirect_stdout(f):
            processor.process_model(rotorcraft)
        output = f.getvalue()
        
        # Check warning was printed
        self.assertIn("Warning", output)
        self.assertIn("unknown_component", output)
        
        # Check default value
        self.assertEqual(processor._value_registry['UNKNOWN_COMPONENT_1'], 0)

    def test_13_write_labels_set(self):
        """Test writing labels.set file."""
        rotorcraft = Rotorcraft(name='TestCraft')
        airframe = Airframe(
            reference_system=ReferenceSystem(name='airframe_1', base_reference='global')
        )
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='global'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        rotorcraft.add_root_component(airframe)
        rotorcraft.add_root_component(rotor)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Write to temporary file
        temp_output = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.set')
        temp_output.close()
        
        try:
            processor.write_labels_set(temp_output.name)
            
            # Read and verify
            with open(temp_output.name, 'r') as f:
                content = f.read()
            
            # Check format
            self.assertIn("#beginpreprocess", content)
            self.assertIn("#endpreprocess", content)
            self.assertIn("ConstMBVar('AIRFRAME_1', 'integer', 80)", content)
            self.assertIn("ConstMBVar('ROTOR_1', 'integer', 10000)", content)
        
        finally:
            os.remove(temp_output.name)

    def test_14_body_and_node_same_value(self):
        """Test that bodies and nodes can share same label value."""
        rotorcraft = Rotorcraft(name='TestCraft')
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='rotor')
        )
        rotorcraft.add_root_component(blade)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Create node and body with same counter
        node0 = processor.create_internal_label(blade, 'NODE', 0)
        body0 = processor.create_internal_label(blade, 'BODY', 0)
        
        # Both should have same value (offset = 0 for both)
        self.assertEqual(
            processor._value_registry['BLADE_1_NODE_0'],
            processor._value_registry['BLADE_1_BODY_0']
        )

    def test_15_complex_hierarchy(self):
        """Test complex component hierarchy with multiple levels."""
        rotorcraft = Rotorcraft(name='TestCraft')
        
        # Root: Airframe
        airframe = Airframe(
            reference_system=ReferenceSystem(name='airframe_1', base_reference='global')
        )
        
        # Level 1: Rotor under Airframe
        rotor = Rotor(
            reference_system=ReferenceSystem(name='rotor_1', base_reference='AIRFRAME_1'),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.25),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05)
        )
        airframe.add_sub_component(rotor)
        
        # Level 2: Hub under Rotor
        hub = Hub(
            reference_system=ReferenceSystem(name='hub', base_reference='AIRFRAME_1_ROTOR_1')
        )
        rotor.add_sub_component(hub)
        
        # Level 2: Blade under Rotor
        blade = Blade(
            reference_system=ReferenceSystem(name='blade', base_reference='AIRFRAME_1_ROTOR_1')
        )
        rotor.add_sub_component(blade)
        
        rotorcraft.add_root_component(airframe)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Check naming hierarchy
        self.assertEqual(airframe.reference_system.mbdyn_label.name, 'AIRFRAME_1')
        self.assertEqual(rotor.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1')  # Strips digit suffix
        self.assertEqual(hub.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1_HUB_1')  # Adds counter
        self.assertEqual(blade.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1_BLADE_1')  # Adds counter
        
        # Check value propagation
        self.assertEqual(processor._value_registry['AIRFRAME_1'], 80)
        # Rotor: 80 + 10000 = 10080
        self.assertEqual(processor._value_registry['AIRFRAME_1_ROTOR_1'], 10080)
        # Hub: 10080 + 200 = 10280
        self.assertEqual(processor._value_registry['AIRFRAME_1_ROTOR_1_HUB_1'], 10280)
        # Blade: 10080 + (1000 * 1) = 11080
        self.assertEqual(processor._value_registry['AIRFRAME_1_ROTOR_1_BLADE_1'], 11080)


if __name__ == '__main__':
    unittest.main(verbosity=2)