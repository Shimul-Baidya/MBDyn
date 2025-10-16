import unittest
import sys
import os

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from RPM.components.rotor import Rotor
from RPM.components.base import Rotorcraft
from RPM.core.datatypes import ReferenceSystem, PhysicalQuantity
from RPM.core.processor import ModelProcessor


class TestRotor(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures including processor and model."""
        base_ref = 'AIRFRAME_1'
        self.reference_system = ReferenceSystem(
            name='rotor_1',
            base_reference=base_ref,
            position_wrt_base=l.Position(relative_position=[0., 0., 0.116], reference=base_ref),
            orientation_wrt_base=l.Position(relative_position=l.eye(), reference=base_ref),
            velocity_wrt_base=l.Position(relative_position=l.null(), reference=base_ref),
            angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference=base_ref),
        )
        
        # Create rotor with required physical parameters
        self.rotor = Rotor(
            reference_system=self.reference_system,
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.250),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05569)
        )
        
        # Create processor and process the model
        self.processor = ModelProcessor()
        self.rotorcraft = Rotorcraft(name="TestRotorcraft")
        self.rotorcraft.add_root_component(self.rotor)
        self.processor.process_model(self.rotorcraft)

    def test_initialization(self):
        """Test that the Rotor component is initialized correctly."""
        self.assertIsInstance(self.rotor, Rotor)
        self.assertEqual(self.rotor.reference_system.name, 'rotor_1')
        self.assertEqual(self.rotor.reference_system.base_reference, 'AIRFRAME_1')
        self.assertEqual(self.rotor.n_blades, 3)
        self.assertEqual(self.rotor.radius.value, 1.250)
        self.assertEqual(self.rotor.radius.unit, 'm')
        self.assertEqual(self.rotor.precone.value, 2.75)
        self.assertEqual(self.rotor.precone.unit, 'deg')
        self.assertEqual(self.rotor.precone_start.value, 0.05569)
        self.assertEqual(self.rotor.precone_start.unit, 'adim')

    def test_label_assignment(self):
        """Test that the processor correctly assigns labels."""
        # After processing, the rotor should have an mbdyn_label
        self.assertIsNotNone(self.rotor.reference_system.mbdyn_label)
        self.assertEqual(self.rotor.reference_system.mbdyn_label.name, 'ROTOR_1')

    def test_create_references(self):
        """Test the _create_references method."""
        references = self.rotor._create_references(self.processor)
        self.assertEqual(len(references), 1, "Rotor should create exactly 1 reference")
        
        ref = references[0]
        self.assertIsInstance(ref, l.Reference)
        
        # The reference should use the label assigned by the processor
        self.assertEqual(ref.idx, self.rotor.reference_system.mbdyn_label)
        self.assertEqual(ref.position, self.reference_system.position_wrt_base)
        self.assertEqual(ref.orientation, self.reference_system.orientation_wrt_base)
        self.assertEqual(ref.velocity, self.reference_system.velocity_wrt_base)
        self.assertEqual(ref.angular_velocity, self.reference_system.angular_velocity_wrt_base)

    def test_create_nodes(self):
        """Test the _create_nodes method returns empty list."""
        nodes = self.rotor._create_nodes(self.processor)
        self.assertEqual(len(nodes), 0, "Rotor is a non-physical container and should not create nodes")
        self.assertIsInstance(nodes, list)

    def test_create_elements(self):
        """Test the _create_elements method returns empty list."""
        elements = self.rotor._create_elements(self.processor)
        self.assertEqual(len(elements), 0, "Rotor is a non-physical container and should not create elements")
        self.assertIsInstance(elements, list)

    def test_generate_mbdyn_entities(self):
        """Test the full entity generation process."""
        self.rotor.generate_mbdyn_entities(self.processor)
        
        # Verify the internal storage
        self.assertEqual(len(self.rotor._references), 1)
        self.assertEqual(len(self.rotor._nodes), 0)
        self.assertEqual(len(self.rotor._elements), 0)

    def test_collect_entities(self):
        """Test that collect_entities returns correct structure."""
        self.rotor.generate_mbdyn_entities(self.processor)
        refs, nodes, elements = self.rotor.collect_entities()
        
        self.assertEqual(len(refs), 1, "Should collect 1 reference")
        self.assertEqual(len(nodes), 0, "Should collect 0 nodes")
        self.assertEqual(len(elements), 0, "Should collect 0 elements")

    def test_rotor_with_sub_components(self):
        """Test that rotor can have sub-components (Hub, Mast, Blades, etc.)."""
        # Import Hub for testing
        from RPM.components.hub import Hub
        
        # Create a hub as a sub-component
        hub = Hub(
            reference_system=ReferenceSystem(
                name='hub',
                base_reference='ROTOR_1'
            )
        )
        
        self.rotor.add_sub_component(hub)
        self.assertEqual(len(self.rotor.sub_components), 1)
        self.assertEqual(self.rotor.sub_components[0], hub)

    def test_multiple_rotors(self):
        """Test creating multiple rotors (e.g., tiltrotor configuration)."""
        rotor2 = Rotor(
            reference_system=ReferenceSystem(
                name='rotor_2',
                base_reference='AIRFRAME_2',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.116], reference='AIRFRAME_2')
            ),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.250),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05569)
        )
        
        rotorcraft = Rotorcraft(name="TiltrotorTest")
        rotorcraft.add_root_component(self.rotor)
        rotorcraft.add_root_component(rotor2)
        
        processor = ModelProcessor()
        processor.process_model(rotorcraft)
        
        # Both rotors should have different labels
        self.assertEqual(self.rotor.reference_system.mbdyn_label.name, 'ROTOR_1')
        self.assertEqual(rotor2.reference_system.mbdyn_label.name, 'ROTOR_2')

    def test_rotor_physical_parameters(self):
        """Test that physical parameters are correctly stored and accessible."""
        # Test radius
        self.assertIsInstance(self.rotor.radius, PhysicalQuantity)
        self.assertEqual(self.rotor.radius.unit, 'm')
        self.assertIsInstance(self.rotor.radius.value, float)
        
        # Test precone
        self.assertIsInstance(self.rotor.precone, PhysicalQuantity)
        self.assertEqual(self.rotor.precone.unit, 'deg')
        self.assertIsInstance(self.rotor.precone.value, float)
        
        # Test precone_start
        self.assertIsInstance(self.rotor.precone_start, PhysicalQuantity)
        self.assertEqual(self.rotor.precone_start.unit, 'adim')
        self.assertIsInstance(self.rotor.precone_start.value, float)
        
        # Test n_blades
        self.assertIsInstance(self.rotor.n_blades, int)
        self.assertGreater(self.rotor.n_blades, 0)

    def test_rotor_as_parameter_holder(self):
        """Test that rotor serves as a parameter container for child components."""
        # The rotor should hold parameters that child components can access
        # This tests the design pattern where Rotor is a non-physical container
        
        # Child components should be able to access parent's parameters
        self.assertTrue(hasattr(self.rotor, 'n_blades'))
        self.assertTrue(hasattr(self.rotor, 'radius'))
        self.assertTrue(hasattr(self.rotor, 'precone'))
        self.assertTrue(hasattr(self.rotor, 'precone_start'))
        

if __name__ == '__main__':
    unittest.main()
