import unittest
import sys
import os

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from RPM.components.airframe import Airframe
from RPM.components.base import Rotorcraft
from RPM.core.datatypes import ReferenceSystem
from RPM.core.processor import ModelProcessor


class TestAirframe(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures including processor and model."""
        base_ref = 'global'
        self.reference_system = ReferenceSystem(
            name='airframe_1',
            base_reference=base_ref,
            position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference=base_ref),
            orientation_wrt_base=l.Position(relative_position=l.eye(), reference=base_ref),
            velocity_wrt_base=l.Position(relative_position=l.null(), reference=base_ref),
            angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference=base_ref),
        )
        self.airframe = Airframe(reference_system=self.reference_system)
        
        # Create processor and process the model
        self.processor = ModelProcessor()
        self.rotorcraft = Rotorcraft(name="TestRotorcraft")
        self.rotorcraft.add_root_component(self.airframe)
        self.processor.process_model(self.rotorcraft)

    def test_initialization(self):
        """Test that the Airframe component is initialized correctly."""
        self.assertIsInstance(self.airframe, Airframe)
        self.assertEqual(self.airframe.reference_system.name, 'airframe_1')
        self.assertEqual(self.airframe.reference_system.base_reference, 'global')
        self.assertEqual(self.airframe.reference_system.position_wrt_base.relative_position, [0., 0., 0.])
        self.assertEqual(self.airframe.reference_system.orientation_wrt_base.relative_position, l.eye())
        self.assertEqual(self.airframe.reference_system.velocity_wrt_base.relative_position, l.null())
        self.assertEqual(self.airframe.reference_system.angular_velocity_wrt_base.relative_position, l.null())

    def test_label_assignment(self):
        """Test that the processor correctly assigns labels."""
        # After processing, the airframe should have an mbdyn_label
        self.assertIsNotNone(self.airframe.reference_system.mbdyn_label)
        self.assertEqual(self.airframe.reference_system.mbdyn_label.name, 'AIRFRAME_1')

    def test_create_references(self):
        """Test the _create_references method."""
        references = self.airframe._create_references(self.processor)
        self.assertEqual(len(references), 1)
        ref = references[0]
        self.assertIsInstance(ref, l.Reference)
        # The reference should use the label assigned by the processor
        self.assertEqual(ref.idx, self.airframe.reference_system.mbdyn_label)
        self.assertEqual(ref.position, self.reference_system.position_wrt_base)
        self.assertEqual(ref.orientation, self.reference_system.orientation_wrt_base)
        self.assertEqual(ref.velocity, self.reference_system.velocity_wrt_base)
        self.assertEqual(ref.angular_velocity, self.reference_system.angular_velocity_wrt_base)

    def test_create_nodes(self):
        """Test the _create_nodes method."""
        nodes = self.airframe._create_nodes(self.processor)
        self.assertEqual(len(nodes), 1)
        node = nodes[0]
        self.assertIsInstance(node, l.DynamicNode)
        # The node should have an internal label created by the processor
        self.assertIsInstance(node.idx, l.MBVar)
        self.assertEqual(node.idx.name, 'AIRFRAME_1_NODE_1')
        self.assertEqual(node.position.relative_position, [0., 0., 0.])
        self.assertEqual(node.orientation.relative_position, l.eye())
        self.assertEqual(node.velocity.relative_position, l.null())
        self.assertEqual(node.angular_velocity.relative_position, l.null())

    def test_create_elements(self):
        """Test the _create_elements method."""
        elements = self.airframe._create_elements(self.processor)
        self.assertEqual(len(elements), 0)

    def test_generate_mbdyn_entities(self):
        """Test the complete entity generation workflow."""
        # Generate all entities
        self.airframe.generate_mbdyn_entities(self.processor)
        
        # Check that entities were created
        self.assertEqual(len(self.airframe._references), 1)
        self.assertEqual(len(self.airframe._nodes), 1)
        self.assertEqual(len(self.airframe._elements), 0)
        
        # Verify the reference
        ref = self.airframe._references[0]
        self.assertEqual(ref.idx, self.airframe.reference_system.mbdyn_label)
        
        # Verify the node
        node = self.airframe._nodes[0]
        self.assertEqual(node.idx.name, 'AIRFRAME_1_NODE_1')


if __name__ == '__main__':
    unittest.main()