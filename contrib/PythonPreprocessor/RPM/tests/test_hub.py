import unittest
import sys
import os

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from RPM.components.hub import Hub
from RPM.components.rotor import Rotor
from RPM.components.airframe import Airframe
from RPM.components.base import Rotorcraft
from RPM.core.datatypes import ReferenceSystem, PhysicalQuantity
from RPM.core.processor import ModelProcessor


class TestHub(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures including processor and model."""
        # Create proper component hierarchy: Airframe -> Rotor -> Hub
        
        # Level 0: Airframe (root component)
        self.airframe = Airframe(
            reference_system=ReferenceSystem(
                name='airframe_1',
                base_reference='global',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference='global'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='global'),
                velocity_wrt_base=l.Position(relative_position=l.null(), reference='global'),
                angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference='global'),
            )
        )
        
        # Level 1: Rotor (sub-component of Airframe)
        self.rotor = Rotor(
            reference_system=ReferenceSystem(
                name='rotor_1',
                base_reference='AIRFRAME_1',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.116], reference='AIRFRAME_1'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='AIRFRAME_1'),
                velocity_wrt_base=l.Position(relative_position=l.null(), reference='AIRFRAME_1'),
                angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference='AIRFRAME_1'),
            ),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.250),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05569)
        )
        
        # Level 2: Hub (sub-component of Rotor)
        self.hub = Hub(
            reference_system=ReferenceSystem(
                name='hub',
                base_reference='parent',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference='parent'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='parent'),
                velocity_wrt_base=l.Position(relative_position=l.null(), reference='parent'),
                angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference='parent'),
            )
        )
        
        # Build hierarchy
        self.airframe.add_sub_component(self.rotor)
        self.rotor.add_sub_component(self.hub)
        
        # Create processor and process the model
        self.processor = ModelProcessor()
        self.rotorcraft = Rotorcraft(name="TestRotorcraft")
        self.rotorcraft.add_root_component(self.airframe)
        self.processor.process_model(self.rotorcraft)

    def test_initialization(self):
        """Test that the Hub component is initialized correctly."""
        self.assertIsInstance(self.hub, Hub)
        self.assertEqual(self.hub.reference_system.name, 'hub')
        self.assertEqual(self.hub.reference_system.base_reference, 'parent')
        
    def test_label_assignment(self):
        """Test that processor assigns the correct label to Hub."""
        # Verify full hierarchy in labels
        self.assertEqual(self.airframe.reference_system.mbdyn_label.name, 'AIRFRAME_1')
        self.assertEqual(self.rotor.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1')
        self.assertEqual(self.hub.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1_HUB_1')
        
        # Verify label values
        self.assertEqual(self.processor._value_registry['AIRFRAME_1'], 80)
        self.assertEqual(self.processor._value_registry['AIRFRAME_1_ROTOR_1'], 10080)  # 80 + 10000
        self.assertEqual(self.processor._value_registry['AIRFRAME_1_ROTOR_1_HUB_1'], 10280)  # 10080 + 200
        
    def test_create_references(self):
        """Test that Hub creates exactly 1 reference."""
        references = self.hub._create_references(self.processor)
        
        self.assertEqual(len(references), 1)
        self.assertIsInstance(references[0], l.Reference)
        self.assertEqual(references[0].idx, self.hub.reference_system.mbdyn_label)
        
    def test_reference_mbdyn_syntax(self):
        """Test that Hub reference generates correct MBDyn syntax."""
        references = self.hub._create_references(self.processor)
        ref = references[0]
        
        # Convert to MBDyn syntax string
        mbdyn_output = str(ref)
        
        # Expected format:
        # reference: AIRFRAME_1_ROTOR_1_HUB_1, 
        #     reference, parent, 0.0, 0.0, 0.0,
        #     reference, parent, eye,
        #     reference, parent, null,
        #     reference, parent, null;
        
        # Verify the output contains key elements
        self.assertIn('reference: AIRFRAME_1_ROTOR_1_HUB_1', mbdyn_output)
        self.assertIn('reference, parent, 0.0, 0.0, 0.0', mbdyn_output)
        self.assertIn('reference, parent, eye', mbdyn_output)
        self.assertIn('reference, parent, null', mbdyn_output)
        
        # Verify proper formatting with tabs and semicolon
        self.assertTrue(mbdyn_output.startswith('reference: '))
        self.assertTrue(mbdyn_output.endswith(';\n'))
        self.assertIn('\t', mbdyn_output)  # Should have tab indentation
        
        # Verify line structure
        lines = mbdyn_output.strip().split('\n')
        self.assertEqual(len(lines), 5)  # 5 lines: idx, position, orientation, velocity, angular_velocity
        
    def test_create_nodes(self):
        """Test that Hub creates exactly 1 node."""
        nodes = self.hub._create_nodes(self.processor)
        
        self.assertEqual(len(nodes), 1)
        self.assertIsInstance(nodes[0], l.DynamicNode)
        
        # Verify node is at origin of hub reference frame
        node = nodes[0]
        self.assertEqual(node.position.relative_position, [0., 0., 0.])
        self.assertEqual(node.position.reference, 'parent')
        
    def test_node_mbdyn_syntax(self):
        """Test that Hub node generates correct MBDyn syntax."""
        nodes = self.hub._create_nodes(self.processor)
        node = nodes[0]
        
        # Convert to MBDyn syntax string
        mbdyn_output = str(node)
        
        # Expected format:
        # structural: AIRFRAME_1_ROTOR_1_HUB_1_NODE_1, dynamic,
        #     reference, parent, 0.0, 0.0, 0.0,
        #     reference, parent, eye,
        #     reference, parent, null,
        #     reference, parent, null;
        
        # Verify the output contains key elements
        self.assertIn('structural: AIRFRAME_1_ROTOR_1_HUB_1_NODE_1', mbdyn_output)
        self.assertIn('dynamic', mbdyn_output)
        self.assertIn('reference, parent, 0.0, 0.0, 0.0', mbdyn_output)
        self.assertIn('reference, parent, eye', mbdyn_output)
        self.assertIn('reference, parent, null', mbdyn_output)
        
        # Verify proper formatting
        self.assertTrue(mbdyn_output.startswith('structural: '))
        self.assertTrue(mbdyn_output.endswith(';\n'))
        self.assertIn('\t', mbdyn_output)  # Should have tab indentation
        
        # Verify line structure (should have 5 lines for structural node)
        lines = mbdyn_output.strip().split('\n')
        self.assertEqual(len(lines), 5)  # 5 lines: idx+type, position, orientation, velocity, angular_velocity
        
        # Verify node type is dynamic
        self.assertIn('dynamic', lines[0])
        
        # Verify no accelerations keyword (since accelerations=None)
        self.assertNotIn('accelerations', mbdyn_output)
        
    def test_complete_mbdyn_output(self):
        """Test complete MBDyn output for Hub reference and node."""
        references = self.hub._create_references(self.processor)
        nodes = self.hub._create_nodes(self.processor)
        
        ref_output = str(references[0])
        node_output = str(nodes[0])
        
        # Expected Reference output
        expected_ref = (
            "reference: AIRFRAME_1_ROTOR_1_HUB_1, \n"
            "\treference, parent, 0.0, 0.0, 0.0,\n"
            "\treference, parent, eye,\n"
            "\treference, parent, null,\n"
            "\treference, parent, null;\n"
        )
        self.assertEqual(ref_output, expected_ref)
        
        # Expected Node output
        expected_node = (
            "structural: AIRFRAME_1_ROTOR_1_HUB_1_NODE_1, dynamic,\n"
            "\treference, parent, 0.0, 0.0, 0.0,\n"
            "\treference, parent, eye,\n"
            "\treference, parent, null,\n"
            "\treference, parent, null;\n"
        )
        self.assertEqual(node_output, expected_node)
        
    def test_create_elements(self):
        """Test that Hub creates 0 elements by default."""
        elements = self.hub._create_elements(self.processor)
        
        self.assertEqual(len(elements), 0)
        
    def test_generate_mbdyn_entities(self):
        """Test the full entity generation process."""
        # Generate all entities
        references = self.hub._create_references(self.processor)
        nodes = self.hub._create_nodes(self.processor)
        elements = self.hub._create_elements(self.processor)
        
        # Verify counts
        self.assertEqual(len(references), 1, "Hub should generate 1 reference")
        self.assertEqual(len(nodes), 1, "Hub should generate 1 node")
        self.assertEqual(len(elements), 0, "Hub should generate 0 elements by default")
        
    def test_collect_entities(self):
        """Test collecting entities from hub."""
        # Use processor to collect all entities
        references = []
        nodes = []
        elements = []
        
        # Collect from hub
        references.extend(self.hub._create_references(self.processor))
        nodes.extend(self.hub._create_nodes(self.processor))
        elements.extend(self.hub._create_elements(self.processor))
        
        self.assertEqual(len(references), 1, "Should collect 1 reference")
        self.assertEqual(len(nodes), 1, "Should collect 1 node")
        self.assertEqual(len(elements), 0, "Should collect 0 elements")     

    # TODO: Some labeling logic in processor.py is not correct for multiple instances of a component. Need to look into that and fix it
    def test_multiple_hubs_different_rotors(self):
        """Test creating hubs under different rotors (tiltrotor configuration)."""
        # Create second airframe for tiltrotor
        airframe2 = Airframe(
            reference_system=ReferenceSystem(
                name='airframe_2',
                base_reference='global',
                position_wrt_base=l.Position(relative_position=[0., -4., 0.], reference='global'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='global')
            )
        )
        
        # Create second rotor (using new naming convention)
        rotor2 = Rotor(
            reference_system=ReferenceSystem(
                name='rotor_2',
                base_reference='AIRFRAME_2',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.116], reference='AIRFRAME_2'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='AIRFRAME_2')
            ),
            n_blades=3,
            radius=PhysicalQuantity(unit='m', value=1.250),
            precone=PhysicalQuantity(unit='deg', value=2.75),
            precone_start=PhysicalQuantity(unit='adim', value=0.05569)
        )
        
        # Create second hub (using new naming convention)
        hub2 = Hub(
            reference_system=ReferenceSystem(
                name='hub',
                base_reference='parent',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference='parent')
            )
        )
        
        airframe2.add_sub_component(rotor2)
        rotor2.add_sub_component(hub2)
        
        # Create new processor and process both airframes (tiltrotor configuration)
        processor2 = ModelProcessor()
        rotorcraft2 = Rotorcraft(name="TiltRotor")
        rotorcraft2.add_root_component(self.airframe)
        rotorcraft2.add_root_component(airframe2)
        processor2.process_model(rotorcraft2)
        
        # First rotor hub: AIRFRAME_1_ROTOR_1_HUB_1 = 10280
        self.assertEqual(processor2._value_registry['AIRFRAME_1_ROTOR_1_HUB_1'], 10280)
        # Second rotor hub: AIRFRAME_2_ROTOR_2_HUB_2 (counters are global, not per-parent)
        # AIRFRAME_2 = 90, ROTOR_2 = 90 + 20000 = 20090, HUB_2 = 20090 + 200 = 20290
        self.assertEqual(processor2._value_registry['AIRFRAME_2_ROTOR_2_HUB_2'], 20290)

if __name__ == '__main__':
    unittest.main()
