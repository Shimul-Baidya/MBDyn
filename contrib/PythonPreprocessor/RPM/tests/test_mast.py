import unittest
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from RPM.components.mast import Mast
from RPM.components.hub import Hub
from RPM.components.rotor import Rotor
from RPM.components.airframe import Airframe
from RPM.components.base import Rotorcraft
from RPM.core.datatypes import ReferenceSystem, PhysicalQuantity
from RPM.core.processor import ModelProcessor


class TestMast(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures including processor and model."""
        # Create proper component hierarchy: Airframe -> Rotor -> Hub -> Mast
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
        
        self.mast = Mast(
            reference_system=ReferenceSystem(
                name='mast',
                base_reference='parent',
                position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference='parent'),
                orientation_wrt_base=l.Position(relative_position=l.eye(), reference='parent'),
                velocity_wrt_base=l.Position(relative_position=l.null(), reference='parent'),
                angular_velocity_wrt_base=l.Position(relative_position=l.null(), reference='parent'),
            )
        )
        
        self.airframe.add_sub_component(self.rotor)
        self.rotor.add_sub_component(self.hub)
        self.hub.add_sub_component(self.mast)
        
        self.processor = ModelProcessor()
        self.rotorcraft = Rotorcraft(name="TestRotorcraft")
        self.rotorcraft.add_root_component(self.airframe)
        self.processor.process_model(self.rotorcraft)

    def test_initialization(self):
        """Test that the Mast component is initialized correctly."""
        self.assertIsInstance(self.mast, Mast)
        self.assertEqual(self.mast.reference_system.name, 'mast')
        
    def test_label_assignment(self):
        """Test that processor assigns the correct label to Mast."""
        self.assertEqual(self.mast.reference_system.mbdyn_label.name, 'AIRFRAME_1_ROTOR_1_HUB_1_MAST_1')
        # MAST offset is 190, so: 10280 + 190 = 10470
        self.assertEqual(self.processor._value_registry['AIRFRAME_1_ROTOR_1_HUB_1_MAST_1'], 10470)
        
    def test_create_references(self):
        """Test that Mast creates exactly 1 reference."""
        references = self.mast._create_references(self.processor)
        self.assertEqual(len(references), 1)
        self.assertIsInstance(references[0], l.Reference)
        
    def test_reference_mbdyn_syntax(self):
        """Test that Mast reference generates correct MBDyn syntax."""
        references = self.mast._create_references(self.processor)
        mbdyn_output = str(references[0])
        self.assertIn('reference: AIRFRAME_1_ROTOR_1_HUB_1_MAST_1', mbdyn_output)
        self.assertTrue(mbdyn_output.endswith(';\n'))
        
    def test_create_nodes(self):
        """Test that Mast creates exactly 1 node."""
        nodes = self.mast._create_nodes(self.processor)
        self.assertEqual(len(nodes), 1)
        self.assertIsInstance(nodes[0], l.DynamicNode)
        
    def test_node_mbdyn_syntax(self):
        """Test that Mast node generates correct MBDyn syntax."""
        nodes = self.mast._create_nodes(self.processor)
        mbdyn_output = str(nodes[0])
        self.assertIn('structural: AIRFRAME_1_ROTOR_1_HUB_1_MAST_1_NODE_1', mbdyn_output)
        self.assertIn('dynamic', mbdyn_output)
        
    def test_complete_mbdyn_output(self):
        """Test complete MBDyn output for Mast."""
        references = self.mast._create_references(self.processor)
        nodes = self.mast._create_nodes(self.processor)
        
        expected_ref = (
            "reference: AIRFRAME_1_ROTOR_1_HUB_1_MAST_1, \n"
            "\treference, parent, 0.0, 0.0, 0.0,\n"
            "\treference, parent, eye,\n"
            "\treference, parent, null,\n"
            "\treference, parent, null;\n"
        )
        self.assertEqual(str(references[0]), expected_ref)
        
        expected_node = (
            "structural: AIRFRAME_1_ROTOR_1_HUB_1_MAST_1_NODE_1, dynamic,\n"
            "\treference, parent, 0.0, 0.0, 0.0,\n"
            "\treference, parent, eye,\n"
            "\treference, parent, null,\n"
            "\treference, parent, null;\n"
        )
        self.assertEqual(str(nodes[0]), expected_node)
        
    def test_create_elements(self):
        """Test that Mast creates 0 elements."""
        elements = self.mast._create_elements(self.processor)
        self.assertEqual(len(elements), 0)


if __name__ == '__main__':
    unittest.main()
