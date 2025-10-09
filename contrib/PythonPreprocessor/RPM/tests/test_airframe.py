import unittest
import sys
import os

# Add the parent directory to the Python path to allow imports from there
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from RPM.components.airframe import Airframe
from RPM.core.datatypes import ReferenceSystem


try:
    import pydantic
except ImportError:
    pydantic = None


class TestAirframe(unittest.TestCase):
    def setUp(self):
        base_ref = 'global'
        self.reference_system = ReferenceSystem(
            name='AIRFRAME_1',
            base_reference=base_ref,
            position_wrt_base=l.Position(relative_position=[0., 0., 0.], reference=base_ref),
            orientation_wrt_base=l.Position(relative_position=[l.eye()], reference=base_ref),
            velocity_wrt_base=l.Position(relative_position=[l.null()], reference=base_ref),
            angular_velocity_wrt_base=l.Position(relative_position=[l.null()], reference=base_ref),
        )
        self.airframe = Airframe(reference_system=self.reference_system)

    def test_initialization(self):
        """Test that the Airframe component is initialized correctly."""
        self.assertIsInstance(self.airframe, Airframe)
        self.assertEqual(self.airframe.reference_system.name, 'AIRFRAME_1')
        self.assertEqual(self.airframe.reference_system.base_reference, 'global')
        self.assertEqual(self.airframe.reference_system.position_wrt_base.relative_position, [0., 0., 0.])
        self.assertEqual(self.airframe.reference_system.orientation_wrt_base.relative_position, [l.eye()])
        self.assertEqual(self.airframe.reference_system.velocity_wrt_base.relative_position, [l.null()])
        self.assertEqual(self.airframe.reference_system.angular_velocity_wrt_base.relative_position, [l.null()])

    def test_create_references(self):
        """Test the _create_references method."""
        references = self.airframe._create_references()
        self.assertEqual(len(references), 1)
        ref = references[0]
        self.assertIsInstance(ref, l.Reference)
        self.assertEqual(ref.idx, 20000)  # Using the hardcoded value for now
        self.assertEqual(ref.position, self.reference_system.position_wrt_base)
        self.assertEqual(ref.orientation, self.reference_system.orientation_wrt_base)
        self.assertEqual(ref.velocity, self.reference_system.velocity_wrt_base)
        self.assertEqual(ref.angular_velocity, self.reference_system.angular_velocity_wrt_base)

    def test_create_nodes(self):
        """Test the _create_nodes method."""
        nodes = self.airframe._create_nodes()
        self.assertEqual(len(nodes), 1)
        node = nodes[0]
        self.assertIsInstance(node, l.DynamicNode)
        self.assertEqual(node.idx, 20000)  # Using the hardcoded value for now
        self.assertEqual(node.position.relative_position, [0., 0., 0.])
        self.assertEqual(node.orientation.relative_position, [l.eye()])
        self.assertEqual(node.velocity.relative_position, [l.null()])
        self.assertEqual(node.angular_velocity.relative_position, [l.null()])

    def test_create_elements(self):
        """Test the _create_elements method."""
        elements = self.airframe._create_elements()
        self.assertEqual(len(elements), 0)


if __name__ == '__main__':
    unittest.main()