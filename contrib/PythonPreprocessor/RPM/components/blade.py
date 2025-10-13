from typing import List, TYPE_CHECKING
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import MBDynLib as l
from .base import RotorcraftComponent
from ..core.datatypes import ReferenceSystem

if TYPE_CHECKING:
    from ..core.processor import ModelProcessor

# This is not final, it's only used for testing ModelProcessor
class Blade(RotorcraftComponent):
    """Represents a simple blade for testing."""
    reference_system: ReferenceSystem

    def _create_references(self, processor: 'ModelProcessor') -> List[l.Reference]: return []
    def _create_nodes(self, processor: 'ModelProcessor') -> List[l.Node]: return []
    def _create_elements(self, processor: 'ModelProcessor') -> List[l.Element]: return []
