import yaml
import json
from typing import Dict
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import MBDynLib as l
from ..components.base import Rotorcraft, RotorcraftComponent

class ModelProcessor:
    """
    Processes a user-defined Rotorcraft model to generate MBDyn labels and entities.
    """
    def __init__(self):
        project_root = os.path.dirname(os.path.dirname(__file__))
        base_path = os.path.join(project_root, 'simulation', 'Labels')

        yml_path = base_path + '.yml'
        json_path = base_path + '.json'

        if os.path.exists(yml_path):
            with open(yml_path, 'r') as f:
                self.label_conventions = yaml.safe_load(f)
        elif os.path.exists(json_path):
            with open(json_path, 'r') as f:
                self.label_conventions = json.load(f)
        else:
            raise FileNotFoundError(f"Labels file not found. Looked for {yml_path} and {json_path}")
        
        self._label_counters: Dict[str, int] = {}
        self._value_registry: Dict[str, int] = {}
        self._mbvar_registry: Dict[str, l.MBVar] = {}

    def process_model(self, rotorcraft: Rotorcraft):
        """Public method to walk the entire rotorcraft model and assign labels.""" 
        for component in rotorcraft.root_components:
            self._process_component(component)

    def _process_component(self, component: RotorcraftComponent, parent_mbdyn_name: str = '', parent_value: int = 0):
        """A recursive helper method to process a component and its children."""
        component_type_name = component.reference_system.name

        count = self._label_counters.get(component_type_name, 0) + 1
        self._label_counters[component_type_name] = count

        # Extract base type for value lookup (e.g., 'rotor' from 'rotor_1')
        if '_' in component_type_name and component_type_name.split('_')[-1].isdigit():
            base_type = component_type_name.rsplit('_', 1)[0]  # 'rotor' from 'rotor_1'
        else:
            base_type = component_type_name

        base_name = component_type_name.upper()
        
        # Check if first instance and has digit suffix
        strip_digit_suffix = (count == 1 and '_' in component_type_name and 
                             component_type_name.split('_')[-1].isdigit())
        
        if parent_mbdyn_name == '':  # Root component
            if strip_digit_suffix:
                mbdyn_name = base_name  # Use as-is: 'ROTOR_1'
            else:
                mbdyn_name = f"{base_name}_{count}"  # Add counter: 'ROTOR_1_2'
        else:  # Sub-component
            if strip_digit_suffix:
                mbdyn_name = f"{parent_mbdyn_name}_{base_name}"  # AIRFRAME_1_ROTOR_1 (no extra counter)
            else:
                mbdyn_name = f"{parent_mbdyn_name}_{base_name}_{count}"  # AIRFRAME_1_HUB_1

        # Get conventions: try full name first, then base type
        conventions = self.label_conventions.get(component_type_name, {})
        if not conventions:
            # If no conventions for full name, try base type (e.g., 'rotor' instead of 'rotor_1')
            conventions = self.label_conventions.get(base_type, {})
        
        is_incremental = False

        base_value = conventions.get(f"{base_name}_{count}")
        if base_value is None:
            base_value = conventions.get(base_name)
        if base_value is None:
            base_value = conventions.get(f"CURR_{base_name}")
            is_incremental = True
        if base_value is None:
            # Try base type conventions
            if base_type != component_type_name:
                base_type_conventions = self.label_conventions.get(base_type, {})
                base_type_upper = base_type.upper()
                base_value = base_type_conventions.get(base_type_upper)
                if base_value is None:
                    base_value = base_type_conventions.get(f"CURR_{base_type_upper}")
                    if base_value is not None:
                        is_incremental = True
        if base_value is None:
            base_value = 0
            print(f"Warning: No label convention found for '{component_type_name}'. Defaulting to 0.")

        if is_incremental:
            final_value = parent_value + (base_value * count)
        else:
            if parent_mbdyn_name == '': # Root component (no parent)
                final_value = base_value
            else: # Sub-component (has parent)
                final_value = parent_value + base_value

        mbvar = l.MBVar(name=mbdyn_name, var_type='integer', expression=final_value)
        component.reference_system.mbdyn_label = mbvar
                
        self._mbvar_registry[mbdyn_name] = mbvar
        self._value_registry[mbdyn_name] = final_value

        for sub_component in component.sub_components:
            self._process_component(sub_component, parent_mbdyn_name=mbdyn_name, parent_value=final_value)

    def create_internal_label(self, owner: RotorcraftComponent, entity_type: str, counter: int) -> l.MBVar:
        """
        Creates a unique label for an internal entity (e.g., a node, or reference)
        owned by a component.
        
        Args:
            owner: The component that owns this entity
            entity_type: Type of entity - 'NODE', 'BODY', 'REF_FEATHERING', 'REF_NEUTRAL', etc.
            counter: Sequential counter for this entity type (starts from 0)
            
        Returns:
            MBVar label for this entity
            
        Notes:
            The legacy system uses specific offsets for different entity types:
            - Nodes, Bodies: component_base + counter (no offset)
            - Feathering axis references: component_base + 40000 + counter
            - Neutral axis references: component_base + 80000 + counter
        """
        owner_label_name = owner.reference_system.mbdyn_label.name
        owner_value = self._value_registry[owner_label_name]

        internal_name = f"{owner_label_name}_{entity_type.upper()}_{counter}"
        
        # Use offsets matching the legacy system
        # Nodes, bodies, and beam elements share the same label space (no offset)
        # Reference systems get dedicated ranges with large offsets
        type_offset = {
            'NODE': 0,           # Nodes: component_base + 0 + counter
            'BODY': 0,           # Bodies: component_base + 0 + counter (same as nodes)
            'BEAM': 0,           # Beams: component_base + 0 + counter
            'REF_FEATHERING': 40000,  # Feathering axis refs: component_base + 40000 + counter
            'REF_NEUTRAL': 80000,     # Neutral axis refs: component_base + 80000 + counter
            'REF_PRECONE': 1,         # Precone ref: component_base + 1 (single ref, not array)
        }.get(entity_type.upper(), 0)
        
        internal_value = owner_value + type_offset + counter

        if internal_name in self._mbvar_registry:
            return self._mbvar_registry[internal_name]

        mbvar = l.MBVar(name=internal_name, var_type='integer', expression=internal_value)
        self._mbvar_registry[internal_name] = mbvar
        self._value_registry[internal_name] = internal_value
        return mbvar

    def write_labels_set(self, output_path: str):
        """Writes the final labels.set file."""
        with open(output_path, 'w') as f:
            f.write("#beginpreprocess\n")
            for name, value in sorted(self._value_registry.items()):
                f.write(f"ConstMBVar('{name}', 'integer', {value})\n")
            f.write("#endpreprocess\n")
            
