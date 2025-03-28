from abc import ABC
from typing import List, Optional, Annotated
from MBDynLib import *

if imported_pydantic:
    from pydantic import BaseModel, ConfigDict, field_validator, validate_call, Field
else:
    class BaseModel:
        """Placeholder for when pydantic is not available"""
        def __init__(self, *args, **kwargs):
            if len(args) > 0:
                raise TypeError(
                    'MBDyn models cannot be initialized using positional arguments')
            for key, value in kwargs.items():
                setattr(self, key, value)

class MBDynModel(MBEntity):
    """
    Main class for holding all blocks of an MBDyn model and generating the complete input file.
    
    Required blocks:
        - data
        - problem (Initial Value)
        - control_data
        - nodes
        - elements
    Optional blocks:
        - drivers
    """
    
    data: Data
    problem: InitialValue  
    control_data: ControlData
    nodes: Annotated[List, Field(arbitrary_type_allowed=True)] #TODO: Replace with proper typing once migration to Node2 is complete
    drivers: Optional[List[FileDriver]] = []
    elements: Annotated[List, Field(arbitrary_type_allowed=True)] #TODO: Replace with proper typing once migration to Element2 is complete
    
    def add_node(self, node: Node) -> None:
        self.nodes.append(node)

    def add_driver(self, driver: FileDriver) -> None:
        self.drivers.append(driver)

    def add_element(self, element: Element) -> None:
        self.elements.append(element)

    def __str__(self) -> str:
        """Generate complete MBDyn input file content"""
        output = []

        # Data block
        output.append(str(self.data))
        output.append("")  # Add blank line

        # Problem block  
        output.append(str(self.problem))
        output.append("")

        # Control data block
        output.append(str(self.control_data))
        output.append("")

        # Nodes block
        output.append("begin: nodes;")
        for node in self.nodes:
            output.append(str(node))
        output.append("end: nodes;\n")

        # Drivers block (optional)
        if self.drivers:
            output.append("begin: drivers;")
            for driver in self.drivers:
                output.append(str(driver))
            output.append("end: drivers;\n")

        # Elements block
        output.append("begin: elements;")
        for element in self.elements:
            output.append(str(element))
        output.append("end: elements;")

        return "\n".join(output)
