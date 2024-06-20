# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

"""Data types for power grid model dataset and component types."""

# This file is automatically generated. DO NOT modify it manually!

import sys
from typing import Literal, Union

# To avoid conflicts with src/power_grid_model/enum.py
# pylint: disable=wrong-import-position

sys_path = sys.path.pop(0)
from enum import Enum

# pylint: enable=wrong-import-position

sys.path.append(sys_path)

# Value names are defined in lower case instead of upper case
# pylint: disable=invalid-name


class DataTypes(Enum):
    """Dataset types."""

    input = "input"

    sym_output = "sym_output"

    asym_output = "asym_output"

    update = "update"

    sc_output = "sc_output"


DataTypesLiteral = Literal[
    "input",
    "sym_output",
    "asym_output",
    "update",
    "sc_output",
]


DataType = Union[DataTypes, DataTypesLiteral]
"""
A DataType is the type of a :class:`BatchDataset`.

- Examples: 

    - DataType.input = "input"
    - DataType.update = "update"
"""


class ComponentTypes(Enum):
    """Grid component types."""

    node = "node"

    line = "line"

    link = "link"

    transformer = "transformer"

    transformer_tap_regulator = "transformer_tap_regulator"

    three_winding_transformer = "three_winding_transformer"

    sym_load = "sym_load"

    sym_gen = "sym_gen"

    asym_load = "asym_load"

    asym_gen = "asym_gen"

    shunt = "shunt"

    source = "source"

    sym_voltage_sensor = "sym_voltage_sensor"

    asym_voltage_sensor = "asym_voltage_sensor"

    sym_power_sensor = "sym_power_sensor"

    asym_power_sensor = "asym_power_sensor"

    fault = "fault"


ComponentTypesLiteral = Literal[
    "node",
    "line",
    "link",
    "transformer",
    "transformer_tap_regulator",
    "three_winding_transformer",
    "sym_load",
    "sym_gen",
    "asym_load",
    "asym_gen",
    "shunt",
    "source",
    "sym_voltage_sensor",
    "asym_voltage_sensor",
    "sym_power_sensor",
    "asym_power_sensor",
    "fault",
]


ComponentType = Union[ComponentTypes, ComponentTypesLiteral]
"""
A ComponentType is the type of a grid component.

- Examples: 

    - ComponentType.node = "node"
    - ComponentType.line = "line"
"""

# pylint: enable=invalid-name
