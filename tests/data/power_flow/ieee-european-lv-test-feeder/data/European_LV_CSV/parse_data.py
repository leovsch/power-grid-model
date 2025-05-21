import pandas as pd
import numpy as np

from power_grid_model import initialize_array

from power_grid_model import LoadGenType, ComponentType, DatasetType

lines = pd.read_csv('Lines.csv', skiprows=1)
line_codes = pd.read_csv('LineCodes.csv', skiprows=1)
bus_coords = pd.read_csv('Buscoords.csv', skiprows=1)
loads = pd.read_csv('Loads.csv', skiprows=2)

n_nodes = bus_coords.shape[0] + 1
n_lines = lines.shape[0]
n_loads = loads.shape[0]

node = initialize_array(DatasetType.input, ComponentType.node, n_nodes)
source_node_id = 9999
node['id'] = np.array(bus_coords["Busname"]).tolist() + [source_node_id]
source_node_voltage = 11000
node["u_rated"] = np.array([230 * np.sqrt(3)] * (n_nodes - 1)).tolist() + [source_node_voltage]

start_line_id = max(node['id']) + 1
lines_with_linecode_df = lines.merge(line_codes, left_on='LineCode', right_on='Name')
line = initialize_array(DatasetType.input, ComponentType.line, n_lines)
line['id'] = np.array(range(start_line_id, start_line_id + n_lines))
line['from_node'] = np.array(lines["Bus1"])
line['to_node'] = np.array(lines["Bus2"])
line['r1'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["R1"] / 1000
line['x1'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["X1"] / 1000
# line['r0'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["R0"] / 1000
# line['x0'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["X0"] / 1000
line['c1'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["C1"] / 1000
# line['c0'] = lines_with_linecode_df["Length"] * lines_with_linecode_df["C0"] / 1000
line['from_status'] = np.array(1)
line['to_status'] = np.array(1)

start_load_id = max(line['id']) + 1
load = initialize_array(DatasetType.input, ComponentType.sym_load, n_loads)
load['id'] = np.array(range(start_load_id, start_load_id + n_loads))
load['node'] = np.array(loads["Bus"])
load['status'] = np.array(1)
load['type'] = np.array([LoadGenType.const_power] * n_loads)

source_id = max(load['id']) + 1
source = initialize_array(DatasetType.input, ComponentType.source, 1)
source['id'] = [source_id]
source['node'] = [source_node_id]
source['status'] = [1]
source['u_ref'] = 1.05
# TODO: add sk attriburte of source

transformer_id = source_id + 1
transformer = initialize_array(DatasetType.input, ComponentType.transformer, 1)
transformer['id'] = [transformer_id]
transformer['from_node'] = [source_node_id]
transformer['to_node'] = [1]
transformer["from_status"] = [1]
transformer["to_status"] = [1]
transformer["u1"] = [11000]
transformer["u2"] = [416]
transformer["sn"] = [0.8e6]
# transformer["uk"] = [0.1]
# transformer["pk"] = [1e3]
# transformer["i0"] = [1.0e-6]
# transformer["p0"] = [0.1]
# transformer["winding_from"] = [2]
# transformer["winding_to"] = [1]
# transformer["clock"] = [5]
# transformer["tap_side"] = [0]
# transformer["tap_pos"] = [3]
# transformer["tap_min"] = [-11]
# transformer["tap_max"] = [9]
# transformer["tap_size"] = [100]



print(node)
