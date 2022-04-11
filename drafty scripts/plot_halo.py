import yt
from yt.units import kpc
import ytree
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera
import unyt
from unyt.unit_registry import default_unit_registry
from unyt import Mpc, kpc
from matplotlib import rc_context
import datetime
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
    
from yt.utilities.nodal_data_utils import get_nodal_data
from yt.units.yt_array import uconcatenate
import numpy as np
from yt.frontends.ytdata.api import YTGridDataset
#yt.enable_parallelism()



# this makes the file name for the covering grids directory
# def make_grid_name(fnumber):
#     grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
#     return grid_name

# # this makes file name for halo trees 
# def make_halo_name(fnumber):
#     halo_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar_{fnumber:04}.h5"
#     return halo_name

# n=int(input('Tree number: '))
#color = str(input("color: "))
# loop through h5 files
# get today's date to save snapshots with different names

today = datetime.datetime.now()

a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
trees = list(a[:])
ds = a.ytds


def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name


# sphere = ds.sphere(ds.domain_center, (5, "Mpc"))
# #print (sphere["halos", "mass"])

# for node in a.get_nodes_from_selection(sphere):
#     print (node["position"])

# ap = ytree.AnalysisPipeline(output_dir="my_analysis")



# for tree in trees:
#     for node in tree["forest"]:
#         ap.process_target(node)

# def get_filename_from_redshift(redshift):
#     for i in range(101):
#         snap_name = make_grid_name(i)
#         ds = yt.load(snap_name)
#         ds_redshift = ds.current_redshift

#         if ds_redshift == redshift:
#             return snap_name


# def say_hello(node):
#     print (f"This is node {node}! I will now be analyzed.")
        
# def calculate_gas_mass(node):
#     sphere = node.sphere
#     node["gas_mass"] = sphere.quantities.total_quantity(("gas", "mass"))

# def print_field_value(node, field, units=None):
#     val = node[field]
#     if units is not None:
#         val.convert_to_units(units)
#     print (f"Value of {field} for node {node} is {val}.")
    
# def minimum_mass(node, value):
#     return node["mass"] >= value

# def delete_attributes(node, attributes):
#     for attr in attributes:
#         if hasattr(node, attr):
#             delattr(node, attr)
            
# def get_yt_dataset(node):
#     # assume you have something like this
#     filename = file_path
#     # attach it to the node for later use
#     node.ds = yt.load(filename)

# def get_yt_sphere(node):
#     # this works if get_yt_dataset has been called first
#     ds = node.ds

#     center = node["position"].to("unitary")
#     radius = node["virial_radius"].to("unitary")
#     node.sphere = ds.sphere((center, "unitary"), (radius, "unitary"))
    
# def calculate_mass(node):
#     sphere = node.sphere
#     node["mass"] = sphere.quantities.total_quantity(("gas", "mass"))

# def gas_mass_recipe(pipeline):
#     pipeline.add_operation(get_yt_dataset)
#     pipeline.add_operation(get_yt_sphere)
#     pipeline.add_operation(calculate_gas_mass)
#     pipeline.add_operation(delete_attributes, ["ds", "sphere"])
    
# ap = ytree.AnalysisPipeline(output_dir="my_analysis")
# ap.add_recipe(gas_mass_recipe)
 

# for node in ytree.parallel_nodes(trees):
#     ap.process_target(node)

def ancestors():
    def max_value(ancestors, field):
        vals = np.array([a[field] for a in ancestors])
        return ancestors[np.argmax(vals)]

    ytree.add_tree_node_selector("max_field_value", max_value)

    a.set_selector("max_field_value", "mass")
    my_tree = a[0]
    print (list(my_tree["prog", "redshift"]))