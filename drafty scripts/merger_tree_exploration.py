import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
import yt
import ytree
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context
from yt.utilities.cosmology import Cosmology

co = Cosmology()
print(co)
def basics():
    a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
    print (f'Box size: {a.box_size}')
    print (f'Cosmological parameters: {a.hubble_constant, a.omega_matter, a.omega_lambda}')
    print (f'Field list: {a.field_list}')
    print (f'Derived field list: {a.derived_field_list}')
    print (f'Index: {a["Snap_idx"]}')
    print (a.field_info['sam_Mvir'])

    my_tree = a[0]
    # AttributeError: 'TreeNode' object has no attribute 'field_list'
    #print (f'Node field list: {my_tree.field_list}')

    # how many trees are there?
    # print (f'Number of trees: {a.size}')

    # # Individual trees can be accessed by indexing the Arbor object.
    # print(f'Individual tree: {a[0]}')
    # print(f'Position: {a["position"]}')

    # # A TreeNode is one halo in a merger tree. The number is the universal identifier associated with halo. 
    # # It is unique to the whole arbor. Fields can be accessed for any given TreeNode in the same dictionary-like fashion.
    # print(f'Mass of the tree: {my_tree["mass"]}')

    # Array slicing can also be used to select multiple TreeNode objects. This will return a generator that can be iterated over or cast to a list.
    # every_second_tree = list(a[::2])
    # print (every_second_tree[0]["mass"])

    '''
    # this will not work
    # a[0].thing = 5
    # print (a[0].thing)
    # Traceback (most recent call last):
    #   File "<stdin>", line 1, in <module>
    # AttributeError: 'TreeNode' object has no attribute 'thing'
    # this will work
    # my_tree = a[0]
    # my_tree.thing = 5
    # print (my_tree.thing)
    # 5
    '''

    # Accessing All Nodes in a Tree

    # print (f'Access nodes: {my_tree["tree"]}')
    # # loop over nodes
    # for my_node in my_tree["tree"]:
    #     print (my_node, my_node["mass"])
        
    #print (f'All nodes in the tree: {list(my_tree["tree"])}')
    #print (f'Some nodes in the tree: {list(my_tree["tree"])[0:5]}')

    # This allows one to analyze the entire dataset using the full range of functionality provided by yt.
    ds = a.ytds

    # slc.set_buff_size(buff_size)
    # slc.set_background_color(("grid", "nbody_mass"))
    # slc.set_cmap(field=("grid", "nbody_mass"), cmap=color)

    # slc.save(f"{color}.png")

    #A list of defined vector fields can be seen by doing:
    # print(f'Vector fields: {a.field_info.vector_fields}')
    # print(f'Position: {a["position"]}')
    # print(f"Mass: {a['mass']}, \n radius: {a['virial_radius']}, \n velocity: {a['velocity']}")
    # print(f"IDs: {a['Tree_root_ID'], a['halo_id']}")
    vel = a['velocity'].to('unitary/s').value
    print(f"Velocity: {vel}")

basics()



#given a frame number, return the scale factor
def scale_factor(node):

    # main progenitor masses
    pmass = node["prog", "mass"]

    mh = 0.5 * node["mass"]
    m50 = pmass <= mh

    if not m50.any():
        ah = node["scale_factor"]
    else:
        pscale = node["prog", "scale_factor"]
        # linearly interpolate
        i = np.where(m50)[0][0]
        slope = (pscale[i-1] - pscale[i]) / (pmass[i-1] - pmass[i])
        ah = slope * (mh - pmass[i]) + pscale[i]

    #node["a50"] = ah
    return ah

a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
my_tree = a[0]
#print(f'Scale factor of first node: {scale_factor(my_tree)}')