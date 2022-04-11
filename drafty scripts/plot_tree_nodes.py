import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
import yt
from yt import YTArray
import unyt
import ytree
from unyt import Mpc, kpc
from matplotlib import rc_context
import graphviz

def my_node(halo):
    prog = list(halo.find_root()['prog', 'uid'])
    if halo['uid'] in prog:
        color = 'red'
    else:
        color = 'black'

    label = \
    """
    id: %d
    mass: %.2e Msun
    """ % (halo['uid'], halo['mass'].to('Msun'))

    my_kwargs = {"label": label, "fontsize": 8,
                 "shape": "square", "color": color}
    return my_kwargs

def main():
    tree_number = int(input('Tree number: '))
    tree = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar{tree_number:04}.h5"
    
    p = ytree.TreePlot(tree, dot_kwargs={'rankdir': "BT"},
                       node_function=my_node)
    p.save('tree_custom_node.png')

    return my_node(halo)
    
main()