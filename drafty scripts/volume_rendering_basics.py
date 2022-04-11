import yt
import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from tqdm import tqdm
import datetime
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name

fnumber = int(input('Grid number: '))
ds = yt.load(make_grid_name(fnumber))
sc = yt.create_scene(ds, field=("grid", "nbody_mass"), lens_type="perspective")

#Note that to render a different field, we would use pass the field name to yt.create_scene using the field argument.
#Now we can look at some information about the Scene we just created using the python print keyword:

print (sc)
#%%
#%%
#This prints out information about the Sources, Camera, and Lens associated with this Scene. Each of these can 
# also be printed individually. For example, to print only the information about the first (and currently, only) 
# Source, we can do:

print (sc.get_source())
 
# Set up a custom transfer function using the TransferFunctionHelper. 
# We use 10 Gaussians evenly spaced logarithmically between the min and max
# field values.

tfh = TransferFunctionHelper(ds)
tfh.set_field(("grid", "nbody_mass"))
tfh.set_log(True)
tfh.set_bounds()
tfh.build_transfer_function()
tfh.tf.add_layers(10)

yt.enable_parallelism()


cam = sc.camera
# save an image at the starting position
frame = 0
sc.save("camera_movement_%04i.png" % frame, sigma_clip = 4)
frame += 1

# Zoom in by a factor of 2 over 5 frames
for _ in cam.iter_zoom(2., 5):
    print('zoom now')
    sc.save("camera_movement`_%04i.png" % frame, sigma_clip = 4, render=False)
    frame += 1

# Move to the position [-10.0, 10.0, -10.0] over 5 frames
pos = ds.arr([-10.0, 10.0, -10.0], "code_length")
for _ in cam.iter_move(pos, 5):
    print('move now')
    sc.save("camera_movement2_%04i.png" % frame, sigma_clip = 4, render=False)
    frame += 1

# Rotate by 180 degrees over 5 frames
for _ in cam.iter_rotate(np.pi, 5):
    print('rotate now')
    sc.save("camera_movement3_%04i.png" % frame, sigma_clip = 4, render=False)
    frame += 1


# Grab the first render source and set it to use the new transfer function
# render_source = sc.get_source()
# render_source.transfer_function = tfh.tf

# sc.render()
# sc.save("{:%Y_%m_%d_%H_%M_%S}".format(today)+"_render", sigma_clip = 4)   # this line takes ages