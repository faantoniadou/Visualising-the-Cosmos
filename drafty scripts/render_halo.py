import yt
import ytree
import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper \
    import TransferFunctionHelper
import datetime
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#yt.enable_parallelism()
today = datetime.datetime.now()

# quantity maximum halo to find:
quantity = str(input('Maximum quantity to search for in halo catalogue: '))

how_large = int(input(f'Order to largest {quantity} (e.g. 1 is maximum and 2 is second to largest): '))
zoom_factor = float(input("Width: "))
no_layers = int(input("Number of colour layers: "))

# Load the dataset
ds = yt.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_099_covering_grid.h5")

# Create a volume rendering
sc = yt.create_scene(ds, field=("grid", "nbody_mass"))
sc.grey_opacity = True

source = sc[0]
cam = sc.camera


# source.transfer_function = yt.ColorTransferFunction(
#     np.log10((1e-30, 1e-23)) )
# def linramp(vals, minval, maxval):
#     return (vals - vals.min()) / (vals.max() - vals.min())
# source.transfer_function.map_to_colormap(
#     np.log10(1e-25), np.log10(8e-24), colormap="cmyt.arbre", scale_func=linramp
# )

# For this low resolution dataset it's very important to use interpolated
# vertex centered data to avoid artifacts. For high resolution data this
# setting may cause a substantial slowdown for marginal visual improvement.
source.set_use_ghost_zones(True)


# HERE WE ATTEMPT TO FIND THE MOST MASSIVE HALO 
a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")

# find most massive halo
 # find halo with the order of maximum of the quantity we want
values = (a[quantity]).tolist()

values.sort()
i = np.where(a[quantity] == (values[-how_large]))[0][0]           #['mass'])

my_tree = a[i]

print(f'Halo position: {my_tree["position"]}')         # gives position now
print(f'Halo mass: {my_tree["mass"].to("Msun")}')         
print(f'Halo radius: {my_tree["virial_radius"].to("kpc")}')         
print((ds.domain_right_edge - ds.domain_left_edge)[0].to("Mpc"))

#print(f'Position of most massive progenitor: {my_tree["prog", "position"]}') # give position of most massive progenitor

#print(f'Redshifts of progenitors: {my_tree["prog", "redshift"]}')

halo_position = np.array(my_tree["position"])
radius = my_tree["virial_radius"].to("unitary")

print(ds.quan(radius * zoom_factor, 'unitary').to('Mpc'))
cam.width = ds.arr(zoom_factor, 'Mpc') #radius * zoom_factor
cam.focus = halo_position

# Set up a custom transfer function using the TransferFunctionHelper.
# We use 10 Gaussians evenly spaced logarithmically between the min and max
# field values.
ad = ds.all_data()
mi, ma = ad.quantities.extrema(("grid", "nbody_mass"))

tfh = TransferFunctionHelper(ds)
tfh.set_field(("grid", "nbody_mass"))
tfh.set_log(True)
tfh.set_bounds()

# tfh.set_bounds((np.log10(mi+10), np.log10(ma)))
tfh.build_transfer_function()

bounds = (10, np.log10(ma))

tf = yt.ColorTransferFunction((np.log10(mi+10), np.log10(ma)))
# source.tfh.bounds = bounds
tfh.tf.add_layers(no_layers, colormap="cmyt.arbre")

# tfh.tf.map_to_colormap(mi+10, ma, colormap='cmyt.arbre')

source.tfh.tf = tf


render_source = sc.get_source()
render_source.transfer_function = tfh.tf

sc.camera.resolution = (1024, 1024)


source.tfh.plot(f"transfer_function{how_large}_{no_layers}.png", profile_field=("grid", "nbody_mass"))
sc.save(f"position_zoom={zoom_factor}_order={how_large}_colours={no_layers}.png", sigma_clip=6.) 
# Draw the grid boundaries
# sc.annotate_grids(ds, alpha=0.01)

# Save the image with the keyword.
# sc.annotate_domain(ds, color=[1, 1, 1, 0.01], alpha=0.01)
# sc.save(f"{ds}_vr_domain_order={how_large}_colours={no_layers}.png", sigma_clip=6.)




#ani.save("{:%Y_%m_%d_%H_%M_%S}".format(today)+"_anim.mp4", writer=Writer)
