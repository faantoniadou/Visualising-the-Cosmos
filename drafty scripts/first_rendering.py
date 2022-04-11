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


zoom = float(input("Zoom factor: "))
ds = yt.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_101_covering_grid.h5")
    
def rendering1():
    
    # Create a volume rendering
    sc = yt.create_scene(ds, field=("grid", "nbody_mass"))#, lens_type="perspective")
    source = sc[0]
    source.set_field(("grid", "nbody_mass"))
    source.set_log(True)
    sc.camera.resolution = (1024, 1024)
    # Set the camera focus to a position that is offset from the center of the domain
    sc.camera.focus = ds.arr([0.3, 0.3, 0.3], "unitary")
    # Move the camera position to the other side of the dataset
    sc.camera.position = ds.arr([0, 0, 0], "unitary")
    sc.grey_opacity = True
    sc.camera.zoom(zoom)

    # sc.annotate_domain(ds, color=[1, 1, 1, 0.0002])
    # sc.annotate_axes(alpha=0.0002)

    tfh = TransferFunctionHelper(ds)
    tfh.set_field(("grid", "nbody_mass"))
    tfh.set_log(True)
    tfh.set_bounds()
    tfh.build_transfer_function()
    tfh.tf.add_layers(10, colormap='cmyt.arbre')

    # Grab the first render source and set it to use the new transfer function
    render_source = sc.get_source()
    render_source.transfer_function = tfh.tf
    # ad = ds.all_data()
    # mi, ma = ad.quantities.extrema(("grid", "nbody_mass"))

    # tfh = TransferFunctionHelper(ds)
    # tfh.set_field(("grid", "nbody_mass"))
    # tfh.set_log(True)
    # tfh.set_bounds()#(np.log10(mi+10), np.log10(ma)))
    # tfh.build_transfer_function()

    # # bounds = (10, np.log10(ma))

    # tf = yt.ColorTransferFunction((1, np.log10(ma)))
    # # source.tfh.bounds = bounds
    # tfh.tf.add_layers(10, colormap="cmyt.arbre")

    # # tfh.tf.map_to_colormap(mi+10, ma, colormap='cmyt.arbre')

    # source.tfh.tf = tf

    sc.save("customvzoom4.png", sigma_clip=4)
    # sc.save("customvzoom6.png", sigma_clip=6, render=False)
    # sc.save("customvzoom2.png", sigma_clip=2, render=False)


    # Grab the first render source and set it to use the new transfer function
    # render_source = sc.get_source()
    # render_source.transfer_function = tfh.tf


    # # Plot the transfer function, along with the CDF of the density field to
    # # see how the transfer function corresponds to structure in the CDF
    # source.tfh.plot(f"transfer_functionnbfsdfjkdsd.png", profile_field=("grid", "nbody_mass"))

    # text_string = f"T = {float(ds.current_time.to('Gyr'))} Gyr"

    # # save an annotated version of the volume rendering including a representation
    # # of the transfer function and a nice label showing the simulation time.
    # sc.save_annotated(
    #     "vol_annotatedsfd2e.png", sigma_clip=0, text_annotate=[[(0.1, 0.91), text_string]], render=False
    # )
    
rendering1()