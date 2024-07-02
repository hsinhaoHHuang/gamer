import yt_libyt
import yt
import gc

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    field = 'density'

    zoom = 2
    box  = ds.box( ds.domain_center-0.5/zoom*ds.domain_width,
                   ds.domain_center+0.5/zoom*ds.domain_width )

    pz_dens = yt.ProjectionPlot( ds, 'z', field,
                                 weight_field=('index','ones'),
                                 center='c',
                                 data_source=box,
                                 width=ds.domain_width[0]/zoom )
    pz_dens.set_axes_unit( 'kpc' )
    pz_dens.set_unit( field, "Msun/kpc**3" )
    pz_dens.set_zlim( field, 3.0e1, 3.0e5   )
    pz_dens.set_cmap( field, 'viridis' )
    pz_dens.set_background_color( field )
    pz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
    # pz_dens.annotate_grids()

    py_dens = yt.ProjectionPlot( ds, 'y', field,
                                 weight_field=('index','ones'),
                                 center='c',
                                 data_source=box,
                                 width=ds.domain_width[0]/zoom )
    py_dens.set_axes_unit( 'kpc' )
    py_dens.set_unit( field, "Msun/kpc**3" )
    py_dens.set_zlim( field, 3.0e1, 3.0e5   )
    py_dens.set_cmap( field, 'viridis' )
    py_dens.set_background_color( field )
    py_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
    # py_dens.annotate_grids()

    if yt.is_root():
        pz_dens.save()
        py_dens.save()

    del(ds)
    del(box)
    del(pz_dens)
    del(py_dens)
    gc.collect()


def yt_inline_inputArg( fields ):
    pass
