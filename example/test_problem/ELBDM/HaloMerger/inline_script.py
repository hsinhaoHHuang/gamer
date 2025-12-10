import yt_libyt
import yt

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    field = 'density'

    pz_dens = yt.ProjectionPlot( ds, 'z', field, center='c' )
    pz_dens.set_axes_unit( 'kpc' )
    pz_dens.set_unit( field, "Msun/kpc**2" )
    pz_dens.set_zlim( field, 1.0e3, 1.0e8   )
    pz_dens.set_cmap( field, 'viridis' )
    pz_dens.set_background_color( field )
    pz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
    # pz_dens.annotate_grids()

    if yt.is_root():
        pz_dens.save()


def yt_inline_inputArg( fields ):
    pass
