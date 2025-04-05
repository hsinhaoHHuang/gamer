import yt_libyt
import yt

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    field = 'density'

    half_width = 3.0
    box  = ds.box( ds.domain_center-ds.arr([half_width, half_width, 0.5*half_width], 'code_length'),
                   ds.domain_center+ds.arr([half_width, half_width, 0.5*half_width], 'code_length') )

    pz_dens = yt.ProjectionPlot( ds, 'z', field, center='c',
                                 data_source=box, width=(2*half_width, 'code_length' )
    pz_dens.set_axes_unit( 'kpc' )
    pz_dens.set_unit( field, "Msun/pc**2" )
    pz_dens.set_zlim( field, 3.0e-2, 3.0e+1 )
    pz_dens.set_cmap( field, 'viridis' )
    pz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )

    if yt.is_root():
        pz_dens.save()


def yt_inline_inputArg( fields ):
    pass
