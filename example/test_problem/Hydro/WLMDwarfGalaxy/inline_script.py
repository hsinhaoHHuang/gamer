import yt_libyt
import yt
import numpy as np

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    field = 'density'

    half_width = 2.25
    box  = ds.box( ds.domain_center-ds.arr([half_width, half_width, half_width], 'code_length'),
                   ds.domain_center+ds.arr([half_width, half_width, half_width], 'code_length') )

    # Gas Projection
    pz_dens = yt.ProjectionPlot( ds, 'z', field, center='c',
                                 data_source=box, width=(2*half_width, 'code_length' ), buff_size=(1024, 1024) )
    pz_dens.set_axes_unit( 'kpc' )
    pz_dens.set_unit( field, "Msun/pc**2" )
    pz_dens.set_zlim( field, 3.0e-2, 3.0e+1 )
    pz_dens.set_cmap( field, 'viridis' )
    pz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )

    if yt.is_root():
        pz_dens.save()


    # Star and SN Projection
    def new_star( pfilter, data ):
        filter = data[ "all", "ParCreTime" ] > 0
        return filter

    yt.add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )
    ds.add_particle_filter( "new_star" )
    pz_dens.annotate_particles( (2*half_width, 'code_length'), ptype='new_star', p_size=0.7, col='w', marker='o' )
    pz_dens.save( '%s_wStars'%ds )

    def exp_SNII( pfilter, data ):
        filter = data[ "all", "ParSNIITime" ] <= 0
        return filter

    yt.add_particle_filter( "exp_SNII", function=exp_SNII, filtered_type="all", requires=["ParSNIITime"] )
    ds.add_particle_filter( "exp_SNII" )
    pz_dens.annotate_particles( (2*half_width, 'code_length'), ptype='exp_SNII', p_size=0.7, col='r', marker='o' )
    pz_dens.save( '%s_wSNeII'%ds )


    # Gas Phase
    x_lim_min = 9.0e-32
    x_lim_max = 1.0e-19
    y_lim_min = 1.0e0
    y_lim_max = 1.0e9
    z_lim_min = 1.0e1
    z_lim_max = 3.0e5
    temp_dens = yt.PhasePlot( box, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'),
                              weight_field=None, x_bins=300, y_bins=300 )
    temp_dens.set_unit( ('gas', 'cell_mass'), 'Msun' )
    temp_dens.set_xlim( x_lim_min, x_lim_max )
    temp_dens.set_ylim( y_lim_min, y_lim_max )
    temp_dens.set_zlim( ('gas', 'cell_mass'), z_lim_min, z_lim_max )
    temp_dens.set_cmap( ('gas', 'cell_mass'), 'viridis' )
    temp_dens.annotate_text( xpos=x_lim_min*10**(0.80*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                             ypos=y_lim_min*10**(0.95*(np.log10(y_lim_max)-np.log10(y_lim_min))),
                             text='$t$ = {:.1f} {:s}'.format( ds.current_time.in_units('Myr').d, 'Myr' ),
                             color='black' )
    temp_dens.save()


def yt_inline_inputArg( fields ):
    pass
