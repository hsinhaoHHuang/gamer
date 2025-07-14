import argparse
from matplotlib import patheffects
import sys
import yt
import matplotlib.pyplot as plt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-c', action='store', required=False, type=str, dest='code',
                     help='simulation code [%(default)s]', default='GAMER' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start    = args.idx_start
idx_end      = args.idx_end
didx         = args.didx
prefix       = args.prefix
code         = args.code

colormap     = {
                 'density'                 :'viridis',
                 'T'                       :'magma',
                 'kinetic_energy_density'  :'plasma',
                 'particle_density_on_grid':'algae',
                 'velocity_magnitude'      :'RdPu',
                 'particle'                :'algae',
                 'resolution_size'         :'cividis_r',
               }
field_unit   = {
                 'density'                 :'Msun/pc**3',
                 'T'                       :'K',
                 'kinetic_energy_density'  :'Msun/pc**3*km**2/s**2',
                 'particle_density_on_grid':'Msun/pc**3',
                 'velocity_magnitude'      :'km/s',
                 'resolution_size'         :'pc',
               }
zlim         = {
                 'density_s'                  :(1.0e-7, 1.0e+0),
                 'density_p'                  :(3.0e-2, 3.0e+1),
                 'T_s'                        :(1.0e+1, 1.0e+8),
                 'T_p'                        :(3.0e+0, 3.0e+6),
                 'kinetic_energy_density_s'   :(1.0e-3, 1.0e+2),
                 'kinetic_energy_density_p'   :(1.0e+2, 1.0e+6),
                 'particle_density_on_grid_s' :(1.0e-5, 1.0e+0),
                 'particle_density_on_grid_p' :(1.0e+0, 1.0e+4),
                 'velocity_magnitude_s'       :(1.0e+1, 1.5e+2),
                 'resolution_size_s'          :(1.0e+0, 1.0e+2),
                 'resolution_size_p'          :(1.0e+0, 1.0e+2),
                 ('all',      'particle_mass'):(1.0e+1, 1.0e+5),
                 ('Halo',     'particle_mass'):(1.0e+4, 1.0e+6),
                 ('Disk',     'particle_mass'):(3.0e+1, 1.0e+3),
                 ('new_star', 'particle_mass'):(3.0e-1, 3.0e+1),
                 ('exp_SNII', 'particle_mass'):(1.0e+1, 1.0e+3),
               }
zoomed_width = {
                 'a': 450.0,
                 'b':  30.0,
                 'c':  16.0,
                 'm':   6.0,
                 'n':   3.0,
               }
dpi          = 150

yt.enable_parallelism()


# load the dataset
if code == 'GAMER':
    ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
elif code == 'GIZMO':

    bbox = [ [-0.5*450.0, 0.5*450.0],
             [-0.5*450.0, 0.5*450.0],
             [-0.5*450.0, 0.5*450.0] ]

    unit_base = {
                  'UnitLength_in_cm'        : 3.085678e+21,
                  'UnitMass_in_g'           : 1.989e+43,
                  'UnitVelocity_in_cm_per_s': 100000,
                }
    ts = yt.DatasetSeries( [ prefix+'/snap_%03d.hdf5'%idx for idx in range(idx_start, idx_end+1, didx) ], unit_base=unit_base, bounding_box=bbox )
else:
    raise RuntimeError('Code %s is NOT supported  !!'%code)


# define the particle types
if code == 'GAMER':
    # the halo particles
    def Halo( pfilter, data ):
        filter = data[ 'all', 'ParType' ] == 2
        return filter
    yt.add_particle_filter( 'Halo', function=Halo, filtered_type='all', requires=['ParType'] )

    # the disk particles
    def Disk( pfilter, data ):
        filter = data[ 'all', 'ParType' ] == 3
        return filter
    yt.add_particle_filter( 'Disk', function=Disk, filtered_type='all', requires=['ParType'] )

    # the newly formed stars
    def new_star( pfilter, data ):
        filter = data[ 'all', 'ParCreTime' ] > 0
        return filter
    yt.add_particle_filter( 'new_star', function=new_star, filtered_type='all', requires=['ParCreTime'] )

    # the exploded SNe
    def exp_SNII( pfilter, data ):
        filter = data[ 'all', 'ParSNIITime' ] <= 0
        return filter
    yt.add_particle_filter( 'exp_SNII', function=exp_SNII, filtered_type='all', requires=['ParSNIITime'] )

    def young_star( pfilter, data ):
        filter = (data[ 'new_star', 'ParCreTime' ] > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
        return filter
    yt.add_particle_filter( 'young_star', function=young_star, filtered_type='new_star', requires=['ParCreTime'] )

    def young_SNII( pfilter, data ):
        filter = (-1.0*data[ 'exp_SNII', 'ParSNIITime' ]*data.ds.units.code_time > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
        return filter
    yt.add_particle_filter( 'young_SNII', function=young_SNII, filtered_type='exp_SNII', requires=['ParSNIITime'] )

elif code == 'GIZMO':
    def Halo( pfilter, data ):
        filter = data[ 'PartType1', 'ParticleIDs' ] > 0
        return filter
    yt.add_particle_filter( 'Halo', function=Halo, filtered_type='PartType1' )

    def Disk( pfilter, data ):
        filter = data[ 'PartType2', 'ParticleIDs' ] > 0
        return filter
    yt.add_particle_filter( 'Disk', function=Disk, filtered_type='PartType2' )

    def new_star( pfilter, data ):
        filter = data[ 'PartType4', 'ParticleIDs' ] > 0
        return filter
    yt.add_particle_filter( 'new_star', function=new_star, filtered_type='PartType4' )

    def young_star( pfilter, data ):
        filter = (data[ 'new_star', 'StellarFormationTime' ] > data.ds.current_time - data.ds.quan(2.5, 'Myr'))
        return filter
    yt.add_particle_filter( 'young_star', function=young_star, filtered_type='new_star', requires=['StellarFormationTime'] )


# add derived fields
def add_derived_fields(ds):

    if code == 'GIZMO':

        def _volume( field, data ):
           return data[('gas', 'mass')] / data[('gas', 'density')]
        ds.add_field( ('gas', 'volume'), function=_volume, sampling_type='particle', units='pc**3' )

    sampling_type = 'cell' if ds.dataset_type == 'gamer' else 'particle'

    field_u = ('gas', 'specific_thermal_energy') if ds.dataset_type == 'gamer' else ('PartType0', 'InternalEnergy')
    ds.mu   = 0.588235294117647  # Fully-ionized, 1/( 2.0*0.76/1 + 3.0*(1-0.76)/4 )

    def _T( field, data ):
       gamma = 5.0/3.0
       return data[field_u] * (gamma - 1.0) * data.ds.mu * data.ds.units.proton_mass / data.ds.units.boltzmann_constant
    ds.add_field( ('gas', 'T'), function=_T, sampling_type=sampling_type, units='K' )

    # define density^2 for calculating the weighted temperature
    def _density_square( field, data ):
        return data[('gas','density')]**2
    ds.add_field( ('gas', 'density_square'), function=_density_square, sampling_type=sampling_type, units='g**2/cm**6' )

    def _resolution_size( field, data ):
        return data[('gas','volume')]**(1./3.)
    ds.add_field( ('gas', 'resolution_size'), function=_resolution_size, sampling_type=sampling_type, units='pc' )

    # add the particle filter
    ds.add_particle_filter( 'Halo'     )
    ds.add_particle_filter( 'Disk'     )
    ds.add_particle_filter( 'new_star' )
    ds.add_particle_filter( 'young_star' )
    if code == 'GAMER':
        ds.add_particle_filter( 'exp_SNII' )
        ds.add_particle_filter( 'young_SNII' )


# main loop
for ds in ts.piter():

    add_derived_fields(ds)

    fields_list  = [ 'density', 'T', 'kinetic_energy_density', 'velocity_magnitude', 'resolution_size' ]
    pfields_list = [ ('all', 'particle_mass'), ('Halo', 'particle_mass'), ('Disk', 'particle_mass') ]

    if ('new_star', 'particle_mass') in ds.derived_field_list:
        pfields_list.append( ('new_star', 'particle_mass') )

    if code == 'GAMER':
        fields_list.append( 'particle_density_on_grid' )
        pfields_list.append( ('exp_SNII', 'particle_mass') )

    # decide the center
    center = ds.domain_center

    for zoom_mode in ['b', 'c', 'm', 'n']:

        # zoom in
        width_kpc = zoomed_width[zoom_mode]

        for direction in ['x', 'z']:
            for field in fields_list:

                # slices
                s = yt.SlicePlot( ds, direction, field, center=center, width=(width_kpc, 'kpc'), buff_size=(1024, 1024) )
                s.set_axes_unit( 'kpc' )
                s.set_unit( field, field_unit[field] )
                s.set_zlim( field, zlim[field+'_s'][0], zlim[field+'_s'][1] )
                s.set_cmap( field, colormap[field] )
                if field == 'velocity_magnitude':
                    s.set_log ( field, False )
                    if direction == 'z':
                        s.annotate_quiver('velocity_x', 'velocity_y', factor=16)
                s.annotate_timestamp( time_unit='Myr', corner='upper_right' )
                s.annotate_text( (0.02, 0.88), '%s'%(field), coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
                s.save( 'fig_%s_%s_Slice_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                if code == 'GAMER':
                    s.annotate_grids( periodic=False )
                    s.save( 'fig_%s_%s_Slice_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )


                if field == 'velocity_magnitude':
                    continue

                weight_field = 'density_square' if field == 'T' else None
                project_unit = '' if field == 'T' or field == 'resolution_size' else '*pc'
                proj_method  = 'min' if field == 'resolution_size' else 'integrate'


                # projections
                p = yt.ProjectionPlot( ds, direction, field, center=center, width=(width_kpc, 'kpc'), method=proj_method, weight_field=weight_field, buff_size=(1024, 1024) )
                p.set_axes_unit( 'kpc' )
                p.set_unit( field, field_unit[field]+project_unit )
                p.set_zlim( field, zlim[field+'_p'][0], zlim[field+'_p'][1] )
                p.set_cmap( field, colormap[field] )
                p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
                p.annotate_text( (0.02, 0.88), '%s'%(field), coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
                try:
                    p.save( 'fig_%s_%s_Projection_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                    if code == 'GAMER':
                        p.annotate_grids( periodic=False )
                        p.save( 'fig_%s_%s_Projection_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                        p.clear_annotations( index=-1 )
                    if field == 'density':
                        if ('new_star', 'particle_mass') in ds.derived_field_list:
                            p.annotate_particles( (width_kpc, 'kpc'), ptype='new_star', p_size=1, col='w', alpha=1.0, marker='o' )
                            p.save( 'fig_%s_%s_Projection_%s_%s_withStars.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                        if code == 'GAMER':
                            p.annotate_particles( (width_kpc, 'kpc'), ptype='exp_SNII', p_size=1, col='r', alpha=1.0, marker='o' )
                            p.save( 'fig_%s_%s_Projection_%s_%s_withSNe.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                            p.clear_annotations( index=-1 )
                        if ('young_star', 'particle_mass') in ds.derived_field_list:
                            p.clear_annotations( index=-1 )
                            p.annotate_particles( (width_kpc, 'kpc'), ptype='young_star', p_size=10, col='w', alpha=1.0, marker='o' )
                            p.save( 'fig_%s_%s_Projection_%s_%s_withYStars.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                        if code == 'GAMER':
                            p.annotate_particles( (width_kpc, 'kpc'), ptype='young_SNII', p_size=10, col='r', alpha=1.0, marker='o' )
                            p.save( 'fig_%s_%s_Projection_%s_%s_withYSNe.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
                except Exception as e:
                    print( e )
                    pass


            # particle plots
            for pfield in pfields_list:

                p = yt.ParticleProjectionPlot( ds, direction, pfield, center=center, width=(width_kpc, 'kpc'), depth=(width_kpc, 'kpc') )
                p.set_unit( pfield, 'Msun' )
                p.set_zlim( pfield, zlim[pfield][0], zlim[pfield][1] )
                p.set_cmap( pfield, colormap['particle'] )
                p.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
                p.save( 'fig_%s_%s_Particles_%s_%s.png'%(ds, zoom_mode, direction, pfield[0]), mpl_kwargs={'dpi':dpi} )

                p = yt.ParticlePlot( ds, (pfield[0], 'particle_velocity_'+direction ), (pfield[0], 'particle_velocity_y'), (pfield[0], 'particle_mass') )
                p.set_unit( pfield, 'Msun' )
                p.set_unit( (pfield[0], 'particle_velocity_'+direction), 'km/s' )
                p.set_unit( (pfield[0], 'particle_velocity_y'),          'km/s' )
                p.set_xlim( -100.0, 100.0 )
                p.set_ylim( -100.0, 100.0 )
                p.set_zlim( pfield, zlim[pfield][0], zlim[pfield][1] )
                p.set_cmap( pfield, colormap['particle'] )
                p.annotate_text( xpos=-96, ypos=76, text='%s\n%s'%(pfield), color='w', path_effects=[patheffects.withStroke(linewidth=2, foreground='k')] )
                p.save( 'fig_%s_%s_Particles_v%s_vy_%s_%s.png'%(ds, zoom_mode, direction, pfield[0], pfield[1]), mpl_kwargs={'dpi':dpi} )
