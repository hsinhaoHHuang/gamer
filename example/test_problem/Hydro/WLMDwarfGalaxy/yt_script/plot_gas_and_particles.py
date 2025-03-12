import argparse
import sys
import yt

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

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap     = {
                 'density'                 :'viridis',
                 'temperature'             :'plasma',
                 'particle_density_on_grid':'algae',
                 'velocity_magnitude'      :'RdPu',
               }
field_unit   = {
                 'density'                 :'Msun/pc**3',
                 'temperature'             :'K',
                 'particle_density_on_grid':'Msun/pc**3',
                 'velocity_magnitude'      :'km/s',
               }
zlim         = {
                 'density_s'                  :(1.0e-5, 3.0e-1),
                 'density_p'                  :(3.0e-2, 3.0e+1),
                 'temperature_s'              :(1.0e+2, 1.0e+6),
                 'temperature_p'              :(3.0e+3, 3.0e+4),
                 'particle_density_on_grid_s' :(1.0e-5, 1.0e+0),
                 'particle_density_on_grid_p' :(1.0e+0, 1.0e+4),
                 'velocity_magnitude_s'       :(1.0e+1, 1.0e+2),
                 ('all',      'particle_mass'):(1.0e+2, 1.0e+6),
                 ('new_star', 'particle_mass'):(1.0e+1, 1.0e+3),
                 ('expSNeII', 'particle_mass'):(1.0e+1, 1.0e+3),
               }
dpi          = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ "all", "ParCreTime" ] > 0
   return filter

yt.add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )

# define the particle filter for the exploded SNe
def expSNeII( pfilter, data ):
   filter = data[ "all", "ParSNIITime" ] <= 0
   return filter

yt.add_particle_filter( "expSNeII", function=expSNeII, filtered_type="all", requires=["ParSNIITime"] )

for ds in ts.piter():

   fields_list  = [ 'density', 'temperature', 'particle_density_on_grid', 'velocity_magnitude' ]
   pfields_list = [ ('all', 'particle_mass'), ('new_star', 'particle_mass'), ('expSNeII', 'particle_mass') ]

#  define density^2 for calculating the weighted temperature
   def _density_square( field, data ):
      return data["density"]**2

   ds.add_field( ('gas', 'density_square'), function=_density_square, sampling_type='cell', units='g**2/cm**6' )

#  add the particle filter
   ds.add_particle_filter( "new_star" )
   ds.add_particle_filter( "expSNeII" )

   for zoom_mode in ['a', 'b', 'c', 'm']:

      # decide the center
      center = ds.domain_center

      # zoom in
      if zoom_mode == 'a':
         width_kpc = 450
      elif zoom_mode == 'b':
         width_kpc = 30
      elif zoom_mode == 'c':
         width_kpc = 16
      elif zoom_mode == 'm':
         width_kpc = 6

      for direction in ['y', 'z']:
         for field in fields_list:

#           slice
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
            s.save( 'fig_%s_%s_Slice_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
            s.annotate_grids( periodic=False )
            s.save( 'fig_%s_%s_Slice_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )

            if field == 'velocity_magnitude':
               continue

            weight_field = 'density_square' if field == 'temperature' else None
            project_unit = '' if field == 'temperature' else '*pc'

#           projection
            p = yt.ProjectionPlot( ds, direction, field, center=center, width=(width_kpc, 'kpc'), method='integrate', weight_field=weight_field, buff_size=(1024, 1024) )
            p.set_axes_unit( 'kpc' )
            p.set_unit( field, field_unit[field]+project_unit )
            p.set_zlim( field, zlim[field+'_p'][0], zlim[field+'_p'][1] )
            p.set_cmap( field, colormap[field] )
            p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
            p.save( 'fig_%s_%s_Projection_%s_%s.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )
            p.annotate_grids( periodic=False )
            p.save( 'fig_%s_%s_Projection_%s_%s_withgrids.png'%(ds, zoom_mode, direction, field), mpl_kwargs={'dpi':dpi} )

         for pfield in pfields_list:
            p = yt.ParticleProjectionPlot( ds, 'z', pfield, center=center, width=(width_kpc, 'kpc'), depth=(width_kpc, 'kpc') )
            p.set_unit( pfield, 'Msun' )
            #p.set_zlim( pfield, zlim[pfield][0], zlim[pfield][1] )
            p.set_cmap( pfield, colormap['particle_density_on_grid'] )
            p.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
            p.save( 'fig_%s_%s_Particles_%s_%s.png'%(ds, zoom_mode, direction, pfield[0]), mpl_kwargs={'dpi':dpi} )
