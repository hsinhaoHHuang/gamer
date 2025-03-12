import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas density-temperature phase diagram' )

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

colormap    = 'viridis'
width_kpc   = 6
nbin        = 300
dpi         = 150

x_lim_min   = 1.0e-30
x_lim_max   = 1.0e-19
y_lim_min   = 1.0e0
y_lim_max   = 1.0e8
z_lim_min   = 1.0e1
z_lim_max   = 1.0e8


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  define center as the location of peak gas density within 1 kpc from the center of gas mass
   v, cen1 = ds.find_max( ('gas', 'density') )
   sp1  = ds.sphere( cen1, (6.0, 'kpc') )
   cen2 = sp1.quantities.center_of_mass( use_gas=True, use_particles=False ).in_units( 'kpc' )
   sp2  = ds.sphere( cen2, (0.2, 'kpc') )
   cen3 = sp2.quantities.max_location( ('gas', 'density') )[1:]         # first value is not position
   cen  = cen3


#  only include the data within a sphere with a radius of width_kpc
   sp = ds.sphere( cen, (0.5*width_kpc, 'kpc') )


#  plot
   temp_dens = yt.PhasePlot( sp, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'),
                             weight_field=None, x_bins=nbin, y_bins=nbin )
   temp_dens.set_unit( 'cell_mass', 'Msun' )
   temp_dens.set_xlim( x_lim_min, x_lim_max )
   temp_dens.set_ylim( y_lim_min, y_lim_max )
   temp_dens.set_zlim( ('gas', 'cell_mass'), z_lim_min, z_lim_max )
   temp_dens.set_cmap( ('gas', 'cell_mass'), colormap )
   temp_dens.set_colorbar_label( ('gas', 'cell_mass'), "Mass ($\mathrm{M}_{\odot}$)" )
   temp_dens.annotate_text( xpos=x_lim_min*10**(0.80*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                            ypos=y_lim_min*10**(0.95*(np.log10(y_lim_max)-np.log10(y_lim_min))),
                            text='$t$ = {:.1f} {:s}'.format( ds.current_time.in_units('Myr').d, 'Myr' ),
                            color='black' )
   temp_dens.save( mpl_kwargs={'dpi':dpi} )


