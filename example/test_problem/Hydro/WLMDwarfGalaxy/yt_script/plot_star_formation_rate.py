import argparse
import sys
import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas density-temperature phase diagram' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-i', action='store', required=True,  type=int, dest='idx',
                     help='data index' )
parser.add_argument( '-c', action='store', required=False, type=str, dest='code',
                     help='simulation code [%(default)s]', default='GAMER' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx         = args.idx
prefix      = args.prefix
code        = args.code


fileout = 'fig__star_formation_rate'
nbin    = 100
dpi     = 150


# load data
if code == 'GAMER':

    ds = yt.load( prefix+'/Data_%06d'%idx )

    # define the particle filter for the newly formed stars
    def new_star( pfilter, data ):
       filter = data[ 'all', 'ParCreTime' ] > 0
       return filter
    add_particle_filter( 'new_star', function=new_star, filtered_type='all', requires=['ParCreTime'] )
    ds.add_particle_filter( 'new_star' )

    # get the mass and creation time of the new stars
    ad            = ds.all_data()
    mass          = ad[ 'new_star', 'ParMass'    ].in_units( 'Msun' )
    creation_time = ad[ 'new_star', 'ParCreTime' ].in_units( 'Myr' )

    # add back the feedback mass
    exploded      = (ad[ 'new_star', 'ParSNIITime'  ] <=0 )
    mass[exploded] += ds.quan( ds.parameters['FB_ResolvedSNeII_EjectMass'], 'code_mass' ).in_units('Msun')

elif code == 'GIZMO':

    bbox = [ [-0.5*450.0, 0.5*450.0],
             [-0.5*450.0, 0.5*450.0],
             [-0.5*450.0, 0.5*450.0] ]

    unit_base = {
                  'UnitLength_in_cm':         3.085678e+21,
                  'UnitMass_in_g':            1.989e+43,
                  'UnitVelocity_in_cm_per_s': 100000,
                }

    ds = yt.load( prefix+'/snap_%03d.hdf5'%idx, unit_base=unit_base, bounding_box=bbox )

    ad            = ds.all_data()
    mass          = ad[ 'PartType4', 'Masses' ].in_units( 'Msun' )
    creation_time = ds.arr( ad[ 'PartType4', 'StellarFormationTime' ].d, 'code_time' ).in_units( 'Myr' )

else:
    raise RuntimeError('Code %s is NOT supported  !!'%code)


# bin the data
t_start   = 0.0
t_end     = ds.current_time.in_units( 'Myr' )
t_bin     = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
upper_idx = np.digitize( creation_time, bins=t_bin, right=True )
time      = 0.5*( t_bin[:-1] + t_bin[1:] )

assert np.all( upper_idx > 0 ) and np.all( upper_idx < len(t_bin) ), 'incorrect upper_idx !!'


# calculate the star formation rate
Myr2yr = 1.0e6
sfr    = np.array(  [ mass[upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j])*Myr2yr )
                    for j in range(len(time)) ]  )
sfr[sfr == 0] = np.nan


# plot
plt.plot( time, sfr )
plt.yscale('log')
plt.xlim( 0.0, 1000 )
plt.ylim( 3.0e-5, 2.0e-2 )
plt.xlabel( '$\mathrm{t\ [Myr]}$',               fontsize='large' )
plt.ylabel( '$\mathrm{SFR\ [M_\odot yr^{-1}]}$', fontsize='large' )


# save figure
plt.savefig( fileout+'.png', bbox_inches='tight', pad_inches=0.05, dpi=dpi )
