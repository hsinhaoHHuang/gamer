import argparse
import sys
import yt
import plot_TempDens_Phase_and_PDF

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
parser.add_argument( '-c', action='store', required=False, type=str, dest='code',
                     help='simulation code [%(default)s]', default='GAMER' )

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
code        = args.code


width_kpc   = 6
nbin        = 300

x_lim_min   = 9.0e-32
x_lim_max   = 2.0e-18
y_lim_min   = 1.0e0
y_lim_max   = 1.0e9


yt.enable_parallelism()

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

for ds in ts.piter():

    if code == 'GIZMO':
        def _volume( field, data ):
           return data[('gas', 'mass')] / data[('gas', 'density')]
        ds.add_field( ('gas', 'volume'), function=_volume, sampling_type='particle', units='pc**3' )


#   define center as the location of peak gas density within 1 kpc from the center of gas mass
    v, cen1 = ds.find_max( ('gas', 'density') )
    sp1  = ds.sphere( cen1, (6.0, 'kpc') )
    cen2 = sp1.quantities.center_of_mass( use_gas=True, use_particles=False ).in_units( 'kpc' )
    cen  = cen2

#   only include the data within a sphere with a radius of width_kpc
    sp = ds.sphere( cen, (0.5*width_kpc, 'kpc') )

    plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, sp, '', 'gas', 'mass',   nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
    plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, sp, '', 'gas', 'volume', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
