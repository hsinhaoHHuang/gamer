import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
import matplotlib.pyplot as plt
from matplotlib import patheffects
import argparse
import sys
import os.path
import gc
import plot_TempDens_Phase_and_PDF

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the galactic outflow rate' )

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
filename_outflow_rate_table = 'Galactic_Outflow_Rate'

dpi         = 150

nbin        = 300
x_lim_min   = 5.0e-33
x_lim_max   = 1.0e-25
y_lim_min   = 5.0e1
y_lim_max   = 5.0e6

plane_z_kpc =  10
delta_z_kpc =   0.1
Delta_z_kpc =  15.0
delta_x_kpc = 450.0
delta_y_kpc = 450.0


def calculate_average_star_formation_rate(ds, dt):

    if ds.dataset_type == 'gamer':
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
        exploded = (ad[ 'new_star', 'ParSNIITime'  ] <=0)
        mass[exploded] += ds.quan( ds.parameters['FB_ResolvedSNeII_EjectMass'], 'code_mass' ).in_units('Msun')

    elif ds.dataset_type == 'gadget_hdf5':
        ad            = ds.all_data()
        mass          = ad[ 'PartType4', 'Masses' ].in_units( 'Msun' )
        creation_time = ds.arr( ad[ 'PartType4', 'StellarFormationTime' ].d, 'code_time' ).in_units( 'Myr' )

    t_end   = ds.current_time.in_units( 'Myr' )
    t_start = t_end - dt

    sfr     = mass[( (creation_time > t_start) & (creation_time <= t_end) )].sum() / dt

    return sfr


def calculate_galactic_outflow_rate_mesh( ds, u_cL, u_cR, l_cL, l_cR ):

    upper_slab   = ds.box( u_cL, u_cR ).cut_region( ["obj['gas', 'velocity_z'].d > 0.0"] )
    lower_slab   = ds.box( l_cL, l_cR ).cut_region( ["obj['gas', 'velocity_z'].d < 0.0"] )

    outflow_slab = upper_slab + lower_slab

    mor = outflow_slab.quantities.total_quantity( ('gas',   'mass_outflow_rate') ).in_units('Msun/yr')
    eor = outflow_slab.quantities.total_quantity( ('gas', 'energy_outflow_rate') ).in_units('erg/yr')

    upper_slab.clear_data()
    lower_slab.clear_data()
    outflow_slab.clear_data()
    gc.collect()

    return mor, eor


def calculate_galactic_outflow_rate_particle( ds, u_cL, u_cR, l_cL, l_cR ):

    def outflow_slab_particle( pfilter, data ):
        filter = ( (
                     (data[ 'gas', 'x' ] > u_cL[0]) & (data[ 'gas', 'x' ] < u_cR[0]) &
                     (data[ 'gas', 'y' ] > u_cL[1]) & (data[ 'gas', 'y' ] < u_cR[1]) &
                     (data[ 'gas', 'z' ] > u_cL[2]) & (data[ 'gas', 'z' ] < u_cR[2]) &
                     (data[ 'gas', 'velocity_z' ] > 0.0)
                   ) |
                   (
                     (data[ 'gas', 'x' ] > l_cL[0]) & (data[ 'gas', 'x' ] < l_cR[0]) &
                     (data[ 'gas', 'y' ] > l_cL[1]) & (data[ 'gas', 'y' ] < l_cR[1]) &
                     (data[ 'gas', 'z' ] > l_cL[2]) & (data[ 'gas', 'z' ] < l_cR[2]) &
                     (data[ 'gas', 'velocity_z' ] < 0.0)
                   )
                 )
        return filter
    yt.add_particle_filter( 'outflow_slab_particle', function=outflow_slab_particle, filtered_type='gas' )
    ds.add_particle_filter( 'outflow_slab_particle' )

    ad = ds.all_data()

    mor = ad['outflow_slab_particle',   'mass_outflow_rate'].in_units('Msun/yr').sum()
    eor = ad['outflow_slab_particle', 'energy_outflow_rate'].in_units('erg/yr').sum()

    ad.clear_data()
    gc.collect()

    return mor, eor


def calculate_galactic_outflow_rate( ds, z_kpc, dx_kpc, dy_kpc, dz_kpc ):

    # define the region
    cen                  = ds.domain_center
    upper_slab_corner_L  = cen + ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) - 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    upper_slab_corner_R  = cen + ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) + 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    lower_slab_corner_L  = cen - ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) - 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    lower_slab_corner_R  = cen - ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) + 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )

    # add derived fields
    sampling_type = 'cell' if ds.dataset_type == 'gamer' else 'particle'

    def _mass_outflow_rate( field, data ):
        return np.abs( data[('gas', 'mass')] * data[('gas', 'velocity_z')] / data.ds.arr( dz_kpc, 'kpc' ) )
    ds.add_field( ('gas', 'mass_outflow_rate'), function=_mass_outflow_rate, sampling_type=sampling_type, units='Msun/yr' )

    def _energy_outflow_rate( field, data ):
        GAMMA = 5.0/3.0
        return np.abs( data[('gas', 'mass')]*( data[('gas', 'velocity_magnitude')]**2 + GAMMA*data[('gas', 'specific_thermal_energy')] )*data[('gas', 'velocity_z')]/data.ds.arr( delta_z_kpc, 'kpc' ) )
    ds.add_field( ('gas', 'energy_outflow_rate'), function=_energy_outflow_rate, sampling_type=sampling_type, units='erg/yr' )

    # compute the outflow rate
    if sampling_type == 'cell':
        mass_outflow_rate, energy_outflow_rate = calculate_galactic_outflow_rate_mesh(     ds, upper_slab_corner_L, upper_slab_corner_R, lower_slab_corner_L, lower_slab_corner_R )
    elif sampling_type == 'particle':
        mass_outflow_rate, energy_outflow_rate = calculate_galactic_outflow_rate_particle( ds, upper_slab_corner_L, upper_slab_corner_R, lower_slab_corner_L, lower_slab_corner_R )

    # compute the loading factor
    E_SN                  = ds.quan( 1.0e51, 'erg' )
    M_SN                  = ds.quan( 100.0, 'Msun' )
    star_formation_rate   = calculate_average_star_formation_rate(ds, ds.quan( 50.0, 'Myr' ) )
    mass_loading_factor   = mass_outflow_rate/ star_formation_rate
    energy_loading_factor = (energy_outflow_rate * M_SN) / (star_formation_rate * E_SN)

    # write to file
    if os.path.exists(filename_outflow_rate_table) == False:
        # create the header
        File = open( filename_outflow_rate_table, 'a' )
        headerstr = '# % 7s'%('DataID') +\
                    '  % 18s'%('Time_cu') +\
                    '  % 18s'%('Time_Myr') +\
                    '  % 18s'%('MassOutflowRate') +\
                    '  % 18s'%('EnergyOutflowRate') +\
                    '  % 18s'%('StarFormationRate') +\
                    '  % 18s'%('MassLoadingFactor') +\
                    '  % 18s'%('EnergyLoadingFactor') +\
                    '\n'
        File.write( headerstr )
        File.close()

    File = open( filename_outflow_rate_table, 'a' )
    printstr = '  % 7d'%(int(str(ds)[5:11])) +\
               '  % 18.8e'%(ds.current_time.in_units('code_time').d) +\
               '  % 18.8e'%(ds.current_time.in_units('Myr').d) +\
               '  % 18.8e'%(mass_outflow_rate.in_units('Msun/yr').d) +\
               '  % 18.8e'%(energy_outflow_rate.in_units('erg/yr').d) +\
               '  % 18.8e'%(star_formation_rate.in_units('Msun/yr').d) +\
               '  % 18.8e'%(mass_loading_factor) +\
               '  % 18.8e'%(energy_loading_factor) +\
               '\n'
    File.write( printstr )
    File.close()


def get_outflow_region( ds, z_kpc, dx_kpc, dy_kpc, dz_kpc ):

    # define the region
    cen            = ds.domain_center
    upper_corner_L = cen + ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) - 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    upper_corner_R = cen + ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) + 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    lower_corner_L = cen - ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) - 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )
    lower_corner_R = cen - ds.arr( [0.0, 0.0, z_kpc], 'kpc' ) + 0.5*ds.arr( [dx_kpc, dy_kpc, dz_kpc], 'kpc' )

    if ds.dataset_type == 'gamer':

        upper_region   = ds.box( upper_corner_L, upper_corner_R ).cut_region( ["obj['gas', 'velocity_z'].d > 0.0"] )
        lower_region   = ds.box( lower_corner_L, lower_corner_R ).cut_region( ["obj['gas', 'velocity_z'].d < 0.0"] )
        outflow_region = upper_region + lower_region

    elif ds.dataset_type == 'gadget_hdf5':

        def outflow_region_particle( pfilter, data ):
            filter = ( (
                         (data[ 'gas', 'x' ] > upper_corner_L[0]) & (data[ 'gas', 'x' ] < upper_corner_R[0]) &
                         (data[ 'gas', 'y' ] > upper_corner_L[1]) & (data[ 'gas', 'y' ] < upper_corner_R[1]) &
                         (data[ 'gas', 'z' ] > upper_corner_L[2]) & (data[ 'gas', 'z' ] < upper_corner_R[2]) &
                         (data[ 'gas', 'velocity_z' ] > 0.0)
                       ) |
                       (
                         (data[ 'gas', 'x' ] > lower_corner_L[0]) & (data[ 'gas', 'x' ] < lower_corner_R[0]) &
                         (data[ 'gas', 'y' ] > lower_corner_L[1]) & (data[ 'gas', 'y' ] < lower_corner_R[1]) &
                         (data[ 'gas', 'z' ] > lower_corner_L[2]) & (data[ 'gas', 'z' ] < lower_corner_R[2]) &
                         (data[ 'gas', 'velocity_z' ] < 0.0)
                       )
                     )
            return filter
        yt.add_particle_filter( 'outflow_region_particle', function=outflow_region_particle, filtered_type='gas' )
        ds.add_particle_filter( 'outflow_region_particle' )
        outflow_region = ds.all_data()

    return outflow_region


def plot_min_temp_location(ds, Min_temp_location, field_type):

    if Min_temp_location[0].in_units('K').d < 1.0e3:
       for direction in ['x', 'z']:
          zoom_out_max = 12 if direction == 'x' else 0
          for field in [(field_type,'density'), (field_type,'T')]:
              for zoom_out in range(0, zoom_out_max+1, 1):
                 s = yt.SlicePlot( ds, direction, field, center=Min_temp_location[1:], width=(3.0*(1.25**zoom_out), 'kpc'), buff_size=(1024, 1024) )
                 s.set_axes_unit( 'kpc' )
                 cmap = 'viridis' if field[1] == 'density' else 'magma'
                 lim_min = x_lim_min if field[1] == 'density' else 3.0e2
                 lim_max = x_lim_max if field[1] == 'density' else 1.0e3
                 s.set_cmap( field, cmap )
                 s.set_zlim( field, lim_min, lim_max*(2.0**zoom_out)  )
                 s.annotate_timestamp( time_unit='Myr', corner='upper_right' )
                 if direction == 'z':
                    s.annotate_quiver(('gas','velocity_x'), ('gas','velocity_y'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
                 if direction == 'x':
                    s.annotate_quiver(('gas','velocity_y'), ('gas','velocity_z'), field_c=('gas','velocity_x'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
                 s.annotate_text( (0.02, 0.08), r'$T_\mathrm{min}$ = %.1f K, $\rho$ = %.1e g cm$^{-3}$'%(Min_temp_location[0].in_units('K').d, ds.point(Min_temp_location[1:])[(field_type,'density')].in_units('g/cm**3').d[0])+
                                                 '\n at [ %.2f, %.2f, %.2f ] kpc'%((Min_temp_location[1]-ds.domain_center[0]).in_units('kpc').d,
                                                                                   (Min_temp_location[2]-ds.domain_center[1]).in_units('kpc').d,
                                                                                   (Min_temp_location[3]-ds.domain_center[2]).in_units('kpc').d),
                                  coord_system='axis', text_args={'color':'w', 'path_effects':[patheffects.withStroke(linewidth=2, foreground='k')]} )
                 s.save( 'fig_%s_galactic_outflow_gas_Slice_%s_%s_%02d.png'%(ds, direction, field[1], zoom_out), mpl_kwargs={'dpi':dpi} )


def plot_evolution():
    DataID, Time_cu, Time_Myr, MassOutflowRate, EnergyOutflowRate, StarFormationRate, MassLoadingFactor, EnergyLoadingFactor = np.loadtxt( filename_outflow_rate_table, skiprows=1, unpack=True )
    _, sorted_indices   = np.unique(DataID[::-1], return_index=True)
    DataID              = DataID             [::-1][sorted_indices]
    Time_cu             = Time_cu            [::-1][sorted_indices]
    Time_Myr            = Time_Myr           [::-1][sorted_indices]
    MassOutflowRate     = MassOutflowRate    [::-1][sorted_indices]
    EnergyOutflowRate   = EnergyOutflowRate  [::-1][sorted_indices]
    MassLoadingFactor   = MassLoadingFactor  [::-1][sorted_indices]
    EnergyLoadingFactor = EnergyLoadingFactor[::-1][sorted_indices]

    def plot_outflow_rate():
        f, ax = plt.subplots( 2, 1, figsize=(6.4, 7.2) )

        # (1) mass outflow rate
        ax[0].plot( Time_Myr, MassOutflowRate, 'r-', lw=1 )
        ax[0].set_yscale( 'log', nonpositive='clip' )
        ax[0].set_xlim(    0.0, 1000.0 )
        ax[0].set_ylim( 1.0e-6, 5.0e-1 )
        ax[0].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax[0].set_ylabel( '$\dot{M}_\mathrm{10kpc}\mathrm{\ [M_{\odot}\ yr^{-1}]}$', fontsize='large' )

        # (2) energy outflow rate
        ax[1].plot( Time_Myr, EnergyOutflowRate/1e51, 'r-', lw=1 )
        ax[1].set_yscale( 'log', nonpositive='clip' )
        ax[1].set_xlim(     0.0, 1000.0 )
        ax[1].set_ylim( 1.0e-11, 3.0e-4 )
        ax[1].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[i for i in range(0, 10, 1)]) )
        ax[1].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax[1].set_ylabel( '$\dot{E}_\mathrm{10kpc}\mathrm{\ [10^{51}erg\ yr^{-1}]}$', fontsize='large' )

        # save figure
        plt.savefig( 'fig__galactic_outflow_rate.png', bbox_inches='tight', pad_inches=0.05, dpi=dpi )

    def plot_outflow_loading_factor():
        f, ax = plt.subplots( 2, 1, figsize=(6.4, 7.2) )

        # (1) mass loading factor
        ax[0].plot( Time_Myr, MassLoadingFactor, 'r-', lw=1 )
        ax[0].set_yscale( 'log', nonpositive='clip' )
        ax[0].set_xlim(    0.0, 1000.0 )
        ax[0].set_ylim( 1.0e-4,  5.0e1 )
        ax[0].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax[0].set_ylabel( '$\eta_{m}$', fontsize='large' )

        # (2) energy loading factor
        ax[1].plot( Time_Myr, EnergyLoadingFactor, 'r-', lw=1 )
        ax[1].set_yscale( 'log', nonpositive='clip' )
        ax[1].set_xlim(    0.0, 1000.0 )
        ax[1].set_ylim( 1.0e-6,  3.0e1 )
        ax[1].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[i for i in range(0, 10, 1)]) )
        ax[1].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax[1].set_ylabel( '$\eta_{e}$', fontsize='large' )

        # save figure
        plt.savefig( 'fig__galactic_outflow_loading_factor.png', bbox_inches='tight', pad_inches=0.05, dpi=dpi )

    plot_outflow_rate()
    plot_outflow_loading_factor()


# load data
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


# loop over all datasets
for ds in ts.piter():

    # define additional fields
    if code == 'GIZMO':
        def _specific_thermal_energy( field, data ):
           return data[('PartType0', 'specific_thermal_energy')]
        ds.add_field( ('gas', 'specific_thermal_energy'), function=_specific_thermal_energy, sampling_type='particle', units='km**2/s**2' )

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

    # main part
    calculate_galactic_outflow_rate( ds, plane_z_kpc, delta_x_kpc, delta_y_kpc, delta_z_kpc )

    # define a larger outflow region
    outflow_region = get_outflow_region( ds, plane_z_kpc, delta_x_kpc, delta_y_kpc, Delta_z_kpc )
    outflow_p_type = 'outflow_region_particle' if ds.dataset_type == 'gadget_hdf5' else 'gas'

    if not outflow_region.quantities.total_quantity( (outflow_p_type, 'volume' ) ) > 0.0:
       continue

    # find the minimum temperature
    if ds.dataset_type == 'gamer':
        Min_temp_location = outflow_region.quantities.min_location( (outflow_p_type, 'T') )
        plot_min_temp_location(ds, Min_temp_location, outflow_p_type)
    elif ds.dataset_type == 'gadget_hdf5':
        idx = outflow_region[(outflow_p_type, 'T')].argmin()
        Min_temp_location = [ outflow_region[(outflow_p_type, 'T')][idx],
                              outflow_region[(outflow_p_type, 'x')][idx],
                              outflow_region[(outflow_p_type, 'y')][idx],
                              outflow_region[(outflow_p_type, 'z')][idx] ]

    # plot the phase diagram
    plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region, '_outflow', outflow_p_type, 'mass',   nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
    plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region, '_outflow', outflow_p_type, 'volume', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )

    # plot different phases separately
    if ds.dataset_type == 'gamer':
        # hot phase
        outflow_region_hot  = outflow_region.cut_region( ["obj['gas', 'T'].in_units('K') > 5.0e5"] )
        if outflow_region_hot.quantities.total_quantity( (outflow_p_type, 'volume' ) ) > 0.0:
            plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region_hot,  '_outflow_hot', outflow_p_type,  'mass',   nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
            plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region_hot,  '_outflow_hot', outflow_p_type,  'volume', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
        outflow_region_hot.clear_data()
        gc.collect()

        # cool phase
        outflow_region_cool = outflow_region.cut_region( ["obj['gas', 'T'].in_units('K') < 2.0e4"] )
        if outflow_region_cool.quantities.total_quantity( (outflow_p_type, 'volume' ) ) > 0.0:
            plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region_cool, '_outflow_cool', outflow_p_type, 'mass',   nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
            plot_TempDens_Phase_and_PDF.plot_PhaseDiagram( ds, outflow_region_cool, '_outflow_cool', outflow_p_type, 'volume', nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max )
        outflow_region_cool.clear_data()
        gc.collect()

    outflow_region.clear_data()
    gc.collect()

plot_evolution()
