import yt
import numpy as np
import matplotlib.pyplot as plt
import os.path


TRangePhase = {
                "all"      :[   0.0, 1.0e10],
                "hot"      :[1.0e+5, 1.0e10],
                "warm-hot" :[4.0e+4, 1.0e+5],
                "warm-cool":[1.0e+4, 4.0e+4],
                "cool"     :[   0.0, 1.0e+4],
                "cold"     :[   0.0, 5.0e+2],
              }


def calculate_average_star_formation_rate(ds, dt):

    ad = ds.all_data()

    if ds.dataset_type == 'gamer':
        mass          = ad[ 'new_star', 'ParMass'    ].in_units( 'Msun' )
        creation_time = ad[ 'new_star', 'ParCreTime' ].in_units( 'Myr' )

        # add back the feedback mass
        exploded = (ad[ 'new_star', 'ParSNIITime'  ] <=0)
        mass[exploded] += ds.quan( ds.parameters['FB_ResolvedSNeII_EjectMass'], 'code_mass' ).in_units('Msun')

    elif ds.dataset_type == 'gadget_hdf5':
        mass          = ad[ 'PartType4', 'Masses' ].in_units( 'Msun' )
        creation_time = ds.arr( ad[ 'PartType4', 'StellarFormationTime' ].d, 'code_time' ).in_units( 'Myr' )

    t_end   = ds.current_time.in_units( 'Myr' )
    t_start = t_end - dt

    sfr     = mass[( (creation_time > t_start) & (creation_time <= t_end) )].sum() / dt

    return sfr


def get_inoutflow_slab( ds, slab_z_kpc, slab_dz_kpc, phase, isIn=False ):

    # define the region
    cen                  = ds.domain_center
    slab_dx_kpc          = ds.domain_width[0].in_units('kpc').d
    slab_dy_kpc          = ds.domain_width[1].in_units('kpc').d
    upper_corner_L = cen + ds.arr( [0.0, 0.0, slab_z_kpc], 'kpc' ) - 0.5*ds.arr( [slab_dx_kpc, slab_dy_kpc, slab_dz_kpc], 'kpc' )
    upper_corner_R = cen + ds.arr( [0.0, 0.0, slab_z_kpc], 'kpc' ) + 0.5*ds.arr( [slab_dx_kpc, slab_dy_kpc, slab_dz_kpc], 'kpc' )
    lower_corner_L = cen - ds.arr( [0.0, 0.0, slab_z_kpc], 'kpc' ) - 0.5*ds.arr( [slab_dx_kpc, slab_dy_kpc, slab_dz_kpc], 'kpc' )
    lower_corner_R = cen - ds.arr( [0.0, 0.0, slab_z_kpc], 'kpc' ) + 0.5*ds.arr( [slab_dx_kpc, slab_dy_kpc, slab_dz_kpc], 'kpc' )

    if ds.dataset_type == 'gamer':

        FlowDir = "obj['gas', 'outflow_z_velocity'].d > 0.0" if not isIn else "obj['gas', 'outflow_z_velocity'].d < 0.0"
        upper_slab   = ds.box( upper_corner_L, upper_corner_R ).cut_region( [FlowDir] )
        lower_slab   = ds.box( lower_corner_L, lower_corner_R ).cut_region( [FlowDir] )
        GasPhase = "(obj['gas', 'T'].in_units('K').d > %.1e) & (obj['gas', 'T'].in_units('K').d < %.1e)"%(TRangePhase[phase][0],TRangePhase[phase][1])
        inoutflow_slab        = upper_slab.cut_region( [GasPhase] ) + lower_slab.cut_region( [GasPhase] )
        inoutflow_slab_p_type = 'gas'

    elif ds.dataset_type == 'gadget_hdf5':

        FlowDir = +1 if not isIn else -1
        def inoutflow_slab_particle( pfilter, data ):
            filter = ( (
                          (
                            (data[ 'gas', 'x' ] > upper_corner_L[0]) & (data[ 'gas', 'x' ] < upper_corner_R[0]) &
                            (data[ 'gas', 'y' ] > upper_corner_L[1]) & (data[ 'gas', 'y' ] < upper_corner_R[1]) &
                            (data[ 'gas', 'z' ] > upper_corner_L[2]) & (data[ 'gas', 'z' ] < upper_corner_R[2]) &
                            (FlowDir*data[ 'gas', 'outflow_z_velocity' ] > 0.0)
                          ) |
                          (
                            (data[ 'gas', 'x' ] > lower_corner_L[0]) & (data[ 'gas', 'x' ] < lower_corner_R[0]) &
                            (data[ 'gas', 'y' ] > lower_corner_L[1]) & (data[ 'gas', 'y' ] < lower_corner_R[1]) &
                            (data[ 'gas', 'z' ] > lower_corner_L[2]) & (data[ 'gas', 'z' ] < lower_corner_R[2]) &
                            (FlowDir*data[ 'gas', 'outflow_z_velocity' ] > 0.0)
                          )
                       ) &
                       (
                          (data['gas', 'T'].in_units('K').d > TRangePhase[phase][0]) &
                          (data['gas', 'T'].in_units('K').d < TRangePhase[phase][1]) &
                          (data['gas', 'density'].in_units('g/cm**3').d > 0.0)
                       )
                     )
            return filter
        p_type_name = 'outflow_slab_particle' if not isIn else 'inflow_slab_particle'
        yt.add_particle_filter( p_type_name, function=inoutflow_slab_particle, filtered_type='gas' )
        ds.add_particle_filter( p_type_name )
        inoutflow_slab        = ds.all_data()
        inoutflow_slab_p_type = p_type_name

    preFlow = 'out' if not isIn else 'in'
    suffix  = phase

    def plot_inoutflow_slab_slice(direction, field, width):
        s = yt.SlicePlot( ds, direction, field, center=ds.domain_center, width=width, buff_size=(1024, 1024), data_source=inoutflow_slab )
        s.set_axes_unit( 'kpc' )
        if field[1] == 'T':
            s.set_cmap( field, 'magma'  )
            s.set_zlim( field, 2.0e3, 2.0e5  )
        elif field[1] == 'density':
            s.set_cmap( field, 'viridis' )
            s.set_zlim( field, 1.0e-31, 1.0e-27 )
        s.annotate_timestamp( time_unit='Myr', corner='upper_right' )
        if direction == 'z':
            s.annotate_quiver(('gas','velocity_x'), ('gas','velocity_y'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
        if direction == 'x':
            s.annotate_quiver(('gas','velocity_y'), ('gas','velocity_z'), field_c=('gas','velocity_z'), factor=32, normalize=True, cmap='bwr_r', clim=(-5e6,5e6), alpha=0.7 )
        if ds.dataset_type != 'gadget_hdf5':
            s.annotate_grids( periodic=False )
        s.save( 'fig_%s_galactic_%sflow_slab_%s_gas_Slice_%s_%s.png'%(ds, preFlow, suffix, direction, field[1]), mpl_kwargs={'dpi':150} )

    def plot_inoutflow_slab_particle_projection(direction, field, width):
        if ds.dataset_type != 'gadget_hdf5':
            return
        p = yt.ParticleProjectionPlot( ds, direction, field, weight_field=(field[0],'ones'), center=ds.domain_center, width=width, depth=width )
        p.set_axes_unit( 'kpc' )
        if field[1] == 'T':
            p.set_cmap( field, 'magma' )
            p.set_zlim( field, 2.0e3, 2.0e5  )
        elif field[1] == 'density':
            p.set_cmap( field, 'viridis' )
            p.set_zlim( field, 1.0e-34, 1.0e-27 )
        p.annotate_timestamp( time_unit='Myr', corner='upper_right' )
        p.save( 'fig_%s_galactic_%sflow_slab_%s_gas_Projection_%s_%s.png'%(ds, preFlow, suffix, direction, field[1]), mpl_kwargs={'dpi':150} )


    # plot_inoutflow_slab_slice(              'x', (inoutflow_slab_p_type,'density'), (2.4*slab_z_kpc, 'kpc'))
    # plot_inoutflow_slab_particle_projection('x', (inoutflow_slab_p_type,'density'), (2.4*slab_z_kpc, 'kpc'))


    return inoutflow_slab, inoutflow_slab_p_type


def calculate_galactic_outflow_rate( ds, outflow_z_kpc, outflow_dz_kpc, outflow_phase ):

    # define the region
    outflow_slab, outflow_slab_p_type = get_inoutflow_slab( ds, outflow_z_kpc, outflow_dz_kpc, outflow_phase )

    # compute the outflow rate
    if ds.dataset_type == 'gamer':
        mass_outflow_rate   = ( outflow_slab.quantities.total_quantity( (outflow_slab_p_type,   'cell_mass_z_outflow_flux') ) / ds.quan( outflow_dz_kpc, 'kpc') ).in_units('Msun/yr')
        energy_outflow_rate = ( outflow_slab.quantities.total_quantity( (outflow_slab_p_type, 'cell_energy_z_outflow_flux') ) / ds.quan( outflow_dz_kpc, 'kpc') ).in_units('erg/s')
    elif ds.dataset_type == 'gadget_hdf5':
        mass_outflow_rate   = ( outflow_slab[outflow_slab_p_type,   'cell_mass_z_outflow_flux'].sum() / ds.quan( outflow_dz_kpc, 'kpc') ).in_units('Msun/yr')
        energy_outflow_rate = ( outflow_slab[outflow_slab_p_type, 'cell_energy_z_outflow_flux'].sum() / ds.quan( outflow_dz_kpc, 'kpc') ).in_units('erg/s')

    # get normalization
    E_SN                = ds.quan( 1.0e51, 'erg' )
    M_SN                = ds.quan( 100.0, 'Msun' )
    star_formation_rate = calculate_average_star_formation_rate(ds, ds.quan( 50.0, 'Myr' ) )

    # compute the loading factor
    mass_loading_factor   = mass_outflow_rate   / star_formation_rate
    energy_loading_factor = energy_outflow_rate / (star_formation_rate * E_SN / M_SN)

    # write to file
    filename_outflow_rate_table = '%s_Galactic_Outflow_Rate_%s_z_%2d_kpc'%(ds, outflow_phase, int(outflow_z_kpc))
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
               '  % 18.8e'%(energy_outflow_rate.in_units('erg/s').d) +\
               '  % 18.8e'%(star_formation_rate.in_units('Msun/yr').d) +\
               '  % 18.8e'%(mass_loading_factor) +\
               '  % 18.8e'%(energy_loading_factor) +\
               '\n'
    File.write( printstr )
    File.close()


def plot_galactic_outflow_rate_evolution(outflow_z_kpc, outflow_phase_list):

    ls_phase = {
                    "all"      :"r-",
                    "hot"      :"m--",
                    "warm-hot" :"y--",
                    "warm-cool":"g--",
                    "cool"     :"c--",
                    "cold"     :"b--",
                  }

    def plot_outflow_rate_phase(ax_or, phase):
        filename_outflow_rate_table = 'Galactic_Outflow_Rate_%s_z_%2d_kpc'%(phase, int(outflow_z_kpc))
        DataID, Time_cu, Time_Myr, MassOutflowRate, EnergyOutflowRate, StarFormationRate, MassLoadingFactor, EnergyLoadingFactor = np.loadtxt( filename_outflow_rate_table, skiprows=1, unpack=True )
        _, sorted_indices   = np.unique(DataID[::-1], return_index=True)
        DataID              = DataID             [::-1][sorted_indices]
        Time_cu             = Time_cu            [::-1][sorted_indices]
        Time_Myr            = Time_Myr           [::-1][sorted_indices]
        MassOutflowRate     = MassOutflowRate    [::-1][sorted_indices]
        EnergyOutflowRate   = EnergyOutflowRate  [::-1][sorted_indices]
        MassLoadingFactor   = MassLoadingFactor  [::-1][sorted_indices]
        EnergyLoadingFactor = EnergyLoadingFactor[::-1][sorted_indices]

        ax_or[0].plot( Time_Myr, MassOutflowRate, ls_phase[phase], label=phase, lw=1 )
        yr_in_s = 3.15569252e7
        ax_or[1].plot( Time_Myr, EnergyOutflowRate/1e51*yr_in_s, ls_phase[phase], lw=1 )

    def plot_outflow_rate():
        f_or, ax_or = plt.subplots( 2, 1, figsize=(6.4, 7.2) )

        for phase in outflow_phase_list:
            plot_outflow_rate_phase(ax_or, phase)

        # (1) mass outflow rate
        ax_or[0].set_yscale( 'log', nonpositive='clip' )
        ax_or[0].set_xlim(    0.0, 825.0 )
        ax_or[0].set_ylim( 1.0e-6, 5.0e-1 )
        ax_or[0].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax_or[0].set_ylabel( '$\dot{M}_\mathrm{10kpc}\mathrm{\ [M_{\odot}\ yr^{-1}]}$', fontsize='large' )
        ax_or[0].legend(loc=2, labelspacing=0.3)

        # (2) energy outflow rate
        ax_or[1].set_yscale( 'log', nonpositive='clip' )
        ax_or[1].set_xlim(     0.0, 825.0 )
        ax_or[1].set_ylim( 1.0e-11, 3.0e-4 )
        ax_or[1].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[i for i in range(0, 10, 1)]) )
        ax_or[1].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax_or[1].set_ylabel( '$\dot{E}_\mathrm{10kpc}\mathrm{\ [10^{51}erg\ yr^{-1}]}$', fontsize='large' )

        # save figure
        f_or.savefig( 'fig__galactic_outflow_rate_z_%2d_kpc.png'%(int(outflow_z_kpc)), bbox_inches='tight', pad_inches=0.05, dpi=150 )

    def plot_outflow_loading_factor_phase(ax_olf, phase):
        filename_outflow_rate_table = 'Galactic_Outflow_Rate_%s_z_%2d_kpc'%(phase, int(outflow_z_kpc))
        DataID, Time_cu, Time_Myr, MassOutflowRate, EnergyOutflowRate, StarFormationRate, MassLoadingFactor, EnergyLoadingFactor = np.loadtxt( filename_outflow_rate_table, skiprows=1, unpack=True )
        _, sorted_indices   = np.unique(DataID[::-1], return_index=True)
        DataID              = DataID             [::-1][sorted_indices]
        Time_cu             = Time_cu            [::-1][sorted_indices]
        Time_Myr            = Time_Myr           [::-1][sorted_indices]
        MassOutflowRate     = MassOutflowRate    [::-1][sorted_indices]
        EnergyOutflowRate   = EnergyOutflowRate  [::-1][sorted_indices]
        MassLoadingFactor   = MassLoadingFactor  [::-1][sorted_indices]
        EnergyLoadingFactor = EnergyLoadingFactor[::-1][sorted_indices]

        ax_olf[0].plot( Time_Myr, MassLoadingFactor,   ls_phase[phase], label=phase, lw=1 )
        ax_olf[1].plot( Time_Myr, EnergyLoadingFactor, ls_phase[phase], lw=1 )

    def plot_outflow_loading_factor():
        f_olf, ax_olf = plt.subplots( 2, 1, figsize=(6.4, 7.2) )

        for phase in outflow_phase_list:
            plot_outflow_loading_factor_phase(ax_olf, phase)

        # (1) mass loading factor
        ax_olf[0].set_yscale( 'log', nonpositive='clip' )
        ax_olf[0].set_xlim(    0.0,  825.0 )
        ax_olf[0].set_ylim( 1.0e-4,  5.0e1 )
        ax_olf[0].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax_olf[0].set_ylabel( '$\eta_\mathrm{m}$',      fontsize='large' )
        ax_olf[0].legend(loc=2, labelspacing=0.3)

        # (2) energy loading factor
        ax_olf[1].set_yscale( 'log', nonpositive='clip' )
        ax_olf[1].set_xlim(    0.0,  825.0 )
        ax_olf[1].set_ylim( 1.0e-6,  3.0e1 )
        ax_olf[1].yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[i for i in range(0, 10, 1)]) )
        ax_olf[1].set_xlabel( '$\mathrm{time\ [Myr]}$', fontsize='large' )
        ax_olf[1].set_ylabel( '$\eta_\mathrm{e}$',      fontsize='large' )

        # save figure
        f_olf.savefig( 'fig__galactic_outflow_loading_factor_z_%2d_kpc.png'%(int(outflow_z_kpc)), bbox_inches='tight', pad_inches=0.05, dpi=150 )

    plot_outflow_rate()
    plot_outflow_loading_factor()
