import yt
import numpy as np
import matplotlib.pyplot as plt

ds0 = yt.load( '../Data_000000' )

Step, Time, CoreDensity, CoreRadius, CoreMass, CoreRadiusCoreMass = np.loadtxt( '../Record__CoreProperties', skiprows=1, unpack=True )

# assign the units
Time               = ds0.arr( Time,               'code_time'             )
CoreDensity        = ds0.arr( CoreDensity,        'code_density'          )
CoreRadius         = ds0.arr( CoreRadius,         'code_length'           )
CoreMass           = ds0.arr( CoreMass,           'code_mass'             )
CoreRadiusCoreMass = ds0.arr( CoreRadiusCoreMass, 'code_length*code_mass' )

step_end = 4648

###################################################################################################################
# density profile
def Soliton_rho_c_from_r_c(m22, r_c):
    return 1.945e7/( m22**2 * r_c**4 )

def Soliton_fitting_analytical_dens(r, m22, r_c):
    rho_c = Soliton_rho_c_from_r_c(m22, r_c)
    x = r / r_c
    return rho_c*( 1.0+9.06e-2*(x**2) )**(-8)
###################################################################################################################

for step in range(0, step_end+1, 1):

      # load the reference profiles
      GAMER_DensProf_r, GAMER_DensProf_dens, GAMER_DensProf_Weight, GAMER_DesnProf_Cells = np.loadtxt( '../DensityProfile_Step_%09d'%step, skiprows=1, unpack=True )
      #YT_DensProf_r, YT_DensProf_dens                               = np.loadtxt( '../plot/Data_%06d_halo#1_density_profile_au'%step, skiprows=1, unpack=True )

      # assign the units
      GAMER_DensProf_r    = ds0.arr( GAMER_DensProf_r,    'code_length'  )
      GAMER_DensProf_dens = ds0.arr( GAMER_DensProf_dens, 'code_density' )

      # parameters for the soliton
      soliton_fitting_m22   = 0.2
      soliton_fitting_r_c   = CoreRadius[(Step==step)][0]
      soliton_fitting_rho_c = ds0.quan( Soliton_rho_c_from_r_c(soliton_fitting_m22, soliton_fitting_r_c.in_units('kpc').d), 'Msun/kpc**3' )

      soliton_densprof_radius  = ds0.arr( np.logspace( np.log10(0.5*np.min(GAMER_DensProf_r.in_units('kpc').d)), np.log10(2.0*np.max(GAMER_DensProf_r.in_units('kpc').d)), num=256 ), 'kpc' )
      soliton_densprof_density = ds0.arr( Soliton_fitting_analytical_dens( soliton_densprof_radius.in_units('kpc').d, soliton_fitting_m22, soliton_fitting_r_c.in_units('kpc').d ), 'Msun/kpc**3' )

      # decide the units for plotting
      #UNIT_L_PLOT = 'code_length'
      #UNIT_D_PLOT = 'code_density'
      UNIT_L_PLOT = 'kpc'
      UNIT_D_PLOT = 'Msun/kpc**3'
      UNIT_M_PLOT = 'Msun'

      # create the figure
      fig = plt.figure()
      ax  = fig.add_subplot(111)

      # plot the profiles
      ax.plot( GAMER_DensProf_r.in_units(UNIT_L_PLOT).d,        GAMER_DensProf_dens.in_units(UNIT_D_PLOT).d,      '.-',  color='b',  label='Simulation'   )
      ax.plot( [0.01*np.min(GAMER_DensProf_r.in_units(UNIT_L_PLOT).d), 100.0*np.max(GAMER_DensProf_r.in_units(UNIT_L_PLOT).d)], [CoreDensity[(Step==step)][0].in_units(UNIT_D_PLOT).d,     CoreDensity[(Step==step)][0].in_units(UNIT_D_PLOT).d      ],  ':',  color='xkcd:sky blue',      label=r'$\rho_{\rm peak}$', zorder=1  )

      # plot the density profile
      ax.plot( soliton_densprof_radius.in_units(UNIT_L_PLOT).d, soliton_densprof_density.in_units(UNIT_D_PLOT).d, '--',  color='r',  label='soliton' )

      # plot some important values for reference
      ax.plot( [0.01*np.min(GAMER_DensProf_r.in_units(UNIT_L_PLOT).d), 100.0*np.max(GAMER_DensProf_r.in_units(UNIT_L_PLOT).d)], [soliton_fitting_rho_c.in_units(UNIT_D_PLOT).d,            soliton_fitting_rho_c.in_units(UNIT_D_PLOT).d             ],  ':',  color='xkcd:pink',          label=r'$\rho_{\rm c}$',    zorder=1  )
      ax.plot( [soliton_fitting_r_c.in_units(UNIT_L_PLOT).d,           soliton_fitting_r_c.in_units(UNIT_L_PLOT).d           ], [0.01*np.min(GAMER_DensProf_dens.in_units(UNIT_D_PLOT).d), 100.0*np.max(GAMER_DensProf_dens.in_units(UNIT_D_PLOT).d) ],  '-.', color='xkcd:light purple',  label=r'$r_{\rm c}$',       zorder=1  )

      # annotate the information
      ax.annotate( r'$r_{\rm 1/2}=r_{\rm c}$ = %.1e kpc'%(CoreRadius[(Step==step)][0].in_units(UNIT_L_PLOT).d)+'\n'+
                   r'$M_{\rm 1/2}$      = %.1e M$_{\odot}$'%(CoreMass[(Step==step)][0].in_units(UNIT_M_PLOT).d)+'\n'+
                   r'$r_{\rm 1/2}M_{\rm 1/2}$ = %.1e M$_{\odot}$kpc'%(CoreRadiusCoreMass[(Step==step)][0].in_units('kpc*Msun').d)+'\n'+
                   r'$\rho_{\rm peak}$     = %.1e M$_{\odot}$/kpc$^3$'%(CoreDensity[(Step==step)][0].in_units(UNIT_D_PLOT).d)+'\n'+
                   r'$\rho_{\rm c}$         = %.1e M$_{\odot}$/kpc$^3$'%(soliton_fitting_rho_c.in_units(UNIT_D_PLOT).d),
                   xy=(0.02,0.02), xycoords='axes fraction' )

      # set the limits and scales
      if step == 0:
         xlim_min =  0.5*np.min(GAMER_DensProf_r.in_units(UNIT_L_PLOT).d)
         xlim_max = 2e2
         ylim_min = 1e2
         ylim_max = 50.0*np.max(GAMER_DensProf_dens.in_units(UNIT_D_PLOT).d)

      ax.set_xlim( xlim_min, xlim_max )
      ax.set_ylim( ylim_min, ylim_max )
      ax.set_xscale('log')
      ax.set_yscale('log')

      # set the labels
      ax.legend(loc='upper right')
      ax.set_xlabel( r'$r$'+r' (${\rm kpc}$)'    )
      ax.set_ylabel( r'$\rho$'+r' ($M_{\odot}/{\rm kpc}^3$)' )
      #ax.set_title( '$t$ = %7.6e Gyr'%ds.current_time.in_units('Gyr') )
      ax.set_title( 'Step = %09d, Time = %14.4e Myr'%(Step[(Step==step)][0], Time.in_units('Myr').d[(Step==step)][0]) )

      # set the grid and ticks
      ax.grid()
      ax.xaxis.set_ticks_position('both')
      ax.yaxis.set_ticks_position('both')
      ax.tick_params( which='both',direction='in' )

      # save the figure
      plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
      fig.savefig( 'fig_DensityProfile_%09d.png'%(step), dpi=150 )
      fig.clear()
      plt.close()
