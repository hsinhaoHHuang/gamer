import yt
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import scipy.integrate

ds0 = yt.load( '../Data_000000' )

m_22   = 0.2
M_h    = 1.88156694e+10
Mmin_0 = 4.91991577e+08

def X_IntersectionInLogSpace( X_1, X_2, F_1, F_2, G_1, G_2  ):
    # Linear interpolation in log-space
    x_1 = np.log10( X_1 )
    x_2 = np.log10( X_2 )
    f_1 = np.log10( F_1 )
    f_2 = np.log10( F_2 )
    g_1 = np.log10( G_1 )
    g_2 = np.log10( G_2 )

    P_X = 10**( (x_1*f_2 - x_2*f_1 - x_1*g_2 + x_2*g_1)/(g_1 - g_2 - f_1+ f_2) )
    return P_X

def Soliton_rho_c_from_r_c( r_c ):   # in Msun/kpc^3
    return 1.9e7/( m_22**2 * r_c**4 )

def Soliton_analytical_dens( r, r_c ):  # in Msun/kpc^3
    rho_c = Soliton_rho_c_from_r_c( r_c )
    return rho_c/( 1.0 + 9.1e-2*( r/r_c )**2 )**8

def Soliton_analytical_shell_density( r, r_c ):
    return (4.0*np.pi*r**2)*Soliton_analytical_dens( r, r_c )

def Soliton_M_c_from_r_c_Ana( r, r_c ):
    rho_c = Soliton_rho_c_from_r_c( r_c )
    gamma = r/r_c

    M_c = 7.376*rho_c*(r_c**3)*( np.arctan(0.3017*gamma) + (1.0+0.091*gamma**2)**(-7)*
                                                                                      ( -3.017e-1*( gamma**(1)  ) +
                                                                                         3.849e-1*( gamma**(3)  ) +
                                                                                         6.656e-2*( gamma**(5)  ) +
                                                                                         6.651e-3*( gamma**(7)  ) +
                                                                                         3.903e-4*( gamma**(9)  ) +
                                                                                         1.255e-5*( gamma**(11) ) +
                                                                                         1.713e-7*( gamma**(13) )  ) )
    return M_c

def Soliton_M_c_from_r_c_Num( r, r_c ):
    M_c = scipy.integrate.quad( lambda x: Soliton_analytical_shell_density( x, r_c ), 0, r )[0]
    return M_c

def Soliton_r_c_M_c_from_r_c_Ana( r_c ): # in Msun*kpc
    return r_c * Soliton_M_c_from_r_c_Ana( r_c, r_c )

def Soliton_r_c_M_c_from_r_c_Num( r_c ): # in Msun*kpc
    return r_c * Soliton_M_c_from_r_c_Num( r_c, r_c )

Soliton_M_c_from_CHrelation        = 0.25*(M_h/Mmin_0)**(1.0/3.0)*Mmin_0
Soliton_r_c_from_CHrelation        = 1.6/m_22*(M_h/1.0e9)**(-1.0/3.0)
Soliton_rho_c_from_CHrelation      = Soliton_rho_c_from_r_c( Soliton_r_c_from_CHrelation )
Soliton_r_c_M_c_from_CHrelation    = Soliton_r_c_from_CHrelation*Soliton_M_c_from_CHrelation

# load the reference profiles
Step, Time, CoreDensity, CoreRadius, CoreMass, CoreRadiusCoreMass = np.loadtxt( '../Record__CoreProperties', skiprows=1, unpack=True )

# decide the units for plotting
#UNIT_L_PLOT = 'code_length'
#UNIT_D_PLOT = 'code_density'
UNIT_T_PLOT  = 'Myr'
UNIT_L_PLOT  = 'kpc'
UNIT_M_PLOT  = 'Msun'
UNIT_D_PLOT  = 'Msun/kpc**3'
UNIT_LM_PLOT = 'Msun*kpc'

# assign the units
Time                = ds0.arr( Time,               'code_time'             ).in_units(UNIT_T_PLOT).d
CoreDensity         = ds0.arr( CoreDensity,        'code_density'          ).in_units(UNIT_D_PLOT).d
CoreRadius          = ds0.arr( CoreRadius,         'code_length'           ).in_units(UNIT_L_PLOT).d
CoreMass            = ds0.arr( CoreMass,           'code_mass'             ).in_units(UNIT_M_PLOT).d
CoreRadiusCoreMass  = ds0.arr( CoreRadiusCoreMass, 'code_length*code_mass' ).in_units(UNIT_LM_PLOT).d

# number of points for interpolation
N_interp                       = 100*(Time.size - 1)                   # 100 times the original data
dt_interp                      = (Time[-1]-Time[0])/(N_interp)
Time_interp                    = np.arange( Time[0], Time[-1], dt_interp )

# interpolation to make sure it is evenly spaced
CoreDensity_interp             = np.interp( Time_interp, Time, CoreDensity        )
CoreRadius_interp              = np.interp( Time_interp, Time, CoreRadius         )
CoreMass_interp                = np.interp( Time_interp, Time, CoreMass           )
CoreRadiusCoreMass_interp      = np.interp( Time_interp, Time, CoreRadiusCoreMass )

# Gaussian sigma for smoothing
Soliton_tau_00                 = 39.9*(np.min(CoreDensity)/1.0e9)**(-0.5) # Soliton wavefunction oscillation period
Gaussian_sigma_to_tau00_ratio  = 8
Gaussian_sigma                 = Gaussian_sigma_to_tau00_ratio*Soliton_tau_00
Gaussian_truncate              = 5.0           # in sigma

# Gaussian smoothing
CoreDensity_smooth_filt        = gaussian_filter1d( CoreDensity_interp,        sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreRadius_smooth_filt         = gaussian_filter1d( CoreRadius_interp,         sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreMass_smooth_filt           = gaussian_filter1d( CoreMass_interp,           sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )
CoreRadiusCoreMass_smooth_filt = gaussian_filter1d( CoreRadiusCoreMass_interp, sigma=Gaussian_sigma/dt_interp, mode='nearest', truncate=Gaussian_truncate )

# define the threshold value for soliton formation
CoreDensity_threshold_ratio        = 2.0
CoreDensity_threshold_value        = CoreDensity_threshold_ratio*CoreDensity[0]
CoreDensity_threshold_array        = np.full( np.shape(Time_interp), CoreDensity_threshold_value )

CoreRadiusCoreMass_threshold_ratio = 0.9
CoreRadiusCoreMass_threshold_value = CoreRadiusCoreMass_threshold_ratio*Soliton_r_c_M_c_from_r_c_Ana( CoreRadius[0] )
CoreRadiusCoreMass_threshold_array = np.full( np.shape(Time_interp), CoreRadiusCoreMass_threshold_value )

# find the intersection of smooth curve and the threshold value
# peak density intersection
CoreDensity_intersection_idx                = np.argwhere( np.diff(np.sign(CoreDensity_threshold_array - CoreDensity_smooth_filt)) ).flatten()[0]
try:
   CoreDensity_intersection_Time            = X_IntersectionInLogSpace( Time_interp[CoreDensity_intersection_idx],
                                                                        Time_interp[CoreDensity_intersection_idx+1],
                                                                        CoreDensity_threshold_array[CoreDensity_intersection_idx],
                                                                        CoreDensity_threshold_array[CoreDensity_intersection_idx+1],
                                                                        CoreDensity_smooth_filt[CoreDensity_intersection_idx],
                                                                        CoreDensity_smooth_filt[CoreDensity_intersection_idx+1] )
except Exception as e:
   print(e)
   CoreDensity_intersection_Time            = np.nan
   pass

# rcmc intersection
CoreRadiusCoreMass_intersection_idx         = np.argwhere( np.diff(np.sign(CoreRadiusCoreMass_threshold_array - CoreRadiusCoreMass_smooth_filt)) ).flatten()[0]
try:
   CoreRadiusCoreMass_intersection_Time     = X_IntersectionInLogSpace( Time_interp[CoreRadiusCoreMass_intersection_idx],
                                                                        Time_interp[CoreRadiusCoreMass_intersection_idx+1],
                                                                        CoreRadiusCoreMass_threshold_array[CoreRadiusCoreMass_intersection_idx],
                                                                        CoreRadiusCoreMass_threshold_array[CoreRadiusCoreMass_intersection_idx+1],
                                                                        CoreRadiusCoreMass_smooth_filt[CoreRadiusCoreMass_intersection_idx],
                                                                        CoreRadiusCoreMass_smooth_filt[CoreRadiusCoreMass_intersection_idx+1] )
except Exception as e:
   print(e)
   CoreRadiusCoreMass_intersection_Time     = np.nan
   pass

###################################################################################################

# create the figure
fig = plt.figure(figsize=(10,7.5))
ax1  = fig.add_subplot(221)
ax2  = fig.add_subplot(222)
ax3  = fig.add_subplot(223)
ax4  = fig.add_subplot(224)

# plot the profiles
ax1.plot( Time, CoreDensity,         color='C2', linewidth=0.5, label='Core Density' )
ax2.plot( Time, CoreRadius,          color='C3', linewidth=0.5, label='Core Radius'  )
ax3.plot( Time, CoreMass,            color='C0', linewidth=0.5, label='Core Mass'    )
ax4.plot( Time, CoreRadiusCoreMass,  color='C1', linewidth=0.5, label='CoreR*CoreM'  )

# from r_c
ax1.plot( Time, Soliton_rho_c_from_r_c(       CoreRadius             ),  color='xkcd:very light brown', linewidth=0.5, label='_from r_c', zorder=0.5 )
ax3.plot( Time, Soliton_M_c_from_r_c_Ana(     CoreRadius, CoreRadius ),  color='xkcd:very light brown', linewidth=0.5, label='_from r_c', zorder=0.5 )
ax4.plot( Time, Soliton_r_c_M_c_from_r_c_Ana( CoreRadius             ),  color='xkcd:very light brown', linewidth=1.5, label='analytical = %.2e'%Soliton_r_c_M_c_from_r_c_Ana( CoreRadius[0] ), zorder=2.1 )

# smooth curves
ax1.plot( Time_interp, CoreDensity_smooth_filt,        '-', color='xkcd:neon pink',   linewidth=1.5, label="smooth, $\sigma$=%.2e Myr"%Gaussian_sigma, zorder=2.2 )
ax2.plot( Time_interp, CoreRadius_smooth_filt,         '-', color='xkcd:neon green',  linewidth=1.5, label="smooth, $\sigma$=%.2e Myr"%Gaussian_sigma, zorder=2.2 )
ax3.plot( Time_interp, CoreMass_smooth_filt,           '-', color='xkcd:neon yellow', linewidth=1.5, label="smooth, $\sigma$=%.2e Myr"%Gaussian_sigma, zorder=2.2 )
ax4.plot( Time_interp, CoreRadiusCoreMass_smooth_filt, '-', color='xkcd:neon purple', linewidth=1.5, label="smooth, $\sigma$=%.2e Myr"%Gaussian_sigma, zorder=2.2 )

# threshold
ax1.plot( Time_interp, CoreDensity_threshold_array,        '--',  color='xkcd:pink',    linewidth=1.5, label="threshold (%.2f)"%CoreDensity_threshold_ratio,        zorder=2.3  )
ax4.plot( Time_interp, CoreRadiusCoreMass_threshold_array, '--',  color='xkcd:purple',  linewidth=1.5, label="threshold (%.2f)"%CoreRadiusCoreMass_threshold_ratio, zorder=2.3  )

# core-halo relation
ax1.plot( [Time_interp[0], Time_interp[-1]], [Soliton_rho_c_from_CHrelation,   Soliton_rho_c_from_CHrelation  ], ':', color='xkcd:dark indigo', linewidth=1.5, label="Core-Halo Relation (%.2e)"%Soliton_rho_c_from_CHrelation,   zorder=2.5 )
ax2.plot( [Time_interp[0], Time_interp[-1]], [Soliton_r_c_from_CHrelation,     Soliton_r_c_from_CHrelation    ], ':', color='xkcd:dark indigo', linewidth=1.5, label="Core-Halo Relation (%.2e)"%Soliton_r_c_from_CHrelation,     zorder=2.5 )
ax3.plot( [Time_interp[0], Time_interp[-1]], [Soliton_M_c_from_CHrelation,     Soliton_M_c_from_CHrelation    ], ':', color='xkcd:dark indigo', linewidth=1.5, label="Core-Halo Relation (%.2e)"%Soliton_M_c_from_CHrelation,     zorder=2.5 )

# intersection time
if not np.isnan(CoreDensity_intersection_Time):
    ax1.plot( [CoreDensity_intersection_Time, CoreDensity_intersection_Time],               [np.min(CoreDensity), np.max(CoreDensity)],               '-.',  color='xkcd:dark pink',   linewidth=1.5,  label='t = %3.2e Myr'%CoreDensity_intersection_Time,        zorder=2.6 )

if not np.isnan(CoreRadiusCoreMass_intersection_Time):
    ax4.plot( [CoreRadiusCoreMass_intersection_Time, CoreRadiusCoreMass_intersection_Time], [np.min(CoreRadiusCoreMass), np.max(CoreRadiusCoreMass)], '-.',  color='xkcd:dark purple', linewidth=1.5,  label='t = %3.2e Myr'%CoreRadiusCoreMass_intersection_Time, zorder=2.6 )

# set the limits and scales
# ax.set_xlim( xlim_min, xlim_max )
# ax.set_ylim( ylim_min, ylim_max )
# ax.set_xscale('log')
# ax.set_yscale('log')

# set the labels
ax1.legend(loc='upper left')
ax2.legend()
ax3.legend()
ax4.legend()

ax1.set_ylabel( r'$\rho_c$'+r' ($M_{\odot}/{\rm kpc}^3$)' )
ax2.set_ylabel( r'$r_c$'+r' (${\rm kpc}$)' )
ax3.set_xlabel( r'$t$'+r' (${\rm Myr}$)'    )
ax3.set_ylabel( r'$M_c$'+r' ($M_{\odot}$)' )
ax4.set_xlabel( r'$t$'+r' (${\rm Myr}$)'    )
ax4.set_ylabel( r'$r_cM_c$'+r' ($M_{\odot}{\rm kpc}$)' )

# set the grid and ticks
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.tick_params( which='both',direction='in' )
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.tick_params( which='both',direction='in' )
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax3.tick_params( which='both',direction='in' )
ax4.xaxis.set_ticks_position('both')
ax4.yaxis.set_ticks_position('both')
ax4.tick_params( which='both',direction='in' )

# save the figure
plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
fig.savefig( 'fig_CoreProperties.png', dpi=150 )
fig.clear()
