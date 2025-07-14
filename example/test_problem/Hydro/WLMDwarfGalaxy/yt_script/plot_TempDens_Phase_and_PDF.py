import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid


###################################################################################################
plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
colormap    = 'inferno'
dpi         = 150
###################################################################################################

def log_bin_width(centers):
    log_centers = np.log10( centers )

    log_edges = np.concatenate( ([log_centers[0]  - 0.5*(log_centers[1]   - log_centers[0] )],
                                                    0.5*(log_centers[:-1] + log_centers[1:]),
                                 [log_centers[-1] + 0.5*(log_centers[-1]  - log_centers[-2])]) )

    return np.diff( log_edges )


def plot_PDF( fig, pos_in_fig, ds, obj, prefix, field_type, field_x, weighting_type, lim_x_min, lim_x_max, nbin, units_x ):

    assert weighting_type == 'volume' or weighting_type == 'mass', 'Incorrect weighting type'

    prof = yt.create_profile( obj, (field_type, field_x), fields=(field_type, weighting_type),
                              weight_field=None,
                              n_bins=nbin,
                              logs={(field_type, field_x): True, (field_type, weighting_type): True},
                              units={(field_type, field_x): units_x},
                              extrema={(field_type, field_x): (lim_x_min, lim_x_max)},
                              fractional=False,
                            )

    arr_x = prof.x.in_units(units_x).d
    arr_y = prof[(field_type, weighting_type)].d

    log_bin_w = log_bin_width(arr_x)
    pdf_y = arr_y / (sum(arr_y) * log_bin_w)

    np.savetxt( '%s%s_%s_weighted_%s_PDF'%(ds, prefix, weighting_type, field_x.capitalize()),
                np.column_stack( (arr_x, pdf_y) ),
                fmt='%24.8e',
                header='%22s %24s'%( field_x.capitalize()+'_in_'+units_x, 'PDF') )

    if pos_in_fig == 'x':
       rect = [0.08,0.74,0.588,0.21]
    elif pos_in_fig == 'y':
       rect = [0.76,0.1,0.21,0.6]

    ax = fig.add_axes(rect=rect)

    W = 'V' if weighting_type == 'volume' else 'M'

    if pos_in_fig == 'x':
        ax.loglog( arr_x, pdf_y )

        ax.set_title( 'PDF weighted by %s'%weighting_type )
        #ax.set_xlabel( field_x.capitalize()+' (%s)'%units_x.replace("**", "^") )
        ax.set_ylabel( r'$\frac{1}{%s}\,\frac{\mathrm{d}\,%s}{\mathrm{d}\,\log_{10}{(\mathrm{\frac{%s}{%s}})}}$'%(W, W, field_x.capitalize(), units_x.replace("**", "^")), fontsize='x-large')

        ax.set_xlim( lim_x_min, lim_x_max )
        ax.set_ylim( 1e-8, 2e1 )

    if pos_in_fig == 'y':
        ax.loglog( pdf_y, arr_x )

        #ax.set_ylabel( field_x.capitalize()+' (%s)'%units_x.replace("**", "^"), rotation=90, labelpad=8.0 )
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_xlabel( r'$\frac{1}{%s}\,\frac{\mathrm{d}\,%s}{\mathrm{d}\,\log_{10}{(\mathrm{\frac{%s}{%s}})}}$'%(W, W, field_x.capitalize(), units_x.replace("**", "^")), rotation=0, fontsize='x-large')
        ax.invert_xaxis()

        plt.setp(ax.get_xticklabels(), rotation=0, va="top",    ha="center")
        plt.setp(ax.get_yticklabels(), rotation=0, va="center", ha="left"  )

        ax.set_ylim( lim_x_min, lim_x_max )
        ax.set_xlim( 1e-8, 2e1 )

    # log scale
    ax.set_xscale('log')
    ax.set_yscale('log')

    minor_ticks = LogLocator( base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10 )
    ax.xaxis.get_ticklocs( minor=True )
    ax.yaxis.get_ticklocs( minor=True )
    ax.minorticks_on()
    ax.xaxis.set_minor_locator( minor_ticks )
    ax.xaxis.set_minor_formatter( NullFormatter() )
    ax.xaxis.set_ticks_position( 'both' )
    ax.yaxis.set_ticks_position( 'both' )
    ax.tick_params( which='both', direction='in' )

    ax.grid()


def plot_PhaseDiagram( ds, obj, prefix, field_type, weighting_type, nbin, x_lim_min, x_lim_max, y_lim_min, y_lim_max ):

#   plot
    fig = plt.figure()

    grid = AxesGrid( fig, rect=[0.08, 0.1, 0.6, 0.6], nrows_ncols=(1, 1),
                     axes_pad=0.05, label_mode="L", share_all=True,
                     cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="0%",
                     aspect=False,
                   )

    temp_dens = yt.PhasePlot( obj, (field_type, 'density'), (field_type, 'T'), (field_type, weighting_type),
                              weight_field=None, x_bins=nbin, y_bins=nbin, fractional=True )
    temp_dens.set_xlim( x_lim_min, x_lim_max )
    temp_dens.set_ylim( y_lim_min, y_lim_max )
    temp_dens.set_zlim( (field_type, weighting_type), 5e-10, 5e-3 )
    temp_dens.set_cmap( (field_type, weighting_type), colormap )
    temp_dens.set_colorbar_label( (field_type, weighting_type), weighting_type.capitalize()+' Fraction' )
    temp_dens.annotate_text( xpos=x_lim_min*10**(0.70*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                             ypos=y_lim_min*10**(0.90*(np.log10(y_lim_max)-np.log10(y_lim_min))),
                             text='$t$ = {:.1f} {:s}'.format( ds.current_time.in_units('Myr').d, 'Myr' ),
                             color='black' )

    plot = temp_dens.plots[field_type, weighting_type]
    plot.figure = fig
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]
    temp_dens.render()

    units_w = 'Msun' if weighting_type == 'mass' else 'kpc**3'
    w_total_obj  = obj.quantities.total_quantity( (field_type, weighting_type ) ).in_units(units_w).d
    grid[0].axes.annotate( 'Total {:s} = {:.2e} {:s}'.format( weighting_type, w_total_obj, units_w.replace("**", "^") ),
                           (x_lim_min*10**(0.70*(np.log10(x_lim_max)-np.log10(x_lim_min))),
                            y_lim_min*10**(0.85*(np.log10(y_lim_max)-np.log10(y_lim_min)))),
                           color='black', fontsize='small' )

    plot_PDF( fig, 'x', ds, obj, prefix, field_type, 'density', weighting_type, x_lim_min, x_lim_max, nbin, 'g/cm**3' )
    plot_PDF( fig, 'y', ds, obj, prefix, field_type, 'T',       weighting_type, y_lim_min, y_lim_max, nbin, 'K'       )
    fig.savefig( 'fig_%s%s_Temp-Dens_%s_weighted_PhaseDiagram.png'%(ds, prefix, weighting_type) )
    plt.close()
