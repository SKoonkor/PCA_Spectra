import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from astropy.cosmology import Planck13, z_at_value
from scipy.interpolate import interp1d
import astropy.units as u
import fsps
import os, sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from composite_stellarpop import build_composite_SFH, SFH_tau, SFH_delayed_tau, SFH_lognormal, SFH_DPL
from composite_stellarpop import SFH_burst, build_composite_SFH
from FSPS_calculations import FSPS_initializer
from visualisations import plot_SED

#############################################################################
# Step 0: Essential functions

# Interpolate table
z_grid = np.linspace(0, 200, 2000)
t_universe = Planck13.age(z_grid).to(u.Gyr).value
t_since_bb = Planck13.age(0).to(u.Gyr).value - t_universe

# Interpolators
t_to_z_interp = interp1d(t_since_bb, z_grid, bounds_error = False, fill_value = (z_grid[-1], 0))
z_to_t_interp = interp1d(z_grid, t_since_bb, bounds_error = False, fill_value = (t_since_bb[0], t_since_bb[-1]))

def t_to_z(t_gyr):
    return t_to_z_interp(t_gyr)
   
def z_to_t(z):
    return z_to_t_interp(z)


#############################################################################
# Step 1. Setup age grid and SFH
print ('\nStep 1: Setup age grid and SFH')

N_grid = 200
t_grid = np.linspace(0, 13.8, N_grid)
t_grid_size = t_grid[1] - t_grid[0]


# Example composite SFH (you can arbitorily adjusted)



#############################################################################
# Step 2. Define SFH component sets for each galaxy
print ('\nStep 2: Define SFH components for each galaxy')

SFH_gal1 = [
        {'type': 'delayed_tau', 'T0': 1.0, 'tau' : 1.5, 'strength': 0.5},
        {'type': 'burst', 'tburst': 1, 'strength' : 0.1, 'sigma': 0.03},
        {'type': 'burst', 'tburst': 2, 'strength' : 0.3, 'sigma': 0.2},
        {'type': 'burst', 'tburst': 5, 'strength' : 0.2, 'sigma': 0.5},
        {'type': 'burst', 'tburst': 8, 'strength' : 0.5, 'sigma': 1.},]
SFH_gal2 = [{'type': 'delayed_tau', 'T0': 1, 'tau': 1.5, 'strength': 1}]
SFH_gal3 = [{'type': 'delayed_tau', 'T0': 0, 'tau': 6, 'strength': 0.01},
            {'type': 'burst', 'tburst': 10, 'strength': 0.8, 'sigma': 0.3}]

SFH_sets = [SFH_gal1, SFH_gal2, SFH_gal3]
print ('The SFH components:')
for i, sfh_i in enumerate(SFH_sets):
    print ('\n                 :::Galaxy {}:::'.format(i+1))
    for comp in sfh_i:
        print (comp)

galaxies = [
        {'name': 'Galaxy A', 'components': SFH_gal1},
        {'name': 'Galaxy B', 'components': SFH_gal2},
        {'name': 'Galaxy B', 'components': SFH_gal3},]

for gal in galaxies:
    gal['SFH'] = build_composite_SFH(t_grid, gal['components'])
    gal['csp'] = FSPS_initializer(sfh = 3, dust_emission = False, neb_emission = -1)
    gal['track'] = {'g-r': [], 'Mr': []}
    gal['mass_formed'] = gal['SFH'] * np.gradient(t_grid)
    gal['cum_mass_formed'] = [np.sum(gal['mass_formed'][:i]) for i in range(200)]
    
#############################################################################
# Step 3. Setup figure
print ('\nStep 3: Setup figure')

# The x-axis ticks and scaling
z_ticks = np.array([0, 0.1, 0.3, 0.6, 1, 1.5, 3, 10, 20])
t_ticks = 13.8 - z_to_t(z_ticks)
print (t_ticks)
print (z_ticks)


fig = plt.figure(figsize = (12, 12))

gs = gridspec.GridSpec(3, 2, width_ratios = [1, 2], height_ratios = [1, 1, 1], figure = fig)


ax_sfh = fig.add_subplot(gs[0, :])
ax_spec = fig.add_subplot(gs[1, :])
ax_cc = fig.add_subplot(gs[2, 0])
ax_mass = fig.add_subplot(gs[2, 1])


colors = ['teal', 'crimson', 'indigo', 'goldenrod', 'olivedrab'] # Extend if needed (more SFHs)
max_sfh_plot = 1.1*np.max([np.max(gal['SFH']) for gal in galaxies])


def animate(frame):
    print (f'\rFrame : {frame}/200', end = '', flush = True)
    ax_sfh.clear()
    ax_spec.clear()
    ax_cc.clear()

    for j, gal in enumerate(galaxies):
        # --- Udpate FSPS with current SFH ---

        gal['csp'].set_tabular_sfh(age = t_grid, sfr = gal['SFH'])
        tage = t_grid[frame]
        wave, spec = gal['csp'].get_spectrum(tage = tage, peraa = True)

        # --- Magnitudes for colour-colour ---
        mags = gal['csp'].get_mags(tage = tage, bands = ['sdss_g', 'sdss_r', 'sdss_i'])
        g_r = mags[0] - mags[1]
        Mr = mags[1]

        gal['track']['g-r'].append(g_r)
        gal['track']['Mr'].append(Mr)

        # --- SFH sumplot ---
        ax_sfh.plot(t_grid, gal['SFH'], color = colors[j], lw = 2)
        ax_sfh.axvline(tage, ls = '--', color = 'k', alpha = 0.5)
        ax_sfh.set_ylabel(r'SFR $M_{sun}/yr]$')
        ax_sfh.set_xlabel('Cosmic Age [Gyr]')
        ax_sfh.set_xlim(0, 14)
        ax_sfh.set_ylim(0, max_sfh_plot)
        ax_sfh.text(7 - 0.45, max_sfh_plot*1.05, 'redshift')
        ax_sfh.text(0.5, 0.9, 'Star Formation Histories', fontsize = 15)

        for tt, t_tick in enumerate(t_ticks):
            ax_sfh.axvline(t_tick, ymin = 0.99, color = 'grey')
    
            if z_ticks[tt] < 10:
                ax_sfh.text(t_tick - 0.16, max_sfh_plot*1.01, f'{z_ticks[tt]:.1f}')
            else:
                ax_sfh.text(t_tick - 0.13, max_sfh_plot*1.01, '{}'.format(int(z_ticks[tt])))

        # --- Spectrum subplot ---
        plot_SED(ax_spec, wave, spec*wave, color = colors[j], lw = 2)
        ax_spec.set_ylabel(r'$\lambda F_\lambda [L_{sun}]$')
        ax_spec.set_xlim(1e2, 1e5)
        ax_spec.set_ylim(1e2, 1e11)
        ax_spec.set_xlabel('wavelength [A]')
        ax_spec.set_xscale('log')
        ax_spec.set_yscale('log')
        ax_spec.text(128, 2.2e9, 'Spectra', fontsize = 15)

        
        # --- Colour-magnitude subplot ---
        ax_cc.scatter(gal['track']['Mr'], gal['track']['g-r'], alpha = 0.6, s = 2,  color = colors[j])
        ax_cc.scatter(Mr, g_r, s = 40, marker = 'o', color = colors[j])
        ax_cc.set_xlabel(r'$M_r$')
        ax_cc.set_ylabel('g - r')
        ax_cc.set_xlim(-13, -19)
        ax_cc.set_ylim(-0.1, 0.9)
        ax_cc.text(-13.8, 0.7, 'CMD', fontsize = 15)

        # --- Mass formed ----
        ax_mass.scatter(t_grid[:frame], gal['cum_mass_formed'][:frame]/np.sum(gal['mass_formed']), color = colors[j], s=1)
        ax_mass.scatter(tage, gal['cum_mass_formed'][frame], color = colors[j], s = 5, marker = 'o')
        ax_mass.set_ylim(0, 1)
        ax_mass.set_xlim(0, 14)
        ax_mass.set_ylabel('Normalised Cumulative Mass Formed')
        ax_mass.set_xlabel('Cosmic Age [Gyr]')
        ax_mass.text(0.6, 0.815, 'Mass formed', fontsize = 15)




##############################################################################
# Step 4. Animation function
print ('\nStep 4: Animation function')
print (f'The age grid size: {t_grid.shape}')

def init():
    for line in spec_lines: 
        line.set_data([], [])
    for marker in sfh_markers:
        marker.set_xdata([0, 0])
    return spec_lines + sfh_markers

################################################################################
# Step 5. Run animation
print ('\nStep 5: Run animation')

ani = animation.FuncAnimation(fig, animate, frames = range(1, len(t_grid)), interval = 100)

# Save as MP4 (requires ffmpeg installed)
ani.save('../outputs/MultiGalaxy_SFH_SED.mp4', writer = 'ffmpeg', dpi = 150)

plt.show()
