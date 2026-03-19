import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.cosmology import FlatLambdaCDM, z_at_value
from scipy.interpolate import interp1d

# =============================================================================
# COSMOLOGY
# =============================================================================

cosmo = FlatLambdaCDM(H0=67.8, Om0=0.307, Ob0=0.0483)


# =============================================================================
# COSMOLOGICAL UTILITIES
# =============================================================================

def z_to_a(z):
    """Redshift -> scale factor"""
    return 1.0 / (1.0 + z)


def a_to_z(a):
    """Scale factor -> redshift"""
    return (1.0 - a) / a


def z_to_age(z):
    """Redshift -> age of the Universe (Gyr)"""
    return cosmo.age(z).to(u.Gyr).value


def age_to_z(age_Gyr):
    """Age of the Universe (Gyr) -> redshift."""
    if np.isscalar(age_Gyr):
        return z_at_value(cosmo.age, age_Gyr * u.Gyr).value
    else:
        return np.array([z_at_value(cosmo.age, a * u.Gyr).value
                         for a in np.asarray(age_Gyr)])


def universe_grid_generator(z_init=127, z_final=0, n_snapshots=271):
    """
    Reproduce the GALFORM snapshot grid (uniform in scale factor).

    Returns
    -------
    age_grid : ndarray, shape (n_snapshots,)
        Age of the Universe at each snapshot (Gyr), ascending.
    """
    a_init  = z_to_a(z_init)
    a_final = z_to_a(z_final)
    a_grid  = np.linspace(a_init, a_final, n_snapshots, endpoint=True)
    z_grid  = a_to_z(a_grid)
    return z_to_age(z_grid)


# =============================================================================
# SFH UTILITIES
# =============================================================================

def sfr_interpolator(sfh_time, sfh_sfr):
    """
    Build a linear interpolator for the SFR(t) curve.

    Parameters
    ----------
    sfh_time : array_like  Cosmic time (Gyr), not required to be sorted.
    sfh_sfr  : array_like  SFR (M_sun / yr) at each time.

    Returns
    -------
    f : callable  f(t) -> SFR; returns 0 outside the tabulated range.
    """
    sfh_time = np.asarray(sfh_time, dtype=float)
    sfh_sfr  = np.asarray(sfh_sfr,  dtype=float)
    order    = np.argsort(sfh_time)
    return interp1d(sfh_time[order], sfh_sfr[order],
                    kind='linear', bounds_error=False, fill_value=0.0)


# =============================================================================
# CSP SED CALCULATOR  (the core routine)
# =============================================================================

def build_csp_sed(sfh_time, sfh_sfr, t_obs,
                  ssp_ages, ssp_seds,
                  min_ssp_age=0.0):
    """
    Compute a composite stellar population (CSP) SED by convolving
    a star formation history with a grid of SSP templates.

    The standard CSP integral is:

        F_CSP(lambda, t_obs) =
            integral_{0}^{t_obs}  SFR(t') * S(lambda, t_obs - t') dt'

    where t' is cosmic time and (t_obs - t') is the lookback age of
    the stars formed at t'.  The integral is evaluated by trapezoidal
    quadrature on the native SFH time grid, extended to t_obs.

    Units
    -----
    * sfh_time  : Gyr  (cosmic time since Big Bang)
    * sfh_sfr   : M_sun / yr
    * ssp_ages  : Gyr  (SSP age = how long ago the stars formed)
    * ssp_seds  : L_sun / Angstrom  per 1 M_sun of stars formed
    * output    : L_sun / Angstrom  (total galaxy luminosity)

    The SFR is in M_sun/yr but time steps are in Gyr, so the trapezoidal
    integration naturally produces M_sun when we apply the 1e9 factor
    (1 Gyr = 1e9 yr).

    Parameters
    ----------
    sfh_time  : array_like, shape (N,)
        Cosmic times at which the SFR is known (Gyr).  Need not be sorted.
    sfh_sfr   : array_like, shape (N,)
        Star formation rate at each time (M_sun / yr).
    t_obs     : float
        Cosmic time of observation / galaxy snapshot (Gyr).
        Stars formed after t_obs are ignored.
    ssp_ages  : array_like, shape (M,)
        Age grid of the SSP templates (Gyr), must be sorted ascending.
    ssp_seds  : array_like, shape (M, n_wave)
        SSP SED templates in L_sun / Angstrom per solar mass, one per age.
    min_ssp_age : float, optional
        Minimum SSP age to use (Gyr).  Bins whose lookback age falls below
        this are assigned the youngest template.  Default 0 (no clipping).

    Returns
    -------
    csp_sed : ndarray, shape (n_wave,)
        CSP SED in L_sun / Angstrom.
    t_grid  : ndarray
        Cosmic-time grid actually used for integration (Gyr).
    sfr_grid : ndarray
        Interpolated SFR on t_grid (M_sun / yr).
    mass_grid : ndarray
        Stellar mass formed in each trapezoidal bin (M_sun).
    """
    sfh_time = np.asarray(sfh_time, dtype=float)
    sfh_sfr  = np.asarray(sfh_sfr,  dtype=float)
    ssp_ages = np.asarray(ssp_ages, dtype=float)
    ssp_seds = np.asarray(ssp_seds, dtype=float)

    # ------------------------------------------------------------------
    # 1. Sort and clip the SFH to [0, t_obs]
    # ------------------------------------------------------------------
    order    = np.argsort(sfh_time)
    t_sorted = sfh_time[order]
    s_sorted = sfh_sfr[order]

    # Keep only times strictly before t_obs, then append t_obs as the
    # final point using linear interpolation of the SFR.
    mask    = t_sorted < t_obs
    t_clip  = t_sorted[mask]
    s_clip  = s_sorted[mask]

    sfr_interp = sfr_interpolator(sfh_time, sfh_sfr)
    t_grid  = np.append(t_clip,  t_obs)
    sfr_grid = np.append(s_clip, sfr_interp(t_obs))

    # ------------------------------------------------------------------
    # 2. Trapezoidal mass in each bin  [M_sun]
    #
    #    delta_M_i = 0.5 * (SFR_{i+1} + SFR_i) * delta_t_i
    #
    #    SFR is in M_sun/yr, delta_t is in Gyr
    #    => multiply by 1e9 to convert Gyr -> yr
    # ------------------------------------------------------------------
    delta_t    = np.diff(t_grid)            # Gyr
    sfr_avg    = 0.5 * (sfr_grid[1:] + sfr_grid[:-1])  # M_sun / yr
    mass_bins  = sfr_avg * delta_t * 1e9    # M_sun

    # ------------------------------------------------------------------
    # 3. SSP age for each bin  [Gyr]
    #
    #    The representative cosmic time of bin i is the bin MIDPOINT:
    #        t_mid_i = 0.5 * (t_{i+1} + t_i)
    #    The SSP age (how old those stars are at t_obs) is:
    #        tau_i   = t_obs - t_mid_i
    #
    #    Using the midpoint is the most accurate single-point quadrature
    #    rule (it is exact for linear SFR variation within the bin).
    #    Front/back are first-order and disagree when bins are wide.
    # ------------------------------------------------------------------
    t_mid     = 0.5 * (t_grid[1:] + t_grid[:-1])  # cosmic time of bin midpoint
    tau_bins  = t_obs - t_mid                       # SSP age = lookback time (Gyr)

    # Clip ages to the SSP template range to avoid extrapolation artefacts
    tau_bins  = np.clip(tau_bins, ssp_ages[0], ssp_ages[-1])

    # ------------------------------------------------------------------
    # 4. Accumulate CSP SED
    #
    #    F_CSP(lambda) = sum_i  delta_M_i  *  SSP(lambda, tau_i)
    # ------------------------------------------------------------------
    n_wave   = ssp_seds.shape[1]
    csp_sed  = np.zeros(n_wave, dtype=float)

    for i, tau in enumerate(tau_bins):
        if mass_bins[i] <= 0.0:
            continue
        ssp_i    = _interpolate_ssp(tau, ssp_ages, ssp_seds)
        csp_sed += mass_bins[i] * ssp_i

    return csp_sed, t_grid, sfr_grid, mass_bins


# =============================================================================
# SSP INTERPOLATION  (internal helper)
# =============================================================================

def _interpolate_ssp(tau, ssp_ages, ssp_seds):
    """
    Linearly interpolate the SSP SED grid at age tau.

    Parameters
    ----------
    tau      : float         Target SSP age (Gyr).
    ssp_ages : ndarray (M,)  Sorted SSP age grid (Gyr).
    ssp_seds : ndarray (M, n_wave)

    Returns
    -------
    ssp_interp : ndarray (n_wave,)
    """
    if tau <= ssp_ages[0]:
        return ssp_seds[0].copy()
    if tau >= ssp_ages[-1]:
        return ssp_seds[-1].copy()

    # Binary search for bracketing indices
    j = np.searchsorted(ssp_ages, tau)   # ssp_ages[j-1] < tau <= ssp_ages[j]
    i = j - 1

    w = (tau - ssp_ages[i]) / (ssp_ages[j] - ssp_ages[i])  # 0 <= w <= 1
    return (1.0 - w) * ssp_seds[i] + w * ssp_seds[j]


# =============================================================================
# SANITY CHECK HELPER
# =============================================================================

def check_total_mass(mass_bins, sfh_time, sfh_sfr, t_obs):
    """
    Cross-check the trapezoidal mass sum against a direct numpy trapz
    integration of the raw SFH.  Prints a comparison.
    """
    sfh_time = np.asarray(sfh_time, dtype=float)
    sfh_sfr  = np.asarray(sfh_sfr,  dtype=float)
    order    = np.argsort(sfh_time)
    t_s = sfh_time[order]
    s_s = sfh_sfr[order]
    mask = t_s <= t_obs
    M_direct = np.trapz(s_s[mask], t_s[mask]) * 1e9   # M_sun

    M_bins   = np.sum(mass_bins)
    print(f"  Mass from bin sum  : {M_bins:.4e} M_sun")
    print(f"  Mass from np.trapz : {M_direct:.4e} M_sun")
    print(f"  Ratio              : {M_bins/M_direct:.6f}  (should be ~1)")


# =============================================================================
# MAIN  –  example usage
# =============================================================================

if __name__ == "__main__":

    # --- Load data -----------------------------------------------------------
    wave = np.genfromtxt('./../data/FSPS_wave.csv')

    ssp_seds = np.load('./../data/FSPS_SED_templates.npy', mmap_mode='r')
    ssp_ages = np.genfromtxt('./../data/FSPS_param_grid_LHS.csv')
    # If ssp_ages is 2-D (e.g. Latin Hypercube parameter table), extract age column:
    # ssp_ages = ssp_ages[:, 0]
    print(f"SSP age grid  : shape={ssp_ages.shape}, "
          f"min={ssp_ages.min():.4f} Gyr, max={ssp_ages.max():.4f} Gyr")
    print(f"SSP SED array : shape={ssp_seds.shape}")
    print(f"Wavelength    : shape={wave.shape}")

    # --- Load SFH ------------------------------------------------------------
    SFH_input = np.genfromtxt('./../data/GALFORM_SFHs/gal_1_SFH_test.csv')
    order     = np.argsort(SFH_input[0])
    SFH_time  = SFH_input[0][order]           # Gyr (cosmic time)
    SFH_sfr   = 10 ** SFH_input[1][order]     # M_sun / yr  (log -> linear)

    # --- Observation epoch ---------------------------------------------------
    Galaxy_lookback_time     = 5.0                              # Gyr
    t_universe_today         = z_to_age(0)                     # ~13.8 Gyr
    t_obs                    = t_universe_today - Galaxy_lookback_time

    print(f"\nUniverse age today : {t_universe_today:.4f} Gyr")
    print(f"t_obs (galaxy)     : {t_obs:.4f} Gyr")

    # --- Compute CSP SED -----------------------------------------------------
    csp_sed, t_grid, sfr_grid, mass_bins = build_csp_sed(
        SFH_time, SFH_sfr, t_obs,
        ssp_ages, ssp_seds
    )

    # --- Sanity check: total stellar mass ------------------------------------
    print("\n--- Mass sanity check ---")
    check_total_mass(mass_bins, SFH_time, SFH_sfr, t_obs)

    # --- Plot ----------------------------------------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))

    # 1. SFH
    ax = axes[0]
    ax.plot(t_grid, sfr_grid, 'k-o', ms=3)
    ax.set_xlabel("Cosmic time (Gyr)")
    ax.set_ylabel(r"SFR ($M_\odot$ yr$^{-1}$)")
    ax.set_title("Star Formation History")
    ax.set_yscale('log')

    # 2. Mass per bin
    ax = axes[1]
    t_mid = 0.5 * (t_grid[1:] + t_grid[:-1])
    ax.bar(t_mid, mass_bins, width=np.diff(t_grid), align='center',
           color='steelblue', edgecolor='k', linewidth=0.4)
    ax.set_xlabel("Cosmic time (Gyr)")
    ax.set_ylabel(r"$\Delta M$ ($M_\odot$)")
    ax.set_title("Stellar Mass Formed per Bin")

    # 3. CSP SED
    ax = axes[2]
    ax.plot(wave, csp_sed, 'b-', lw=1.2, label='CSP SED (midpoint)')
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.set_ylabel(r"Luminosity ($L_\odot\,\AA^{-1}$)")
    ax.set_title("Composite Stellar Population SED")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()

    plt.tight_layout()
    plt.savefig('./../data/csp_sed_output.png', dpi=150)
    plt.show()
    print("\nDone. Figure saved.")
