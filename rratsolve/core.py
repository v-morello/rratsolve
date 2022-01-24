import logging
import time
from dataclasses import dataclass
from typing import List

import numpy as np
from numpy import log, exp, log10


@dataclass
class TrialGrid:
    """ Stores trial grid parameters """
    size: int
    period_min: float
    period_max: float

@dataclass
class Result:
    """ Stores fit results """
    toas: List[float]
    toa_uncertainties: List[float]
    grid: TrialGrid
    rotation_indices: List[int]
    period: float
    period_uncertainty: float
    formatted_period: str
    scaled_toa_uncertainties: List[float]
    residuals: List[float]
    solve_time: float


DAY_SECONDS = 86400.0


logger = logging.getLogger('rratsolve')


def dot3(a, b, c):
    return np.dot(a, np.dot(b, c))[0, 0]


def format_uncertain_quantity(q, u):
    decimals = -int(np.floor(log10(u)))
    qr = round(q, decimals)
    ur = int(round(u * 10**decimals))
    return f"{qr:.{decimals}f}({ur})"


def rratsolve(toas, toa_uncertainty, max_grid_size=None):
    """
    Find the longest spin period that fits a sparse set of single pulse TOAs.

    Parameters
    ----------
    toas : list or ndarray
        Pulse arrival MJDs.
    toa_uncertainty : float
        The estimated TOA RMS uncertainty in second.
    max_grid_size : int or None
        Maximum allowed number of points in the trial grid. If None, no limit is enforced.
        If not None and the limit is exceeded, raise ValueError.

    Returns
    -------
    p : float
        The best-fit period.
    delta_p : float
        The estimated 1-sigma period uncertainty.
    """
    start_time = time.time()

    n = len(toas)
    toa_uncertainties = np.repeat(toa_uncertainty, n)
    iref = 0
    #logger.debug(f"Using TOA #{iref} as reference")

    tref = toas[iref]
    toas = np.asarray(toas)
    T = (toas - tref) * DAY_SECONDS

    T = np.delete(T, iref).reshape(-1, 1)
    sigma = np.delete(toa_uncertainties, iref)
    #logger.debug(f"Time intervals: {T.ravel()}")

    C = np.diag(sigma**2) + sigma[iref]**2
    M = np.linalg.inv(C)

    pmin = 10 * sigma.max()
    pmax = 1.1 * abs(T).min() + 10 * sigma.max()
    delta_logp = n**0.5 * dot3(T.T, M, T) ** -0.5

    pgrid = exp( np.arange(log(pmin), log(pmax), delta_logp) )
    #logger.debug(f"Min trial period: {pgrid[0]:.6f}")
    #logger.debug(f"Max trial period: {pgrid[-1]:.6f}")
    #logger.debug(f"Period grid size: {pgrid.size:,}")

    if max_grid_size is not None and pgrid.size > max_grid_size:
        raise ValueError(
            "Trial grid size would exceed allowed maximum. "
            "Try reducing the time span of the TOAs or increasing the estimated TOA uncertainty"
        )

    trial_grid = TrialGrid(size=pgrid.size, period_min=pgrid[0], period_max=pgrid[-1])

    # Estimated fractional turn counts
    D = T / pgrid

    # Estimated integer turn counts
    K = D.round().astype(int)

    # Phase residuals
    R = D - K

    # NOTE: This is q = Q / P^2
    # Q is a chi2 with n-1 degrees of freedom
    q = (np.dot(M, R) * R).sum(axis=0)

    iopt = q.argmin()
    popt = pgrid[iopt]
    Kopt = K[:, iopt].reshape(-1, 1)
    #logger.debug(f"Initial period guess: {popt:.6f}")

    pstar = dot3(T.T, M, Kopt) / dot3(Kopt.T, M, Kopt)
    #logger.debug(f"Refined period guess: {pstar:.6f}")

    Dstar = T / pstar
    Kstar = Dstar.round().astype(int)
    Rstar = Dstar - Kstar
    Qstar = (np.dot(M, Rstar) * Rstar).sum() * pstar**2

    time_residuals = (Rstar * pstar).ravel()
    #logger.debug(f"Time intervals residuals: {time_residuals}")

    # Uncertainty scaling factor such that Qstar = n - 1
    uscale = (Qstar / (n - 1)) ** 0.5
    #logger.debug(f"Uncertainty scaling factor: {uscale:.4f}")
    #logger.debug(f"Scaled TOA uncertainties: {toa_uncertainties * uscale}")

    # 1-sigma uncertainty on Pstar
    pstar_uncertainty = pstar * dot3(T.T, M, T)**-0.5 * uscale

    sol_str = format_uncertain_quantity(pstar, pstar_uncertainty)
    #logger.debug(f"Best-fit period: {sol_str} s")

    # NOTE: must cast to int from np.int64 to avoid JSON serialization problems later
    # Also, the rotation index of the first TOA is always 0
    rotation_indices = [0] + list(map(int, Kopt.ravel()))
    time_residuals = [0] + list(time_residuals)

    end_time = time.time()
    #logger.debug(f"Run time: {end_time - start_time:.3f} s")

    result = Result(
        toas=list(toas),
        toa_uncertainties=list(toa_uncertainties),
        grid=trial_grid,
        rotation_indices=rotation_indices,
        period=pstar,
        period_uncertainty=pstar_uncertainty,
        formatted_period=sol_str,
        scaled_toa_uncertainties=list(toa_uncertainties * uscale),
        residuals=time_residuals,
        solve_time=end_time - start_time
    )
    return result
