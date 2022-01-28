import time
from dataclasses import dataclass
from typing import List, Iterable, Optional

import numpy as np
from numpy import log, exp, log10


@dataclass
class TrialGrid:
    """Stores trial grid parameters"""

    size: int
    period_min: float
    period_max: float


@dataclass
class Result:
    """ Stores fit results.

    Attributes
    ----------
    period : float
        Best-fit period, in seconds
    period_uncertainty : float
        Uncertainty on best-fit period, in seconds
    toas : list
        List of TOAs that were provided as input
    toa_uncertainties : list
        List of TOA uncertainties that were provided as input
    grid : TrialGrid
        Parameters of the period trial grid that was applied
    rotation_indices : list
        List of fitted integer rotation indices, one per input TOA. The rotation index
        associated with the first provided TOA is 0.
    formatted_period : str
        Best-fit period formatted following the usual pulsar astronomy convention, e.g. 1.05(1)
    scaled_toa_uncertainties : list
        The input TOA uncertainties scaled by a constant factor so that they
        would yield a reduced chi-square of 1.
    residuals : list
        Fit residuals (in seconds) on each time interval between each TOA and the first one.
        The fit residual on the first TOA is always 0, because the first TOA is used as the
        time reference.
    solve_time : float
        Total time spent to get a solution, in seconds.
    """

    period: float
    period_uncertainty: float
    toas: List[float]
    toa_uncertainties: List[float]
    grid: TrialGrid
    rotation_indices: List[int]
    formatted_period: str
    scaled_toa_uncertainties: List[float]
    residuals: List[float]
    solve_time: float


DAY_SECONDS = 86400.0


def dot3(a, b, c) -> float:
    """ Returns the matrix product A.B.C, implicitly assuming that the result is a scalar """
    return np.dot(a, np.dot(b, c))[0, 0]


def format_uncertain_quantity(quantity: float, uncertainty: float) -> str:
    """ Format a number with an associated uncertainty following the usual pulsar astronomy convention.

    >>> format_uncertain_quantity(1.05, 0.01)
    1.05(1)
    """
    decimals = -int(np.floor(log10(uncertainty)))
    qr = round(quantity, decimals)
    ur = int(round(uncertainty * 10 ** decimals))
    return f"{qr:.{decimals}f}({ur})"


def rratsolve(
    toas: Iterable[float], toa_uncertainty: float, max_grid_size: Optional[int] = 100_000_000,
) -> Result:
    """ Find the longest spin period that fits a sparse set of single pulse TOAs.

    The method is based on finding a common periodicity that divides all the time intervals
    between the first and all subsequent TOAs. Only period is fitted; initial phase is ignored.
    Step 1: Find the optimal spacing between consecutive trial periods
    Step 2: Try periods between:
        * Period min: 10 times the TOA uncertainty and
        * Period max: the smallest interval between all TOA pairs plus 10%,
        and pick the period that yields the smallest RMS phase residuals
    Step 3: Refine the solution assuming that the rotation counts associated to each TOA
        have been inferred correctly. In this case there is an analytical solution, which
        the one returned.

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
    result : Result
        Result object (dataclass) which wraps all the outputs as attributes

    See Also
    --------
    Result : storage class for the outputs of this function
    """
    start_time = time.time()

    n = len(toas)
    if not n >= 3:
        raise ValueError("Need at least 3 TOAs")

    toa_uncertainties = np.repeat(toa_uncertainty, n)
    iref = 0
    tref = toas[iref]
    toas = np.asarray(toas)
    T = (toas - tref) * DAY_SECONDS

    T = np.delete(T, iref).reshape(-1, 1)
    sigma = np.delete(toa_uncertainties, iref)

    C = np.diag(sigma ** 2) + sigma[iref] ** 2
    M = np.linalg.inv(C)

    pmin = 10 * sigma.max()
    pmax = 1.1 * abs(T).min() + 10 * sigma.max()
    delta_logp = n ** 0.5 * dot3(T.T, M, T) ** -0.5

    pgrid = exp(np.arange(log(pmin), log(pmax), delta_logp))

    if max_grid_size is not None and pgrid.size > max_grid_size:
        raise ValueError(
            "Trial grid size would exceed allowed maximum. "
            "Try reducing the time span of the TOAs "
            "or increasing the estimated TOA uncertainty"
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

    # Best-fit trial grid point
    iopt = q.argmin()
    Kopt = K[:, iopt].reshape(-1, 1)

    # Assuming the turn counts are correct, we can now refine the best-fit period
    pstar = dot3(T.T, M, Kopt) / dot3(Kopt.T, M, Kopt)
    Dstar = T / pstar
    Kstar = Dstar.round().astype(int)
    Rstar = Dstar - Kstar
    Qstar = (np.dot(M, Rstar) * Rstar).sum() * pstar ** 2
    time_residuals = (Rstar * pstar).ravel()

    # Uncertainty scaling factor such that Qstar = n - 1
    uscale = (Qstar / (n - 1)) ** 0.5

    # 1-sigma uncertainty on Pstar
    # NOTE: On some artificial inputs, uncertainty can be exactly zero, and we make sure that
    # does not happen
    pstar_uncertainty = max(pstar * dot3(T.T, M, T) ** -0.5 * uscale, np.finfo(float).eps)
    formatted_period = format_uncertain_quantity(pstar, pstar_uncertainty)

    # NOTE: must cast to int from np.int64 to avoid JSON serialization problems later
    # Also, the rotation index of the first TOA is always 0
    rotation_indices = [0] + list(map(int, Kopt.ravel()))
    time_residuals = [0] + list(time_residuals)

    end_time = time.time()

    result = Result(
        period=pstar,
        period_uncertainty=pstar_uncertainty,
        toas=list(toas),
        toa_uncertainties=list(toa_uncertainties),
        grid=trial_grid,
        rotation_indices=rotation_indices,
        formatted_period=formatted_period,
        scaled_toa_uncertainties=list(toa_uncertainties * uscale),
        residuals=time_residuals,
        solve_time=end_time - start_time,
    )
    return result
