# rratsolve

A brute-force but robust pulsar periodicity solver that needs only a few (but at least 3) sparse TOAs to infer a period. Has proven to be quite useful for getting an initial period guess for newly discovered radio sources in single pulse searches.

### Web App

A handy web version of RRATSolve is [currently hosted here](https://rratsolve.herokuapp.com/).

### Installation

The recommended method is to clone the repository and install in development mode, preferably in a virtualenv/conda environment:
```
pip install -e .
```

### Usage

The module exposes a single function called `rratsolve()` which is fully documented. Example:

```python
>>> from rratsolve import rratsolve

# Some actual topocentric TOAs of PSR J1819-1458
>>> toas = [
    59607.3866952756,
    59607.3930600505,
    59607.4015472327,
    59607.4185223751,
    59607.4286862792,
    59607.4470411329,
    59607.4591301745,
    59607.4689014586
]

>>> result = rratsolve(toas, toa_uncertainty=10.0e-3)
>>> print(result.period, result.period_uncertainty)
4.263244257764385 3.036291589650858e-05

>>> print(result.formatted_period)
4.26324(3)
```

This is consistent with the published _barycentric_ period of `P = 4.26321348` seconds (see [Bhattacharyya et al., 2018](https://arxiv.org/abs/1803.10277))

### Method

The solver attempts to find the **largest** period that is an integer divisor of all time intervals between the first and all subsequent input TOAs. The "initial phase term" is therefore ignored. It can be shown that the trial grid should be geometrically spaced; the multiplicative factor between consecutive trial periods is chosen so that at least one grid point falls within the expected 1-sigma confidence interval for the true period. This "expected" confidence interval depends explicitly on the `toa_uncertainty` input parameter (the smaller the uncertainty, the more tightly packed the grid).

After testing the whole trial grid, the solver picks the best-fit points and assumes it is close enough to the true period so that the inferred (integer) rotation counts are correct. From there there is an analytical solution, that the code returns.
