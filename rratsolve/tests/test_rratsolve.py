from pytest import raises
from rratsolve import rratsolve


def test_check_toa_number():
    """ Check that an error is raised if there are less than 3 TOAs """
    with raises(ValueError):
        rratsolve([], 10.0e-3)

    with raises(ValueError):
        rratsolve([0.0], 10.0e-3)

    with raises(ValueError):
        rratsolve([0.0, 0.01], 10.0e-3)


def test_max_grid_size():
    """ Check that the max grid size parameter is handled properly """
    toas = [0.0, 0.1, 0.2]  # days
    with raises(ValueError):
        rratsolve(toas, 10.0e-3, max_grid_size=1000)
    rratsolve(toas, 10.0e-3, max_grid_size=None)


def test_output():
    """ Test on real topocentric TOAs from PSR J1819-1458 """
    toas = [
        59607.3866952756,
        59607.3930600505,
        59607.4015472327,
        59607.4185223751,
        59607.4286862792,
        59607.4470411329,
        59607.4591301745,
        59607.4689014586
    ]
    bary_period = 4.26321348  # Barycentric period from PSRCAT
    result = rratsolve(toas, 10.0e-3)

    # NOTE: here we cannot check to higher precision because the TOAs above are topocentric
    # instead of barycentric
    assert abs(result.period - bary_period) < 1.0e-4
