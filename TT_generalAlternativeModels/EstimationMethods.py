import numpy as np
import statsmodels.api as sm


def detrend_quadratic(timeseries):
    fit = np.polyfit(range(len(timeseries)), timeseries, 2)
    timeseries = timeseries - [x**2 * fit[0] + x * fit[1] + fit[2] for x in range(len(timeseries))]
    return timeseries


def phi_gls(timeseries):
    timeseries = np.array(timeseries)
    timeseries = detrend_quadratic(timeseries)
    diff = timeseries[1:] - timeseries[:-1]
    model = sm.GLSAR(diff, timeseries[:-1], rho=1)
    result = model.iterative_fit(maxiter=100)
    return result.params[0]+1
    