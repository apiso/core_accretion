"""
array creation routines beyond numpy's linspace and logspace:

Functions
---------
mylogspace : log spaced samples by specifying array values (not their log).
linlogspace : a mix of linear and log spacing with continuous stepsizes.
loglinspace : equivalent to linlogspace, but with arguments in different order.
"""

from scipy.optimize import brentq
import numpy as np
import sys

def mylogspace(start, stop, num = 50, endpoint = True):
    """
    Return log spaced numbers over a specified interval.  In this alternate to numpy's logspace
    `start` and `stop` are given as numbers in the array not their logaritm.  Also handles
    negative values of `start` and `stop`.
    
    The endpoint of the interval can optionally be excluded.

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The end value of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        log spaced samples, so that `stop` is excluded.  Note that the step
        size changes when `endpoint` is False.
    num : int, optional
        Number of samples to generate. Default is 50.
    endpoint : bool, optional
        If True, `stop` is the last sample. Otherwise, it is not included.
        Default is True.

    Returns
    -------
    samples : ndarray
        There are `num` equally spaced samples in the closed interval
        ``[start, stop]`` or the half-open interval ``[start, stop)``
        (depending on whether `endpoint` is True or False).
    """

    if start > 0 and stop > 0:
        samples = np.logspace(np.log10(start), np.log10(stop), num, endpoint = endpoint)
    elif start < 0 and stop < 0:
        samples = np.logspace(np.log10(-start), np.log10(-stop), num, endpoint = endpoint)
        samples = -samples
    else:
        sys.exit("start and stop must both be positive or both be negative")

    samples[0]  = start #correct for numerical errors in log and exponentiation
    if endpoint:
        samples[-1] = stop
        
    return samples
    
def linlogspace(linend, logend, num = 50, nlin = 10):
    """
    Return numbers with a mix of linear and log spacing.  The first log step
    (after switching from linear spacing) has the same size as the linear steps.

    Returns `num` samples, calculated over the interval [`linend`, `logend` ].
    There are `nlin` linearly spaced samples and the remainder are evenly spaced.
    Choosing `linend` > `logend` is allowed and the array will be returned in
    ascending order (starting from smaller of `logend` or  `linend`).

    Parameters
    ----------
    linend : scalar
        The terminal value at the linear end.  If too negative, no solution is possible. 
    logend : scalar
        The terminal value at the log end, must be > 0.
    num : int, optional
        Number of samples to generate. Default is 50.
    nlin : int, optional
        Number of linear samples to generate. Default is 10.
        
    Returns
    -------
    samples : ndarray
        There are `num` samples in the closed interval ``[logend, linend]`` or
        ``[linend, logend]`` (depending on which gives ascending order).
    """
    return loglinspace(logend,linend,num,num-nlin)

def loglinspace(logend, linend, num = 50, nlog = 10):
    """
    Return numbers with a mix of log and linear spacing.  The final log step
    (before switching to linear spacing) has the same size as the linear steps.

    Returns `num` samples, calculated over the interval [`logend`, `linend` ].
    There are `nlog` log spaced samples and the remainder are evenly spaced.
    Choosing `logend` > `linend` is allowed and the array will be returned in
    ascending order (starting from smaller of `logend` or  `linend`).

    Parameters
    ----------
    logend : scalar
        The terminal value at the log end, must be > 0.
    linend : scalar
        The terminal value at the linear end.  If too negative, no solution is possible. 
    num : int, optional
        Number of samples to generate. Default is 50.
    nlog : int, optional
        Number of log samples to generate. Default is 10.
        
    Returns
    -------
    samples : ndarray
        There are `num` samples in the closed interval ``[logend, linend]`` or
        ``[linend, logend]`` (depending on which gives ascending order).
    """
    nlin = num - nlog
    if nlin < 1:
        sys.exit("at least one linear point")
    
    if logend <= 0:
        sys.exit("positive values at logend required.")
    if type(nlog) != int or type(num) != int:
        sys.exit("integer numbers of array elements!")
    if nlog < 2:
        sys.exit("at least 2 log points")
    
    samples = np.empty(num)
    
    M , N = nlog, nlin #to match notation in defining functions
    if linend == 0:
        switch = logend * (nlin/(nlin+1.))**(nlog - 1)
    elif linend > 0:
        switch = rootswitch(logend, linend, M, N) #
    else: #linend < 0
        testmin = _rootfnmin(logend, M, N)
        if _rootfn(testmin, logend, linend, M,N) > 0:
            sys.exit("No solution possible, choose more positive linend or lower nlog")
        switch = rootswitch(logend, linend, M, N, linsearch = testmin)
    samples[0] = logend
    for i in range(1,nlog-1):
        samples[i] = samples[i-1] * (switch/logend)**(1./(nlog - 1))
    samples[nlog-1 : ] = np.linspace(switch,linend,nlin+1)
    if logend < linend:
        return samples
    else:
        return samples[::-1]

def _rootfn(switch,logend,linend,M,N):
    """
    Used to find the numerical value to switch from lin to log spacing as follows.
    Consider two arrays:
        -a length M log-spaced array from [logend,switch] (with both logend and switch > 0)
        -a length N+1 linearly spaced array from [switch,linend]
    This rootfn returns zero when the stepsize from switch to the next value in both arrays is equal.
    Note that the concatenation of both arrays has length M+N.

    To derive this equation let:
        -`logend` = x[0] (first element of log array)
        -`switch` = x[M-1] = y[0] (last of log and first of lin)
        -`linend` = y[N] (last element of lin)
    equating dy = (y[N] - y[0])/N to dx[M-2] = x[M-1] - x[M-2] = y[0]*(1 - (x[0]/y[0])**(1/(M-1.)))
    gives the desired result
    """
    return switch  - N/(N+1.) * logend**(1./(M-1)) * switch**((M-2.)/(M-1)) - linend/(N + 1.)

def _rootfnmin(logend, M, N):
    """
    The switch value that minimizes rootfn.
    This value is needed to determine if solution exists when linend < 0.
    A solution only exists if rootfn is negative at this switch value."""
    return logend * ( N/(N + 1.) * (M-2.)/(M-1.) )**(M-1)
      
def _rootfnarray(logend,linend,M,N,nsample = 100):
    """for testing the behavior of the root fn (or positive values)"""
    if linend < 0:
        linendsample = 1e-9 #arbitrary small number
    else:
        linendsample = linend
    switcharr = np.linspace(logend,linendsample,nsample)
    return switcharr, _rootfn(switcharr, logend, linend, M, N)
        
def rootswitch(logend,linend,M,N,linsearch = None):
    if linsearch == None:
        linsearch = linend
    swargs = (logend, linend, M, N)
    switch = brentq(_rootfn, logend, linsearch, args = swargs)
    return switch

