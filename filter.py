import math
import operator
import warnings

import numpy
import numpy as np
from numpy import (atleast_1d, poly, polyval, roots, real, asarray,
                   resize, pi, absolute, logspace, r_, sqrt, tan, log10,
                   arctan, arcsinh, sin, exp, cosh, arccosh, ceil, conjugate,
                   zeros, sinh, append, concatenate, prod, ones, full, array,
                   mintypecode)
from numpy.polynomial.polynomial import polyval as npp_polyval
from numpy.polynomial.polynomial import polyvalfromroots

from scipy import special, optimize, fft as sp_fft
from scipy.special import comb
from scipy._lib._util import float_factorial

def _relative_degree(z, p):
    """
    Return relative degree of transfer function from zeros and poles
    """
    degree = len(p) - len(z)
    if degree < 0:
        raise ValueError("Improper transfer function. "
                         "Must have at least as many poles as zeros.")
    else:
        return degree

def lp2lp_zpk(z, p, k, wo=1.0):
    r"""
    Transform a lowpass filter prototype to a different frequency.

    Return an analog low-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency,
    using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : array_like
        Zeros of the analog filter transfer function.
    p : array_like
        Poles of the analog filter transfer function.
    k : float
        System gain of the analog filter transfer function.
    wo : float
        Desired cutoff, as angular frequency (e.g., rad/s).
        Defaults to no change.

    Returns
    -------
    z : ndarray
        Zeros of the transformed low-pass filter transfer function.
    p : ndarray
        Poles of the transformed low-pass filter transfer function.
    k : float
        System gain of the transformed low-pass filter.

    See Also
    --------
    lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
    lp2lp

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s}{\omega_0}

    .. versionadded:: 1.1.0

    """
    # print(p,type(p))
    # z = atleast_1d(z)
    # p = atleast_1d(p)
    # print(p,type(p))
    wo = float(wo)  # Avoid int wraparound

    degree = _relative_degree(z, p)
    
    # Scale all points radially from origin to shift cutoff frequency
    z_lp = wo * z
    p_lp = wo * p

    # Each shifted pole decreases gain by wo, each shifted zero increases it.
    # Cancel out the net change to keep overall gain the same
    k_lp = k * wo**degree

    return z_lp, p_lp, k_lp


def buttap(N):
    """Return (z,p,k) for analog prototype of Nth-order Butterworth filter.

    The filter will have an angular (e.g., rad/s) cutoff frequency of 1.

    See Also
    --------
    butter : Filter design function using this prototype

    """
    # if abs(int(N)) != N:
    #     raise ValueError("Filter order must be a nonnegative integer")
    z = numpy.array([])
    m = numpy.arange(-N+1, N, 2)
    # Middle value is 0 to ensure an exactly real pole
    p = -numpy.exp(1j * pi * m / (2 * N))
    k = 1
    # print("pole",p) pass
    return z, p, k

def bilinear_zpk(z, p, k, fs):
    r"""
    Return a digital IIR filter from an analog one using a bilinear transform.

    Transform a set of poles and zeros from the analog s-plane to the digital
    z-plane using Tustin's method, which substitutes ``(z-1) / (z+1)`` for
    ``s``, maintaining the shape of the frequency response.

    Parameters
    ----------
    z : array_like
        Zeros of the analog filter transfer function.
    p : array_like
        Poles of the analog filter transfer function.
    k : float
        System gain of the analog filter transfer function.
    fs : float
        Sample rate, as ordinary frequency (e.g., hertz). No prewarping is
        done in this function.

    Returns
    -------
    z : ndarray
        Zeros of the transformed digital filter transfer function.
    p : ndarray
        Poles of the transformed digital filter transfer function.
    k : float
        System gain of the transformed digital filter.

    Examples
    --------
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> fs = 100
    >>> bf = 2 * np.pi * np.array([7, 13])
    >>> filts = signal.lti(*signal.butter(4, bf, btype='bandpass', analog=True,
    ...                                   output='zpk'))
    >>> filtz = signal.lti(*signal.bilinear_zpk(filts.zeros, filts.poles,
    ...                                         filts.gain, fs))
    >>> wz, hz = signal.freqz_zpk(filtz.zeros, filtz.poles, filtz.gain)
    >>> ws, hs = signal.freqs_zpk(filts.zeros, filts.poles, filts.gain,
    ...                           worN=fs*wz)
    >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hz).clip(1e-15)),
    ...              label=r'$|H_z(e^{j \omega})|$')
    >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hs).clip(1e-15)),
    ...              label=r'$|H(j \omega)|$')
    >>> plt.legend()
    >>> plt.xlabel('Frequency [Hz]')
    >>> plt.ylabel('Magnitude [dB]')
    >>> plt.grid()
    """
    # z = atleast_1d(z)
    # p = atleast_1d(p)

    degree = _relative_degree(z, p)

    fs2 = 2.0*fs
    # print("fs2",fs2)
    # Bilinear transform the poles and zeros
    z_z = (fs2 + z) / (fs2 - z)
    p_z = (fs2 + p) / (fs2 - p)
    # print("pz",p_z)
    # Any zeros that were at infinity get moved to the Nyquist frequency
    z_z = append(z_z, -ones(degree))

    # Compensate for gain change
    k_z = k * real(prod(fs2 - z) / prod(fs2 - p))

    # print("bi",k_z,prod(fs2 - z),prod(fs2 - p),k)

    return z_z, p_z, k_z

def zpk2tf(z, p, k):
    """
    Return polynomial transfer function representation from zeros and poles

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    """
    # z = atleast_1d(z)
    # k = atleast_1d(k)
    if len(z.shape) > 1:
        temp = poly(z[0])
        b = np.empty((z.shape[0], z.shape[1] + 1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]] * z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * poly(z[i])
    else:
        b = k * poly(z)
    a = atleast_1d(poly(p))
    # print("k",k)
    # print("ploly",poly(z),z)
    # print("ploly",poly(p),p)
    return b, a

def iirfilter(N, Wn, rp=None, rs=None, btype='lowpass', analog=False,
              ftype='butter', output='ba', fs=None):
    """
    Parameters
    ----------
    N : int
        The order of the filter.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.

        For digital filters, `Wn` are in the same units as `fs`. By default,
        `fs` is 2 half-cycles/sample, so these are normalized from 0 to 1,
        where 1 is the Nyquist frequency. (`Wn` is thus in
        half-cycles / sample.)

        For analog filters, `Wn` is an angular frequency (e.g., rad/s).
    rp : float, optional
        For Chebyshev and elliptic filters, provides the maximum ripple
        in the passband. (dB)
    rs : float, optional
        For Chebyshev and elliptic filters, provides the minimum attenuation
        in the stop band. (dB)
    btype : {'bandpass', 'lowpass', 'highpass', 'bandstop'}, optional
        The type of filter.  Default is 'bandpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    ftype : str, optional
        The type of IIR filter to design:

            - Butterworth   : 'butter'
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    output : {'ba', 'zpk', 'sos'}, optional
        Filter form of the output:

            - second-order sections (recommended): 'sos'
            - numerator/denominator (default)    : 'ba'
            - pole-zero                          : 'zpk'

        In general the second-order sections ('sos') form  is
        recommended because inferring the coefficients for the
        numerator/denominator form ('ba') suffers from numerical
        instabilities. For reasons of backward compatibility the default
        form is the numerator/denominator form ('ba'), where the 'b'
        and the 'a' in 'ba' refer to the commonly used names of the
        coefficients used.

        Note: Using the second-order sections form ('sos') is sometimes
        associated with additional computational costs: for
        data-intense use cases it is therefore recommended to also
        investigate the numerator/denominator form ('ba').

    fs : float, optional
        The sampling frequency of the digital system.

        .. versionadded:: 1.2.0

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    """
    Wn = 2*Wn/fs

    z, p, k = buttap(N)

    # Pre-warp frequencies for digital filter design
    fs = 2.0
    warped = 2 * fs * tan(pi * Wn / fs)

    z, p, k = lp2lp_zpk(z, p, k, wo=warped)
    # print("after lp2lp",p,z,k) pass
    z, p, k = bilinear_zpk(z, p, k, fs=fs)
    # print("after bi",p,z,k) 
    # Transform to proper out type (pole-zero, state-space, numer-denom)
    if output == 'zpk':
        return z, p, k
    elif output == 'ba':
        return zpk2tf(z, p, k)
    # elif output == 'sos':
    #     return zpk2sos(z, p, k)


def butter(N, Wn, btype='low', analog=False, output='ba', fs=None):
    """
    N : int
        The order of the filter.
    Wn : array_like
        The critical frequency or frequencies. For lowpass and highpass
        filters, Wn is a scalar; for bandpass and bandstop filters,
        Wn is a length-2 sequence.

        For a Butterworth filter, this is the point at which the gain
        drops to 1/sqrt(2) that of the passband (the "-3 dB point").

        For digital filters, `Wn` are in the same units as `fs`.  By default,
        `fs` is 2 half-cycles/sample, so these are normalized from 0 to 1,
        where 1 is the Nyquist frequency. (`Wn` is thus in
        half-cycles / sample.)

        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba' for backwards
        compatibility, but 'sos' should be used for general-purpose filtering.
    fs : float, optional
        The sampling frequency of the digital system.

        .. versionadded:: 1.2.0

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    buttord, buttap

    Notes
    -----
    The Butterworth filter has maximally flat frequency response in the
    passband.

    The ``'sos'`` output parameter was added in 0.16.0.

    If the transfer function form ``[b, a]`` is requested, numerical
    problems can occur since the conversion between roots and
    the polynomial coefficients is a numerically sensitive operation,
    even for N >= 4. It is recommended to work with the SOS
    representation.
    """
    return iirfilter(N, Wn, btype=btype, analog=analog,
                     output=output, ftype='butter', fs=fs)

def lfilter(b,a,x,sync=False,cycle=0,INsyncdata=[],OUTsyncdata=[]):
    n=len(a)
    x=list(x)
    sync=sync and cycle
    if sync:
        x=list(INsyncdata)+x
        out = np.zeros(len(x))
        out[:n-1]=OUTsyncdata
    else:
        out = np.zeros(len(x))
    
    for i in range(n-1,len(x)):
        su,sv=0,0
        for k in range(n):
            su+=b[k]*x[i-k]
            # print("i-k",i-k)
            # print("temp1 ",i,b[k],x[i-k],su)
        for k in range(1,n):
            sv+=a[k]*out[i-k]
            # print("temp2 ",i,a[k],out[i-k],sv)
        out[i]=(su-sv)
        # print(i,su,sv)
    if sync:
        # print("sync",len(out[n-1:]))
        return out[n-1:]
    # print("unsync",len(out))
    return out


filter_dict = {'butter': [buttap]
               }
# nyq = 0.5 * (50*24)
# normal_cutoff =50 / nyq
if __name__=="__main__":
    b,a =butter(N=8,Wn=50,btype='low',analog=False,fs=50*25)
    fs=50*25
    dt=1/fs
    phase=30
    t=np.arange(0,2,dt)
    sig=np.sin(2*math.pi*50*t+math.radians(phase))+0.06*np.sin(2*math.pi*2*50*t+math.radians(phase/2))
    # for i,j in enumerate(sig):
    #     print(i," ", j)
    filtered=lfilter(b,a,sig)
    filtered_c=[]
    with  open("filtered_c.txt") as ifile:
        for i in ifile.readlines():
            filtered_c.append(float(i)) 
    import matplotlib.pyplot as plt
    # fig,axe=plt.subplot(2)
    # plt.plot(t,filtered)
    plt.plot(filtered_c)
    # plt.plot(t,np.sin(2*math.pi*50*t+math.radians(phase)))
    plt.show()
