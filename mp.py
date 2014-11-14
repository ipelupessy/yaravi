import math

BACKEND="mpmath"

#matches mpmath
def prec_to_dps(n):
    """Return number of accurate decimals that can be represented
    with a precision of n bits."""
    return max(1, int(round(int(n)/3.3219280948873626)-1))

def dps_to_prec(n):
    """Return the number of bits required to represent n decimals
    accurately."""
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

def repr_dps(n):
    """Return the number of decimal digits required to represent
    a number with n-bit precision so that it can be uniquely
    reconstructed from the representation."""
    dps = prec_to_dps(n)
    if dps == 15:
        return 17
    return dps + 3

if BACKEND=="mpmath":
    from mpmath import *
        
    def set_dps(new_dps):
      global dps
      mp.dps=new_dps
      dps=mp.dps
      
    def set_precision(new_precision):
      mp.prec=new_precision

    dps=mp.dps
  
elif BACKEND=="gmpy2":
    from gmpy2 import *

    mpf=mpfr
        
    def set_dps(new_dps):
      global dps
      get_context().precision = dps_to_prec(new_dps)    
      dps=new_dps

    def set_precision(new_precision):
      get_context().precision = new_precision
      dps=prec_to_dps(get_context().precision)

    dps=prec_to_dps(get_context().precision)

  
