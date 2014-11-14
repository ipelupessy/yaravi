import math

BACKEND="gmpy2"

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
      get_context().precision = int(log2(10.) * int(new_dps))    
      dps=new_dps

    def set_precision(new_precision):
      get_context().precision = new_precision
      dps=int(new_precision/log2(10.))

    dps=int(get_context().precision/log2(10.))

  
