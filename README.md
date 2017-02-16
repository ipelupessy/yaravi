# Yaravi: Yet AnotheR Arithmetic-precision Varying Integrator #
    
Yaravi is a code to solve the astrophysical N-body problem 
in arbitrary precision arithmatic, with a default target of 
the python float( 64-bit floating point) precision. It 
uses an python implementation of a Bulirsch-Stoer integrator with 
adaptive timestepping. It has crude support for parallelization. 

It can be used through the AMUSE interface (included) or seperately if
desired.

Note that (re)commit_particles calls generate a checkpoint in the 
code: the solution found is guaranteed to floating point precision 
from the last check point.

    .. [#] integrator: Boekholt & Portegies Zwart ,Computational Astrophysics and Cosmology, Volume 2, article id.2, 21 pp.
    .. [#] python implementation: F.I. Pelupessy