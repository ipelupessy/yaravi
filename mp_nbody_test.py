from mpmath import mp
import numpy
import copy
import cPickle

from matplotlib import pyplot

import itertools

import mp_nbody
from mp_nbody import bulirschStoer,particle,total_energy,error_controlled_BS
from mp_nbody import Local_Processor,pp_Processor,AmuseProcessor,MultiProcessor

def plummer(N):
    from amuse.ic.plummer import new_plummer_model
 
    numpy.random.seed(123456)
    
    plum=new_plummer_model(N)
    
    particles=[]
    for p in plum:
        particles.append(particle(p.mass.number, p.x.number,p.y.number,p.z.number,
                                      p.vx.number,p.vy.number,p.vz.number))
    return particles                                   
    
def BS_test():
    import time

    N=10
    nproc=4
    nbunch=N/nproc

#    mp_nbody.pproc=MultiProcessor(nbunch=nbunch,pre_pickle=True,nproc=nproc)
    mp_nbody.pproc=AmuseProcessor(hosts=["emphyrio"]*nproc,nbunch=nbunch,preamble="from mpmath import mp",pre_pickle=True)
#    mp_nbody.pproc=pp_Processor(nbunch=nbunch,pre_pickle=True)
#    mp_nbody.pproc=Local_Processor()

    mp.dps=64

#    parts=pythagorean()
#    parts=two(mratio=1.e9)
    parts=plummer(N) 
    
#    print parts[0].__module__
#    raise
    
    t=0
    tend=mp.mpf(1.)/8
    dt=mp.mpf(1)/mp.mpf(8)
    dt_param=mp.mpf('1.')

    integrator=bulirschStoer(mp.mpf('1.e-16'),dt_param)

    e0=total_energy(parts)

    x=[[part.x.__float__() for part in parts]]
    y=[[part.y.__float__() for part in parts]]
    z=[[part.z.__float__() for part in parts]]

    t1=time.time()
    while t<tend:
      t+=dt
      integrator.evolve(parts,dt)
      print "t,de:",t,((total_energy(parts)-e0)/e0).__float__()
      x.append([part.x.__float__() for part in parts])
      y.append([part.y.__float__() for part in parts])
      z.append([part.z.__float__() for part in parts])
    t2=time.time()
  
    e=total_energy(parts)
    print 'de/e0:', ((e-e0)/e0).__float__()
  
    print 'runtime:', t2-t1
    print integrator.nsteps,integrator.jcount,integrator.jcount/float(integrator.nsteps)
    print integrator.nkickstep    
    
    x=numpy.array(x)
    y=numpy.array(y)
    z=numpy.array(z)

    c=['r','g','b','y','m','c','p']
    for i in range(len(parts)):
      pyplot.plot(x[:,i],y[:,i],c[i%7])
    pyplot.xlim(-3,3)
    pyplot.ylim(-3,3)
    pyplot.savefig('test.png')
    
#    print x[-1],y[-1]

    print parts[0].x

def ec_BS_test():
    import time


    mp_nbody.pproc=MultiProcessor(nbunch=16,pre_pickle=True)
#    mp_nbody.pproc=AmuseProcessor(hosts=["emphyrio"]*4,nbunch=16,preamble="from mpmath import mp",pre_pickle=True)
#    mp_nbody.pproc=pp_Processor(nbunch=4,pre_pickle=True)
#    mp_nbody.pproc=Local_Processor()

    mp.dps=64

#    parts=pythagorean()
#    parts=two(mratio=1.e9)
    parts=plummer(10) 
    
#    print parts[0].__module__
#    raise
    
    t=0
    tend=mp.mpf(.125)
    dt=mp.mpf(1)/mp.mpf(8)

    integrator=error_controlled_BS(mp.mpf('1.e-16'))
    integrator.particles=parts
    integrator.commit_particles()
    
    e0=total_energy(parts)

    x=[[part.x.__float__() for part in parts]]
    y=[[part.y.__float__() for part in parts]]
    z=[[part.z.__float__() for part in parts]]

    t1=time.time()
    while t<tend:
      t+=dt
      integrator.evolve(t)
      x.append([part.x.__float__() for part in integrator.particles])
      y.append([part.y.__float__() for part in integrator.particles])
      z.append([part.z.__float__() for part in integrator.particles])
      t2=time.time()
      print "t,de:",t,((total_energy(integrator.particles)-e0)/e0).__float__(),t2-t1
  
    e=total_energy(integrator.particles)
    print 'de/e0:', ((e-e0)/e0).__float__()
  
    print 'runtime:', t2-t1
#    print integrator.nsteps,integrator.jcount,integrator.jcount/float(integrator.nsteps)
#    print integrator.nkickstep    
    
    x=numpy.array(x)
    y=numpy.array(y)
    z=numpy.array(z)

    c=['r','g','b','y','m','c','p']
    for i in range(len(integrator.particles)):
      pyplot.plot(x[:,i],y[:,i],c[i%7])
    pyplot.xlim(-3,3)
    pyplot.ylim(-3,3)
    pyplot.savefig('test.png')
    
#    print x[-1],y[-1]

    print parts[0].x

    
if __name__=="__main__":
#    import cProfile
#    cProfile.run('BS_test()','prof')
#    from mp_nbody_test import BS_test
     ec_BS_test()

