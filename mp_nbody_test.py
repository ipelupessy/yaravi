from mpmath import mp
import numpy
import copy
import cPickle

from matplotlib import pyplot

import itertools

import mp_integrator
from mp_integrator import bulirschStoer,particle,total_energy#,error_controlled_BS
from mp_integrator import Local_Processor,pp_Processor,AmuseProcessor,MultiProcessor

def plummer(N):
    from amuse.ic.plummer import new_plummer_model
 
    numpy.random.seed(123456)
    
    plum=new_plummer_model(N)
    
    particles=[]
    for p in plum:
        particles.append(particle(p.mass.number, p.x.number,p.y.number,p.z.number,
                                      p.vx.number,p.vy.number,p.vz.number))
    return particles                                   
    
def BS_test(N=16,processor="local"):
    import time

    nproc=4
    nslices=nproc

    if processor=="multi":
      mp_integrator.pproc=MultiProcessor(nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=["emphyrio"]*nproc,nslices=nslices,preamble="from mpmath import mp",pre_pickle=True)
    elif processor=="pp":
      mp_integrator.pproc=pp_Processor(nslices=nslices,pre_pickle=True)
    else:
      mp_integrator.pproc=Local_Processor()

    mp.dps=64

    parts=plummer(N) 
        
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
    
    """
    x=numpy.array(x)
    y=numpy.array(y)
    z=numpy.array(z)

    c=['r','g','b','y','m','c','p']
    for i in range(len(parts)):
      pyplot.plot(x[:,i],y[:,i],c[i%7])
    pyplot.xlim(-3,3)
    pyplot.ylim(-3,3)
    pyplot.savefig('test.png')
    """
    
#    print x[-1],y[-1]

    print parts[0].x
    return total_energy(parts)

def check_BS_test(N=16):
  hashes={ 16:-1805142856, 32: -1690459273}
  for p in ["multi","amuse","pp","local"]:
    h=hash(BS_test(N=N,processor=p))
    if hashes[N]==h:
      print p+" ok"
    else:
      print p+" nook:", h


def ec_BS_test():
    import time


    mp_integrator.pproc=MultiProcessor(nslices=4,pre_pickle=True)
#    mp_integrator.pproc=AmuseProcessor(hosts=["emphyrio"]*4,nslices=4,preamble="from mpmath import mp",pre_pickle=True)
#    mp_integrator.pproc=pp_Processor(nslices=4,pre_pickle=True)
#    mp_integrator.pproc=Local_Processor()

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

def time_kick(N=16,processor="local"):
    import time
    from mp_integrator import kick

    nproc=4
    nslices=nproc

    if processor=="multi":
      mp_integrator.pproc=MultiProcessor(nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=["emphyrio"]*nproc,nslices=nslices,
       preamble="from mpmath import mp",pre_pickle=True,
       channel_type="mpi")
    elif processor=="pp":
      mp_integrator.pproc=pp_Processor(nslices=nslices,pre_pickle=True)
    else:
      mp_integrator.pproc=Local_Processor()

    mp.dps=64
    mp_integrator.pproc.exec_("mp.dps="+str(mp.dps))
    
    parts=plummer(N)

    dt=mp.mpf(1)/mp.mpf(8)
    
    t1=time.time()
    kick(parts,parts,dt)
    t2=time.time()

    print "time:",t2-t1
    return hash(parts[-1].vz)
    
def check_kick(N=16):
    hashes={ 16: 1762445124, 50:3277754040,150:-1946492655,64: -3541374424,256:-2250164681,
       512: 2466512669}
    for p in ["multi","amuse","pp","local"]:
      h=hash(time_kick(N=N,processor=p))
      if hashes[N]==h:
        print p+" ok"
      else:
        print p+" nook:", h

if __name__=="__main__":
#    import cProfile
#    cProfile.run('BS_test()','prof')
#    check_BS_test(N=16)
#    from mp_integrator_test import BS_test
     BS_test(N=128,processor="amuse")
#     check_kick(N=512)
#     time_kick(N=256,processor="amuse")
