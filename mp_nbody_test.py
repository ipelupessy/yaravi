from mpmath import mp
import numpy
import copy
import cPickle

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

import itertools

import mp_integrator
from mp_integrator import bulirschStoer,particle,total_energy#,error_controlled_BS
from mp_integrator import Local_Processor,pp_Processor,AmuseProcessor,MultiProcessor,Processor

def plummer(N):
    from amuse.ic.plummer import new_plummer_model
 
    numpy.random.seed(123456)
    
    plum=new_plummer_model(N)
    
    particles=[]
    for p in plum:
        particles.append(particle(p.mass.number, p.x.number,p.y.number,p.z.number,
                                      p.vx.number,p.vy.number,p.vz.number))
    return particles                                   
    
def BS_test(N=16,processor="local",tend=1./8,prec='1.e-16',res="energy"):
    import time

    nproc=4
    nslices=nproc

    if processor=="multi":
      mp_integrator.pproc=MultiProcessor(preamble="from mpmath import mp",nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=[None]*nproc,nslices=nslices,preamble="from mpmath import mp",pre_pickle=True)
    elif processor=="pp":
      preamble="from mpmath import mp;from mp_integrator import _kick,_potential,_timestep,particle,jparticle"
      mp_integrator.pproc=pp_Processor(preamble=preamble,nslices=nslices,pre_pickle=True)
    elif processor=="proc":
      mp_integrator.pproc=Processor(preamble="from mpmath import mp")
    else:
      mp_integrator.pproc=Local_Processor()

    mp.dps=64

    parts=plummer(N) 
        
    t=0
    tend=mp.mpf(tend)
    dt=mp.mpf(1)/mp.mpf(8)
    dt_param=mp.mpf('1.')

    integrator=bulirschStoer(mp.mpf(prec),dt_param)

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
    if res=="energy":
      return total_energy(parts)
    else:
      return str(parts[0].x)[0:32]

def check_BS_test(N=16):
  hashes={ 16:-1805142856, 32: -1690459273}
  for p in ["multi","amuse","pp","local","proc"]:
    h=hash(BS_test(N=N,processor=p))
    if hashes[N]==h:
      print p+" ok"
    else:
      print p+" nook:", h

def check_BS_test_long(N=16):
  hashes={ 8: 7787870482253759781}
  for p in ["amuse","local"]:
    h=hash(BS_test(N=N,processor=p,tend=0.5,prec='1.e-32',res="x"))
    if hashes[N]==h:
      print p+" ok"
    else:
      print p+" nook:", h

def ec_BS_test():
    import time


    mp_integrator.pproc=MultiProcessor(nslices=4,pre_pickle=True)
#    mp_integrator.pproc=AmuseProcessor(hosts=[None]*4,nslices=4,preamble="from mpmath import mp",pre_pickle=True)
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
      mp_integrator.pproc=MultiProcessor(preamble="from mpmath import mp",nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=[None]*nproc,nslices=nslices,preamble="from mpmath import mp",pre_pickle=True)
    elif processor=="pp":
      preamble="from mpmath import mp;from mp_integrator import _kick,_potential,_timestep,particle,jparticle"
      mp_integrator.pproc=pp_Processor(preamble=preamble,nslices=nslices,pre_pickle=True)
    elif processor=="proc":
      mp_integrator.pproc=Processor(preamble="from mpmath import mp")
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
       512: 2466512669,128:-2819487303}
    for p in ["multi","amuse","pp","local","proc"]:
      h=hash(time_kick(N=N,processor=p))
      if hashes[N]==h:
        print p+" ok"
      else:
        print p+" nook:", h

if __name__=="__main__":
#    import cProfile
#    cProfile.run('BS_test()','prof')
#    from mp_integrator_test import BS_test
#     BS_test(N=8,processor="local",tend=2.)
     check_kick(N=128)
#     time_kick(N=256,processor="amuse")
#    check_BS_test(N=16)
