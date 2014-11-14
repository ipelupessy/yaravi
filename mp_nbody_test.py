import mp
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
    
def BS_test(N=16,processor="local",tend=1./8,prec='1.e-16',res="energy",dps=64):
    import time

    nproc=4
    nslices=nproc

    if processor=="multi":
      mp_integrator.pproc=MultiProcessor(preamble="import mp",nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=[None]*nproc,nslices=nslices,preamble="import mp",pre_pickle=True)
    elif processor=="pp":
      preamble="import mp;from mp_integrator import _kick,_potential,_timestep,particle,jparticle"
      mp_integrator.pproc=pp_Processor(preamble=preamble,nslices=nslices,pre_pickle=True)
    elif processor=="proc":
      mp_integrator.pproc=Processor(preamble="import mp")
    else:
      mp_integrator.pproc=Local_Processor()

    mp.set_dps(dps)
    mp_integrator.pproc.exec_("mp.set_dps("+str(mp.dps)+")")
    
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
  
    print N,processor+' runtime:', t2-t1
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

#    print parts[0].x,total_energy(parts)
    if res=="energy":
      return total_energy(parts)
    else:
      return str(parts[0].x)[0:32]

def check_BS_test(Ns=None):
  dps=64
  results={16: "-0.2097788238899595006650143725272014434855335456003604589938999139641"}
  if Ns is None: Ns=results.keys()
  check=True
  for N in Ns:
    for p in ["multi","amuse","pp","local","proc"]:
        h=BS_test(N=N,processor=p,dps=dps+6)
        hr=h.__repr__()
        if hr.find(results[N])>=0:
          pass
#          print p+" ok"
        else:
          d=abs(h-mp.mpf(results[N]))/h
          print d
          print N,p+" mismatch, got:", hr,len(hr),float(abs(mp.log10(d)))
          if abs(d)>10**-dps: check=False
  return check
      
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
#    mp_integrator.pproc=AmuseProcessor(hosts=[None]*4,nslices=4,preamble="import mp",pre_pickle=True)
#    mp_integrator.pproc=pp_Processor(nslices=4,pre_pickle=True)
#    mp_integrator.pproc=Local_Processor()

    mp.set_dps(64)
    mp_integrator.pproc.exec_("mp.set_dps("+str(mp.dps)+")")

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

def time_kick(N=16,processor="local",dps=64):
    import time
    from mp_integrator import kick

    nproc=4
    nslices=nproc

    if processor=="multi":
      mp_integrator.pproc=MultiProcessor(preamble="import mp",nslices=nslices,pre_pickle=True,nproc=nproc)
    elif processor=="amuse":
      mp_integrator.pproc=AmuseProcessor(hosts=[None]*nproc,nslices=nslices,preamble="import mp",pre_pickle=True)
    elif processor=="pp":
      preamble="import mp;from mp_integrator import _kick,_potential,_timestep,particle,jparticle"
      mp_integrator.pproc=pp_Processor(preamble=preamble,nslices=nslices,pre_pickle=True)
    elif processor=="proc":
      mp_integrator.pproc=Processor(preamble="import mp")
    else:
      mp_integrator.pproc=Local_Processor()

    mp.set_dps(dps)
    mp_integrator.pproc.exec_("mp.set_dps("+str(mp.dps)+")")
    
    parts=plummer(N)

    dt=mp.mpf(1)/mp.mpf(8)
    
    t1=time.time()
    kick(parts,parts,dt)
    t2=time.time()

    print N,processor+" time:",t2-t1
    return parts[-1].vz
    
def check_kick(Ns=None):
    dps=64
    results={ 16: '0.5280968268833149585365524469736931751421937461245719231626600116459',
              50: '0.9570351506440379530178266489555425661480746175752162677540570048262',
              64: '-0.7561747845694916564830825081707286817639452371097957709612294389958',
              150: "-0.09344415268620983521705314736452944414180898050922914284511780006364",
              256: "-0.6293985589493140191821290296902117972778559732348390605690993188598",
              350: "0.735469970603869347834820970508731867404741985456426980852100139499",
              512: "0.9517024891436366308665374005000831926596595154320406804913855436092"}
    check=True
    if Ns is None: Ns=results.keys()
    for N in Ns:
      for p in ["multi","amuse","pp","local","proc"]:
        h=time_kick(N=N,processor=p,dps=dps)
        hr=h.__repr__()
        if hr.find(results[N])>=0:
          pass
#          print p+" ok"
        else:
          d=abs(h-mp.mpf(results[N]))/h
          print N,p+" mismatch, got:", hr,len(hr),float(abs(mp.log10(d)))
          if d>10**-dps: check=False
    return check
      
if __name__=="__main__":
     assert check_BS_test()
#     assert check_kick()
#    BS_test(N=5,processor="pp",tend=10.,prec='1.e-6',res="energy")
#    import cProfile
#    cProfile.run('BS_test()','prof')
#    from mp_integrator_test import BS_test
#     BS_test(N=8,processor="local",tend=2.)
#     time_kick(N=256,processor="amuse")
#    check_BS_test(N=32)
#    BS_test(N=256,processor="amuse",tend=1.,prec='1.e-6',res="energy")
