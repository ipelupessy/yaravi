import numpy
import time

import hashlib

from amuse.units import nbody_system

from amuse.community.yaravi.interface import Yaravi

from amuse.ic.plummer import new_plummer_model

def test_run(N=4,tend=1. | nbody_system.time,dt=0.125 | nbody_system.time, parameters={}):

  numpy.random.seed(1233231)

  mp=Yaravi()
    
  for key in parameters:
    setattr(mp.parameters,key,parameters[key])
    
  print mp.parameters

  plum=new_plummer_model(N,do_scale=True)

  mp.particles.add_particles(plum)

  t1=time.time()
  while mp.model_time < tend-dt/2:
    tnext=mp.model_time+dt
    if mp.model_time+dt>=tend-dt/2:
      tnext=tend 
    mp.evolve_model(tnext)
  t2=time.time()

  print t2-t1,(tend-mp.model_time).number

  return mp.particles.copy()

class TestYaravi(object):

  def test1(self):
    """ N=4, test different dt """
    p1=test_run(N=4,tend=1. | nbody_system.time,dt=(1./3) | nbody_system.time)
    p2=test_run(N=4,tend=1. | nbody_system.time,dt=1. | nbody_system.time)
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 1 sucessful"

  def test2(self):
    """ N=8, test different dt """
    p1=test_run(N=8,tend=1. | nbody_system.time,dt=(1./3) | nbody_system.time)
    p2=test_run(N=8,tend=1. | nbody_system.time,dt=1. | nbody_system.time)
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 2 sucessful"

  def test3(self):
    """ test factor """
    p1=test_run(parameters=dict(factor=1.e5))
    p2=test_run(parameters=dict(factor=1.e10))
    p3=test_run(parameters=dict(factor=1.e15))
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
    sha3=hashlib.sha1()
    sha3.update(p3.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest() == sha3.hexdigest()
    print "test 3 sucessful"

  def test4(self):
    """ test timestep_paramter """
    p1=test_run(parameters=dict(timestep_parameter=1.))
    p2=test_run(parameters=dict(timestep_parameter=.1))
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 4 sucessful"

  def test5(self):
    """ test timestep_paramter """
    p1=test_run(parameters=dict(initial_dps=20))
    p2=test_run(parameters=dict(initial_dps=40))
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 5 sucessful"
    
  def test6(self):
    """ test timestep_paramter """
    p1=test_run(parameters=dict(initial_target_error=1.e-16))
    p2=test_run(parameters=dict(initial_target_error=1.e-40))
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 6 sucessful"    

  def test7(self):
    """ N=16, test different dt, long """
    p1=test_run(N=16,tend=3. | nbody_system.time,dt=(1./7) | nbody_system.time)
    p2=test_run(N=16,tend=3. | nbody_system.time,dt=1. | nbody_system.time)
  
    sha1=hashlib.sha1()
    sha1.update(p1.position.number.data)
    sha2=hashlib.sha1()
    sha2.update(p2.position.number.data)
  
    assert sha1.hexdigest() == sha2.hexdigest()
    print "test 7 sucessful"
  
if __name__=="__main__":
  test=TestYaravi()
  test.test1()
  test.test2()
  test.test3()
  test.test4()
  test.test5()
  test.test6()
  test.test7()
