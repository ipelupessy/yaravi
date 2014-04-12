from mpmath import mp
import numpy
import copy
import cPickle
from itertools import izip
import processors
import time

max_timestep=1000.

# naieve function to find biggest factor smaller than root
def fsmrt(N):
  n=int(N**0.5)
  while N%n:
    n-=1
  return n

ceil=lambda x,y: (x/y+(x%y>0))

class reducers(object):
  @staticmethod
  def _kick(result,jobresult):
    if result is None:
      result=jobresult
    else:
      result[0]+=jobresult[0]
      result[1]+=jobresult[1]
      result[2]+=jobresult[2]
    return result  
  @staticmethod
  def _potential(result,jobresult):
    if result is None:
      result=jobresult
    else:
      result+=jobresult
    return result  
  @staticmethod
  def _timestep(result,jobresult):
    if result is None:
      result=jobresult
    else:
      if jobresult<result:
        result=jobresult
    return result

class NbodyProcessor(object):
  def __init__(self,nslices=1,nblocks=None):
    self.nslices=nslices
    if nblocks is None: 
      nblocks=nslices
    n=fsmrt(nblocks)
    ijblocks=(n,nblocks/n)
    self.nblocks=nblocks
    self.ijblocks=ijblocks
  def evaluate(self,func, iparts,jparts,*args, **kwargs ):
    nbunch=kwargs.get('nbunch', ceil(len(iparts),self.nslices))
    if self.pre_pickle:
      jparts=cPickle.dumps(jparts,-1)
      func=eval("pickled"+func.__name__)
    i=0
    while i<len(iparts): 
      job = self.submit_job(func,(iparts[i:i+nbunch],jparts)+args,{})       
      job.range=(i,min(i+nbunch,len(iparts)))
      i=i+nbunch
    result=[None]*len(iparts)
    while self.wait():
      job=self.last_finished_job
      result[job.range[0]:job.range[1]]=job.result
    return result
  def evaluate2(self,func, iparts,jparts,*args, **kwargs ):

    jobreduce=eval("reducers."+func.__name__)

    ibunch=ceil(len(iparts),self.ijblocks[1])
    jbunch=ceil(len(jparts),self.ijblocks[0])
    
    if self.pre_pickle:
      func=eval("pickled2"+func.__name__)
    
    ipart_sets=dict()
    jpart_sets=dict()

    if self.pre_pickle:
      i,j=0,0
      while i<len(iparts):
        ipart_sets[i]=cPickle.dumps(iparts[i:i+ibunch],-1)
        i=i+ibunch
      while j<len(jparts):
        jpart_sets[j]=cPickle.dumps(jparts[j:j+jbunch],-1)
        j=j+jbunch
    else:
      i,j=0,0
      while i<len(iparts):
        ipart_sets[i]=iparts[i:i+ibunch]
        i=i+ibunch
      while j<len(jparts):
        jpart_sets[j]=jparts[j:j+jbunch]
        j=j+jbunch

    i=0
    while i<len(iparts):
      j=0
      while j<len(jparts):
        job = self.submit_job(func,(ipart_sets[i],jpart_sets[j])+args,{})       
        job.range=(i,min(i+ibunch,len(iparts)))
        j=j+jbunch
      i=i+ibunch

    result=[None]*len(iparts)
    while self.wait():
      job=self.last_finished_job
      for i in range(job.range[0],job.range[1]):
        result[i]=jobreduce(result[i],job.result[i-job.range[0]])
    return result

class Processor(processors.Processor,NbodyProcessor):
  def __init__(self,nslices=None,nblocks=None,preamble="pass"):
    processors.Processor.__init__(self,preamble=preamble)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks)

class MultiProcessor(processors.MultiProcessor,NbodyProcessor):
  def __init__(self,preamble="pass",pre_pickle=True,nproc=4,nslices=None,nblocks=None):
    processors.MultiProcessor.__init__(self,preamble=preamble,pre_pickle=pre_pickle,nproc=nproc)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks)

class AmuseProcessor(processors.AmuseProcessor,NbodyProcessor):
  def __init__(self,hosts=[],preamble="pass", pre_pickle=True,channel_type="mpi",verbose=False,
                      nslices=None,nblocks=None):
    processors.AmuseProcessor.__init__(self,hosts=hosts,preamble=preamble, pre_pickle=pre_pickle,
                                         channel_type=channel_type,verbose=verbose)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks)
    
class pp_Processor(processors.pp_Processor,NbodyProcessor):
  def __init__(self,preamble="pass",pre_pickle=True,ppservers=(),nslices=None,nblocks=None):    
    processors.pp_Processor.__init__(self,preamble=preamble,pre_pickle=pre_pickle,
                                       ppservers=ppservers)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks)

class Local_Processor(object):
  def exec_(self,arg):
    pass
  def evaluate(self,func, iparts,*args, **kwargs):
    return func(iparts,*args)
  evaluate2=evaluate

global pproc
pproc=None

class particle(object):
  def __init__(self,m=0,x=0,y=0,z=0,vx=0,vy=0,vz=0):
    self.m=mp.mpf(m)
    self.x=mp.mpf(x)
    self.y=mp.mpf(y)
    self.z=mp.mpf(z)
    self.vx=mp.mpf(vx)
    self.vy=mp.mpf(vy)
    self.vz=mp.mpf(vz)

class jparticle(object):
  def __init__(self,**kwargs):
    for key,value in kwargs.items():
      setattr(self,key,mp.mpf(value))
  
def Particles(N):
    return [particle() for i in range(N)]

def copy_particles(parts):
    copy=Particles(len(parts))
    for i in range(len(parts)):
      copy[i].m=mp.mpf(parts[i].m)
      copy[i].x=mp.mpf(parts[i].x)
      copy[i].y=mp.mpf(parts[i].y)
      copy[i].z=mp.mpf(parts[i].z)
      copy[i].vx=mp.mpf(parts[i].vx)
      copy[i].vy=mp.mpf(parts[i].vy)
      copy[i].vz=mp.mpf(parts[i].vz)
    return copy  

def kinetic_energy(parts):
    return reduce(lambda e,p: e+0.5*p.m*(p.vx**2+p.vy**2+p.vz**2),parts,0.)

def potential_energy(parts):
    global pproc
    pproc.exec_("mp.dps="+str(mp.dps))

    potential(parts,parts)
    return reduce(lambda e,p: e+p.m*p.pot,parts,0.)/2

def total_energy(parts):
    return kinetic_energy(parts)+potential_energy(parts)

def _potential(iparts,jparts):
    result=[]
    for ipart in iparts:
      phi=0.
      for jpart in jparts:
        dx=ipart.x-jpart.x  
        dy=ipart.y-jpart.y  
        dz=ipart.z-jpart.z
        dr2=dx**2+dy**2+dz**2
        if dr2>0:
          dr=mp.sqrt(dr2)
          phi-=jpart.m/dr
      result.append(phi)       
    return result

def pickled_potential(iparts,jparts):
  return _potential(iparts,cPickle.loads(jparts))
def pickled2_potential(iparts,jparts):
  return _potential(cPickle.loads(iparts),cPickle.loads(jparts))

def potential(iparts,jparts):
  result=pproc.evaluate2(_potential, iparts, jparts)
  for ipart,pot in izip(iparts,result):
    ipart.pot=pot
    
def drift(parts,dt):
    for part in parts:
      part.x=part.x+part.vx*dt
      part.y=part.y+part.vy*dt
      part.z=part.z+part.vz*dt
      part.timestep=max_timestep
    return parts  

def _kick(iparts,jparts,dt):
    result=[]
    for ipart in iparts:
      ax=0.
      ay=0.
      az=0.
      for jpart in jparts:
        dx=ipart.x-jpart.x  
        dy=ipart.y-jpart.y  
        dz=ipart.z-jpart.z
        dr2=dx**2+dy**2+dz**2
        if dr2>0:
          dr=mp.sqrt(dr2)
          dr3=dr*dr2
          acci=jpart.m/dr3
          ax-=dx*acci
          ay-=dy*acci
          az-=dz*acci
      result.append([dt*ax,dt*ay,dt*az])
    return result

def pickled_kick(iparts,jparts,dt):
  return _kick(iparts,cPickle.loads(jparts),dt)

def pickled2_kick(iparts,jparts,dt):
  return _kick(cPickle.loads(iparts),cPickle.loads(jparts),dt)

def kick(iparts,jparts,dt):
  iparts_=[jparticle(x=x.x,y=x.y,z=x.z) for x in iparts]
  jparts_=[jparticle(m=x.m,x=x.x,y=x.y,z=x.z) for x in jparts]
  result=pproc.evaluate2(_kick, iparts_, jparts_,dt)
  for ipart,dv in izip(iparts, result):
    ipart.vx+=dv[0]
    ipart.vy+=dv[1]
    ipart.vz+=dv[2]

def _timestep(iparts,jparts, dt_param, rarvratio=1.,max_timestep=1000.):
    sqrt2=mp.sqrt(2.)
    result=[]
    for ipart in iparts:
      timestep=max_timestep
      for jpart in jparts:
        dx=ipart.x-jpart.x  
        dy=ipart.y-jpart.y
        dz=ipart.z-jpart.z
        dr2=dx**2+dy**2+dz**2
        if dr2>0:
          dr=mp.sqrt(dr2)
          dr3=dr*dr2
          dvx=ipart.vx-jpart.vx  
          dvy=ipart.vy-jpart.vy  
          dvz=ipart.vz-jpart.vz
          vdotdr2=(dvx*dx+dvy*dy+dvz*dz)/dr2
          dv2=dvx**2+dvy**2+dvz**2
          mu=ipart.m+jpart.m
          tau=rarvratio*dt_param/sqrt2*mp.sqrt(dr3/mu)
          dtau=3*tau*vdotdr2/2.
          if dtau>1.:
            dtau=1.
          tau/=(1-dtau/2)          
          if tau< timestep:
            timestep=tau
          if dv2>0:
            tau=dt_param*dr/mp.sqrt(dv2)
            dtau=tau*vdotdr2*(1+mu/(dv2*dr))  
            if dtau>1.:
              dtau=1.
            tau/=(1-dtau/2)          
            if tau< timestep:
              timestep=tau
      result.append(timestep)
    return result

def pickled_timestep(iparts,jparts, dt_param, rarvratio=1.,max_timestep=1000.):
  return _timestep(iparts,cPickle.loads(jparts), dt_param, rarvratio=rarvratio,max_timestep=max_timestep)
def pickled2_timestep(iparts,jparts, dt_param, rarvratio=1.,max_timestep=1000.):
  return _timestep(cPickle.loads(iparts),cPickle.loads(jparts), dt_param, rarvratio=rarvratio,max_timestep=max_timestep)

def timestep(iparts,jparts,dt_param, rarvratio=1.,max_timestep=max_timestep):
  result=pproc.evaluate2(_timestep, iparts, jparts,dt_param,rarvratio,max_timestep)
  for ipart,timestep in izip(iparts,result):
    ipart.timestep=timestep

def global_timestep(parts):
  return reduce(lambda x,y: min(x,y.timestep),parts,max_timestep)
      
class bulirschStoer(object):  
  def __init__(self,target_error, dt_param=mp.mpf('1.'), coeff=[1.], 
                 MAXLEVEL=64,jmax=64,rhombus=False,fixed_j=0):
    self.dt_param=dt_param
    self.dtmin=target_error
    self.target_error=target_error
    self.MAXLEVEL=MAXLEVEL
    self.nkickstep=0
    self.nkicks=0
    self.ndriftstep=0
    self.ndrift=0
    self.ntimecalls=0
    self.ntimes=0
    self.nstage=2*len(coeff)-1
    self.coeff=[mp.mpf(x) for x in (coeff+coeff[-2::-1])]
    self.jmax=jmax
    self.jcount=0
    self.nsteps=0
    self.rhombus=rhombus
    self.fixed_j=fixed_j


    self.drift_coeff=[]
    self.kick_coeff=[]
    for j in range(jmax):
      n=self.nsequence(j+1)
      c=n*map(lambda x: x/n,self.coeff)
      d=copy.deepcopy(c)
    
      k=[c[0]/2]
      for i in range(len(c)-1):
        k.append( (c[i]+c[i+1])/2)
      k.append( c[-1]/2)
    
      self.kick_coeff.append(k)
      self.drift_coeff.append(d)
  
  def nsequence(self,j):
    return 2*j
        
  def kdk(self,parts,dt,kick_coeff,drift_coeff):
    kick(parts,parts,kick_coeff[0]*dt)
    for i in range(len(drift_coeff)):
      drift(parts,drift_coeff[i]*dt)  
      kick(parts,parts,kick_coeff[i+1]*dt)    
    self.nkicks+=len(kick_coeff)*len(parts)**2
    self.ndrift+=len(drift_coeff)*len(parts)
    self.nkickstep+=len(kick_coeff)
    self.ndriftstep+=len(drift_coeff)

  def aitkenneville(self,j,k,parts_jk,parts_j1k):
    parts=Particles(len(parts_jk))
    fac=1/((self.nsequence(j)/mp.mpf(self.nsequence(j-k)))**2-1)
    for i in range(len(parts_jk)):
      parts[i].m=parts_jk[i].m
      parts[i].x=parts_jk[i].x+fac*(parts_jk[i].x-parts_j1k[i].x)
      parts[i].y=parts_jk[i].y+fac*(parts_jk[i].y-parts_j1k[i].y)
      parts[i].z=parts_jk[i].z+fac*(parts_jk[i].z-parts_j1k[i].z)
      parts[i].vx=parts_jk[i].vx+fac*(parts_jk[i].vx-parts_j1k[i].vx)
      parts[i].vy=parts_jk[i].vy+fac*(parts_jk[i].vy-parts_j1k[i].vy)
      parts[i].vz=parts_jk[i].vz+fac*(parts_jk[i].vz-parts_j1k[i].vz)
    return parts  

  def rhombusrule(self,j,k,parts_jk,parts_j1k,parts_j1k1):
    parts=Particles(len(parts_jk))
    fac=(self.nsequence(j)/mp.mpf(self.nsequence(j-k)))**2
    for i in range(len(parts_jk)):
      parts[i].m=parts_jk[i].m
      parts[i].x=parts_jk[i].x+(parts_jk[i].x-parts_j1k[i].x)/ \
        (fac*(1-(parts_jk[i].x-parts_j1k[i].x)/(parts_jk[i].x-parts_j1k1[i].x))-1)
      parts[i].y=parts_jk[i].y+(parts_jk[i].y-parts_j1k[i].y)/ \
        (fac*(1-(parts_jk[i].y-parts_j1k[i].y)/(parts_jk[i].y-parts_j1k1[i].y))-1)
      parts[i].z=parts_jk[i].z+(parts_jk[i].z-parts_j1k[i].z)/ \
        (fac*(1-(parts_jk[i].z-parts_j1k[i].z)/(parts_jk[i].z-parts_j1k1[i].z))-1)
      parts[i].vx=parts_jk[i].vx+(parts_jk[i].vx-parts_j1k[i].vx)/ \
        (fac*(1-(parts_jk[i].vx-parts_j1k[i].vx)/(parts_jk[i].vx-parts_j1k1[i].vx))-1)
      parts[i].vy=parts_jk[i].vy+(parts_jk[i].vy-parts_j1k[i].vy)/ \
        (fac*(1-(parts_jk[i].vy-parts_j1k[i].vy)/(parts_jk[i].vy-parts_j1k1[i].vy))-1)
      parts[i].vz=parts_jk[i].vz+(parts_jk[i].vz-parts_j1k[i].vz)/ \
        (fac*(1-(parts_jk[i].vz-parts_j1k[i].vz)/(parts_jk[i].vz-parts_j1k1[i].vz))-1)
    return parts

  def error_function(self,parts1,parts2):
      maxdiv=0.
      for p1,p2 in izip(parts1,parts2):
        maxdiv=max( [maxdiv,  abs(p1.x-p2.x), 
                              abs(p1.y-p2.y),
                              abs(p1.z-p2.z),
                              abs(p1.vx-p2.vx),
                              abs(p1.vy-p2.vy),
                              abs(p1.vz-p2.vz) ])
      return maxdiv

  def evolve_BS(self,parts,dt):

    self.nsteps+=1

#    if dt<self.dtmin:
#      raise Exception

    error=1+2*self.target_error*dt
    j=0
    jline=[]
    j1line=[]
    while error/dt > self.target_error:
      j=j+1
      if j>self.jmax:
        print "fail"
        del jline,j1line
        self.evolve_BS(parts,dt/2)
        self.evolve_BS(parts,dt/2)
        return

      self.jcount+=1

      del j1line
      j1line=jline
      jline=[Particles(len(parts))]
      
      parts1=copy_particles(parts)
      self.kdk(parts1,dt,self.kick_coeff[j-1],self.drift_coeff[j-1])
      jline.append(parts1)
      for k in range(1,j):
        if not self.rhombus:
          jline.append( self.aitkenneville(j,k,jline[k],j1line[k]) )
        else:
          jline.append( self.rhombusrule(j,k,jline[k],j1line[k],j1line[k-1]) )
      if j1line==[]:
        error=1+2*self.target_error*dt
      else:
        error=self.error_function(jline[-1],j1line[-1])
#        error=self.error_function(jline[-1],jline[-2])
      if j==self.fixed_j:
        break
    parts[:]=jline[-1]
#    print 'error:', error,error/dt
    
  def init_evolve(self,parts):
    global pproc
    self.clevel=0
    pproc.exec_("mp.dps="+str(mp.dps))


  def evolve_adapt(self,parts,dt,calctimestep=True):
    self.clevel+=1
    if self.clevel> self.MAXLEVEL:
      print "clevel > MAXLEVEL"
      raise Exception
    if calctimestep:
      for part in parts:
        part.timestep=max_timestep
      timestep(parts,parts,dt_param=self.dt_param)
      self.ntimes+=len(parts)**2
      self.ntimecalls+=1
    dtsys=global_timestep(parts)
    if dtsys < dt:
      self.evolve_adapt(parts,dt/2,calctimestep=False)
      self.evolve_adapt(parts,dt/2,calctimestep=True)
    else:
      self.evolve_BS(parts,dt)
    self.clevel-=1  

  def evolve(self,parts,dt):
    self.init_evolve(parts)    
    self.evolve_adapt(parts,dt)

class floating_point_exact_BS(object):
  def __init__(self, target_error=mp.mpf("1.e-16"),factor=1000,dps=20,**kwargs):
    self.dps_safety=6
    print "fac:",factor
    self.error_factor=factor
    self.initial_target_error=target_error
    self.initial_dps=dps
    self.time=mp.mpf(0.)
    self.particles=[]
    self.checkpoint=[]
    self.checkpoint_time=self.time
    self.kwargs=kwargs
  def commit_particles(self):

    self.target_error1=self.initial_target_error/self.error_factor
    self.target_error2=self.initial_target_error
    mp.dps=self.initial_dps
    while -mp.log10(self.target_error1) > mp.dps - self.dps_safety:
      mp.dps*=2
      print "(re)setting, target error:",float(-mp.log10(self.target_error1)),
      print "digits:",mp.dps  
    self.integrator1=bulirschStoer(self.target_error1,**self.kwargs)
    self.integrator2=bulirschStoer(self.target_error2,**self.kwargs)

    self.checkpoint=copy_particles(self.particles)
    self.particles1=self.particles
    self.particles2=copy_particles(self.particles1)
    self.checkpoint_time=self.time
    
    self.wallclock1=0
    self.wallclock2=0
    
  def recommit_particles(self):
    self.commit_particles()  
  def evolve(self,tend):
    dt=tend-self.time
    
    if dt<=0:
      return

    t1=time.time()
    self.integrator1.evolve(self.particles1,dt)
    t2=time.time()
    self.wallclock1+=t2-t1
    self.integrator2.evolve(self.particles2,dt)
    t3=time.time()
    self.wallclock2+=t3-t2

    print "timing:", t2-t1,t3-t2,", ratio:",(t2-t1)/(t3-t2)

    error=self.error_function(self.particles1,self.particles2)

    while error:

      self.integrator2=self.integrator1
      self.particles2=self.particles1
      self.target_error2=self.target_error1
      self.wallclock2=self.wallclock1
 
      self.target_error1=self.target_error1/self.error_factor
      while -mp.log10(self.target_error1) > mp.dps - self.dps_safety:
        mp.dps*=2
        print "extending, target error:",float(-mp.log10(self.target_error1)),
        print "digits:",mp.dps  
      self.integrator1=bulirschStoer(self.target_error1,**self.kwargs)
      self.particles1=copy_particles(self.checkpoint)
      t1=time.time()
      self.integrator1.evolve(self.particles1,tend-self.checkpoint_time)
      t2=time.time()
      self.wallclock1+=t2-t1
      print "timing (totals):", self.wallclock1,self.wallclock2,self.wallclock1/self.wallclock2
      
      
      error=self.error_function(self.particles1,self.particles2)

    self.particles=self.particles1
    self.time=tend


  def error_function(self,parts1,parts2):
      for p1,p2 in izip(parts1,parts2):
        if ( float(p1.x)!=float(p2.x) or 
             float(p1.y)!=float(p2.y) or
             float(p1.z)!=float(p2.z) or
             float(p1.vx)!=float(p2.vx) or
             float(p1.vy)!=float(p2.vy) or
             float(p1.vz)!=float(p2.vz) ): return True
      return False
