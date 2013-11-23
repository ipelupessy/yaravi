from mpmath import mp
import numpy
import copy
import cPickle

from matplotlib import pyplot

import itertools

max_timestep=1000.

class MultiProcessor(object):
  def __init__(self,preamble=None,nbunch=24,pre_pickle=True,nproc=4):
    from multiprocessing import Pool
    self.preamble=preamble
    self.nbunch=nbunch
    self.pre_pickle=pre_pickle
    self.pool=Pool(processes=nproc)
    self.nproc=nproc
  def exec_(self,arg):
    from multiprocessing import Pool
    self.pool.close()
    self.pool.join()
    exec arg
    self.pool=Pool(processes=self.nproc)
  def evaluate(self,func, iparts,jparts,*args, **kwargs ):
    i=0
    nbunch=kwargs.get('nbunch',self.nbunch)
    if self.pre_pickle:
      jparts=cPickle.dumps(jparts,-1)
      orgfunc=(func,)
      func=eval("pickled"+func.__name__)
    else:
      orgfunc=()  
    jobs=[]
    while i<len(iparts): 
      job = self.pool.apply_async(func,(iparts[i:i+nbunch],jparts)+args)       
      job.range=(i,min(i+nbunch,len(iparts)))
      jobs.append(job)
      i=i+nbunch
    for job in jobs:
      iparts[job.range[0]:job.range[1]]=job.get()
    return iparts

class AmuseProcessor(object):
  def __init__(self,hosts=[],preamble=None,nbunch=24,pre_pickle=True,channel_type="mpi",verbose=False):
    from amuse.ext.job_server import JobServer
    self.preamble=preamble
    self.nbunch=nbunch
    self.pre_pickle=pre_pickle
    self.amuse_servers=hosts
    self.job_server=JobServer(self.amuse_servers, channel_type=channel_type,preamble=self.preamble,verbose=verbose)
  def exec_(self,arg):
    self.job_server.exec_(arg)
  def evaluate(self,func, iparts,jparts,*args, **kwargs ):
    i=0
    nbunch=kwargs.get('nbunch',self.nbunch)
    if self.pre_pickle:
      jparts=cPickle.dumps(jparts,-1)
      orgfunc=(func,)
      func=eval("pickled"+func.__name__)
    else:
      orgfunc=()  
    while i<len(iparts): 
      job = self.job_server.submit_job(func,(iparts[i:i+nbunch],jparts)+args,{})       
      job.range=(i,min(i+nbunch,len(iparts)))
      i=i+nbunch
    while self.job_server.wait():
      job=self.job_server.last_finished_job
      iparts[job.range[0]:job.range[1]]=job.result
    return iparts
    

class pp_Processor(object):
  def __init__(self,preamble="pass",nbunch=24,pre_pickle=True):
    import pp
    self.ppservers=()#("galaxy",32),("koppoel",4),("biesbosch",4))#,("gaasp",3))
    if len(self.ppservers)>0:
      from paramiko import SSHClient
    self.openclients()
    self.job_server = pp.Server(ncpus=4,ppservers=tuple(map(lambda x:x[0],self.ppservers)))
    print "Starting pp with", self.job_server.get_ncpus(), "workers and", len(self.ppservers),"clients"
    print "for a total of", self.job_server.get_ncpus()+reduce(lambda x,y: x+y[1],self.ppservers,0)," cpus"
    self.preamble=preamble
    self.nbunch=nbunch
    self.pre_pickle=pre_pickle
  def openclients(self):
    self.clients=[]    
    for ppserver,ncpu in self.ppservers:
      client=SSHClient()
      client.load_system_host_keys()
      client.connect(ppserver)
      stdin,stdout,stderr=client.exec_command('killall -KILL ppserver.py')
      stdin,stdout,stderr=client.exec_command('source /disks/paddegat2/pelupes/amuse/setdev; ppserver.py -w %d'%ncpu)
      self.clients.append(client)
  def __del__(self):
    self.closeclients()
  def closeclients(self):  
    for client in self.clients:
      client.close()  
  def exec_(self,arg):
    self.preamble=self.preamble+";"+arg
  def evaluate(self,func, iparts,jparts,*args, **kwargs ):
    i=0
    nbunch=kwargs.get('nbunch',self.nbunch)
    jobs=[]
    if self.pre_pickle:
      jparts=cPickle.dumps(jparts,-1)
      orgfunc=(func,)
      func=eval("pickled"+func.__name__)
    else:
      orgfunc=()  
    while i<len(iparts): 
      job = self.job_server.submit(func,(iparts[i:i+nbunch],jparts)+args,orgfunc+(particle,jparticle),
                               ("cPickle","mpmath","from mpmath import mp;"+self.preamble)) 
      job.range=(i,min(i+nbunch,len(iparts)))
      jobs.append(job)
      i=i+nbunch
    for job in jobs:   
      result=job()
      iparts[job.range[0]:job.range[1]]=result
      i=i+nbunch
    return iparts
    
class Local_Processor(object):
  def exec_(self,arg):
    pass
  def evaluate(self,func, iparts,*args, **kwargs):
    return func(iparts,*args)

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
  def __init__(self,m=0,x=0,y=0,z=0):
    self.m=mp.mpf(m)
    self.x=mp.mpf(x)
    self.y=mp.mpf(y)
    self.z=mp.mpf(z)
  
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
      ipart.pot=phi       
    return iparts

def pickled_potential(iparts,jparts):
  return _potential(iparts,cPickle.loads(jparts))

def potential(iparts,jparts):
  return pproc.evaluate(_potential, iparts, jparts)
    
def drift(parts,dt):
    for part in parts:
      part.x=part.x+part.vx*dt
      part.y=part.y+part.vy*dt
      part.z=part.z+part.vz*dt
      part.timestep=max_timestep
    return parts  

def _kick(iparts,jparts,dt):
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
      ipart.vx+=dt*ax       
      ipart.vy+=dt*ay       
      ipart.vz+=dt*az       
    return iparts

def pickled_kick(iparts,jparts,dt):
  return _kick(iparts,cPickle.loads(jparts),dt)

def kick(iparts,jparts,dt):
  jparts=[jparticle(x.m,x.x,x.y,x.z) for x in jparts]
  return pproc.evaluate(_kick, iparts, jparts,dt)

def _timestep(iparts,jparts, dt_param, rarvratio=1.,max_timestep=1000.):
    sqrt2=mp.sqrt(2.)
    for ipart in iparts:
      timestep=max_timestep
      for jpart in jparts:
        if ipart==jpart:
          continue
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
      if timestep < ipart.timestep:
        ipart.timestep=timestep
    return iparts

def pickled_timestep(iparts,jparts, dt_param, rarvratio=1.,max_timestep=1000.):
  return _timestep(iparts,cPickle.loads(jparts), dt_param, rarvratio=rarvratio,max_timestep=max_timestep)


def timestep(iparts,jparts,dt_param, rarvratio=1.,max_timestep=max_timestep):
  return pproc.evaluate(_timestep, iparts, jparts,dt_param,rarvratio,max_timestep)

def global_timestep(parts):
  return reduce(lambda x,y: min(x,y.timestep),parts,max_timestep)

class shared_2nd_kdk(object):

  def __init__(self,dt_param,MAXLEVEL=64):
    self.dt_param=dt_param
    self.MAXLEVEL=MAXLEVEL
    self.nkickstep=0
    self.nkicks=0
    self.ndriftstep=0
    self.ndrift=0
    self.ntimecalls=0
    self.ntimes=0
    
  def kdk(self,parts,dt):
    kick(parts,parts,dt/2)
    drift(parts,dt)
    kick(parts,parts,dt/2)
    self.nkicks+=2*len(parts)*len(parts)
    self.ndrift+=len(parts)
    self.nkickstep+=2
    self.ndriftstep+=1

  def init_evolve(self,parts):
    self.clevel=0
    for part in parts:
      part.timestep=max_timestep

  def evolve_shared2(self,parts,dt,calctimestep=True):
    self.clevel+=1
    if self.clevel> self.MAXLEVEL:
      print "clevel > MAXLEVEL"
      raise Exception
    if calctimestep:
      timestep(parts,parts,dt_param=self.dt_param)
      self.ntimes+=len(parts)**2
      self.ntimecalls+=1
    dtsys=global_timestep(parts)
    if dtsys < dt:
      self.evolve_shared2(parts,dt/2,calctimestep=False)
      self.evolve_shared2(parts,dt/2,calctimestep=True)
    else:
      self.kdk(parts,dt)
    self.clevel-=1  
  def evolve(self,parts,dt):
    self.init_evolve(parts)    
    self.evolve_shared2(parts,dt)

class shared_4nd_kdk(object):

  K1=(642+mp.sqrt(471))/3924 
  K2=121*(12 - mp.sqrt( 471 ) )/3924 
  K3=1-2*(K1+K2)
  D1=mp.mpf(6)/11 
  D2=0.5-D1  

  def __init__(self,dt_param,MAXLEVEL=64):
    self.dt_param=dt_param
    self.MAXLEVEL=MAXLEVEL
    self.nkickstep=0
    self.nkicks=0
    self.ndriftstep=0
    self.ndrift=0
    self.ntimecalls=0
    self.ntimes=0
    
  def kdk4(self,parts,dt):
    kick(parts,parts,self.K1*dt)
    drift(parts,self.D1*dt)
    kick(parts,parts,self.K2*dt)
    drift(parts,self.D2*dt)
    kick(parts,parts,self.K3*dt)
    drift(parts,self.D2*dt)
    kick(parts,parts,self.K2*dt)
    drift(parts,self.D1*dt)
    kick(parts,parts,self.K1*dt)
    self.ndrift+=4*len(parts)
    self.nkicks+=5*len(parts)*len(parts)
    self.nkickstep+=5
    self.ndriftstep+=4

  def init_evolve(self,parts):
    self.clevel=0
    for part in parts:
      part.timestep=max_timestep

  def evolve_shared4(self,parts,dt,calctimestep=True):
    self.clevel+=1
    if self.clevel> self.MAXLEVEL:
      print "clevel > MAXLEVEL"
      raise Exception
    if calctimestep:
      timestep(parts,parts,dt_param=self.dt_param)
      self.ntimes+=len(parts)**2
      self.ntimecalls+=1
    dtsys=global_timestep(parts)
    if dtsys < dt:
      self.evolve_shared4(parts,dt/2,calctimestep=False)
      self.evolve_shared4(parts,dt/2,calctimestep=True)
    else:
      self.kdk4(parts,dt)
    self.clevel-=1  

  def evolve(self,parts,dt):
    self.init_evolve(parts)    
    self.evolve_shared4(parts,dt)

class symmetric_composition(object):

  def __init__(self,dt_param, coeff, MAXLEVEL=64):
    self.dt_param=dt_param
    self.MAXLEVEL=MAXLEVEL
    self.nkickstep=0
    self.nkicks=0
    self.ndriftstep=0
    self.ndrift=0
    self.ntimecalls=0
    self.ntimes=0
    
    self.drift_coeff=[]
    self.nstage=2*len(coeff)-1
    for i in range(self.nstage):
      self.drift_coeff.append( coeff[min(i,self.nstage-i-1)] )
    
    self.kick_coeff=[coeff[0]/2]
    for i in range(self.nstage-1):
      self.kick_coeff.append( (coeff[min(i,self.nstage-i-1)]+coeff[min(i+1,self.nstage-i-2)])/2)
    self.kick_coeff.append( coeff[0]/2)
    
    print len(self.drift_coeff),sum(self.drift_coeff)
    print len(self.kick_coeff),sum(self.kick_coeff)
        
  def kdk(self,parts,dt):
    kick(parts,parts,self.kick_coeff[0]*dt)
    for i in range(len(self.drift_coeff)):
      drift(parts,self.drift_coeff[i]*dt)  
      kick(parts,parts,self.kick_coeff[i+1]*dt)    
    self.nkicks+=len(self.kick_coeff)*len(parts)**2
    self.ndrift+=len(self.drift_coeff)*len(parts)
    self.nkickstep+=len(self.kick_coeff)
    self.ndriftstep+=len(self.drift_coeff)

  def init_evolve(self,parts):
    self.clevel=0
    for part in parts:
      part.timestep=max_timestep

  def evolve_adapt(self,parts,dt,calctimestep=True):
    self.clevel+=1
    if self.clevel> self.MAXLEVEL:
      print "clevel > MAXLEVEL"
      raise Exception
    if calctimestep:
      timestep(parts,parts,dt_param=self.dt_param)
      self.ntimes+=len(parts)**2
      self.ntimecalls+=1
    dtsys=global_timestep(parts)
    if dtsys < dt:
      self.evolve_adapt(parts,dt/2,calctimestep=False)
      self.evolve_adapt(parts,dt/2,calctimestep=True)
    else:
      self.kdk(parts,dt)
    self.clevel-=1  

  def evolve(self,parts,dt):
    self.init_evolve(parts)    
    self.evolve_adapt(parts,dt)
      
class bulirschStoer(object):  
  def __init__(self,target_error, dt_param=mp.mpf('.125'), coeff=[1.], MAXLEVEL=64,jmax=64,rhombus=False,fixed_j=0):
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
    
#      print len(d),sum(d)
#      print len(k),sum(k)
      self.kick_coeff.append(k)
      self.drift_coeff.append(d)
  
  def nsequence(self,j):
    return 2*j
#[1,2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
# 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163
#, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251
#, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349
#, 353, 359, 367, 373, 379, 383, 389, 397][j]
        
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
      for i in range(len(parts1)):
        maxdiv=max( [maxdiv, abs(parts1[i].x-parts2[i].x),
                             abs(parts1[i].y-parts2[i].y),
        abs(parts1[i].z-parts2[i].z),
        abs(parts1[i].vx-parts2[i].vx),
        abs(parts1[i].vy-parts2[i].vy),
        abs(parts1[i].vz-parts2[i].vz)])
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

class error_controlled_BS(object):
  def __init__(self,target_error,factor=1000,**kwargs):
    self.global_target_error=target_error
    self.error_factor=factor
    self.target_error1=target_error
    self.target_error2=target_error*factor
    self.integrator1=bulirschStoer(self.target_error1,**kwargs)
    self.integrator2=bulirschStoer(self.target_error2,**kwargs)
    self.time=mp.mpf(0.)
    self.particles=[]
    self.checkpoint=[]
    self.checkpoint_time=self.time
    self.dps_safety=10
    self.kwargs=kwargs
  def commit_particles(self):
    self.checkpoint=copy_particles(self.particles)
    self.particles1=self.particles
    self.particles2=copy_particles(self.particles1)
    self.checkpoint_time=self.time
  def recommit_particles(self):
    self.commit_particles()  
  def evolve(self,tend):
    dt=tend-self.time
    
    if dt<=0:
      return
    self.integrator1.evolve(self.particles1,dt)
    self.integrator2.evolve(self.particles2,dt)

    error=self.error_function(self.particles1,self.particles2)

    while error > self.global_target_error:

      self.integrator2=self.integrator1
      self.particles2=self.particles1
      self.target_error2=self.target_error1
 
      self.target_error1=self.target_error1/self.error_factor
      if -mp.log10(self.target_error1) > mp.dps - self.dps_safety:
        mp.dps*=2
      print "extending precision:",-mp.log10(self.target_error1),mp.dps  
      self.integrator1=bulirschStoer(self.target_error1,**self.kwargs)
      self.particles1=copy_particles(self.checkpoint)
      self.integrator1.evolve(self.particles1,tend-self.checkpoint_time)
      
      error=self.error_function(self.particles1,self.particles2)

    self.particles=self.particles1
    self.time=tend

  def error_function(self,parts1,parts2):
      maxdiv=0.
      for i in range(len(parts1)):
        maxdiv=max( [maxdiv, abs(parts1[i].x-parts2[i].x),
                             abs(parts1[i].y-parts2[i].y),
        abs(parts1[i].z-parts2[i].z),
        abs(parts1[i].vx-parts2[i].vx),
        abs(parts1[i].vy-parts2[i].vy),
        abs(parts1[i].vz-parts2[i].vz)])
      return maxdiv
