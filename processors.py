from collections import deque
from functools import partial

"""
abstraction of job server functionality:
local (not parallel), multiprocessing, amuse and parallel python versions 

example:
processor=Processor()

processor.exec_("<set some global var>")

processor.submit_job( f, args, kwargs)
..

while processor.wait():
  result=processor.last_finished_job.result
  <do something with result>


"""

#trivial local processor
class Processor(object):
  def __init__(self,preamble="pass",pre_pickle=False):
    self.nproc=1
    self._jobs=deque()
    self.pre_pickle=pre_pickle
    self.context=dict()
    exec preamble in self.context
  def submit_job(self,f,args=(),kwargs={}):
    job=partial(f,*args,**kwargs)
    self._jobs.append(job)
    return job
  def exec_(self,arg):
    exec arg in self.context
  def wait(self):
    if self._jobs:
      job=self._jobs.popleft()
      self.context["_func"]=job
      job.result=eval("_func()",self.context)
      self._last_finished_job=job
      return True
    else:
      return False 
  @property
  def last_finished_job(self):
    return self._last_finished_job

class MultiProcessor(object):
  def __init__(self,preamble="pass",pre_pickle=True,nproc=4):
    from multiprocessing import Pool
    self.preamble=preamble
    self.pre_pickle=pre_pickle
    exec self.preamble
    self.pool=Pool(processes=nproc)
    self.nproc=nproc
    self._jobs=deque()
    self._last_finished_job=None
  def exec_(self,arg):
    from multiprocessing import Pool
    self.pool.close()
    self.pool.join()
    exec self.preamble
    exec arg
    self.pool=Pool(processes=self.nproc)
  def submit_job(self,f,args=(),kwargs={}):
    job=self.pool.apply_async(f,args)
    self._jobs.append(job)
    return job
  def wait(self):
    if self._jobs:
      job=self._jobs.popleft()
      job.result=job.get()
      self._last_finished_job=job
      return True
    else:
      return False 
  @property
  def last_finished_job(self):
    return self._last_finished_job       

class AmuseProcessor(object):
  def __init__(self,hosts=[],preamble="pass", pre_pickle=True,channel_type="mpi",verbose=False):
    from amuse.ext.job_server import JobServer
    self.preamble=preamble
    self.pre_pickle=pre_pickle
    self.amuse_servers=hosts
    self.job_server=JobServer(self.amuse_servers, channel_type=channel_type,
      preamble=self.preamble,verbose=verbose,no_wait=True)
  @property
  def nproc(Self):
    return self.job_server.number_available_codes  
  def exec_(self,arg):
    self.job_server.exec_(arg)
  def submit_job(self,f,args=(),kwargs={}):
    return self.job_server.submit_job(f,args,kwargs)
  def wait(self):
    return self.job_server.wait()
  @property
  def last_finished_job(self):
    return self.job_server.last_finished_job       
    
# ppservers should be tuple with hostnames and #cpu:
# ppserver=(("galaxy",32),("koppoel",4),("biesbosch",4),("gaasp",3))
# note hardcoded strw environment
class pp_Processor(object):
  def __init__(self,preamble="pass",pre_pickle=True,ppservers=(),depfuncs=()):
    import pp
    self.ppservers=()
    if len(self.ppservers)>0:
      from paramiko import SSHClient
    self.openclients()
    self.job_server = pp.Server(ncpus=4,ppservers=tuple(map(lambda x:x[0],self.ppservers)))
    ncpu=self.job_server.get_ncpus()+reduce(lambda x,y: x+y[1],self.ppservers,0)
    print "Starting pp with", self.job_server.get_ncpus(), "workers and", len(self.ppservers),"clients"
    print "for a total of", ncpu," cpus"
    self.nproc=ncpu
    self.preamble=preamble
    self.pre_pickle=pre_pickle
    self._jobs=deque()
    self._last_finished_job=None
    self.depfuncs=depfuncs
    self.modules=("cPickle","mp","mpmath","gmpy2")
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
    self.preamble_=self.preamble+";"+arg
  def submit_job(self,f,args=(),kwargs={}):
    job=self.job_server.submit(f,args,self.depfuncs,self.modules+(self.preamble_,))
    self._jobs.append(job)
    return job
  def wait(self):
    if self._jobs:
      job=self._jobs.popleft()
      job.result=job()
      self._last_finished_job=job
      return True
    else:
      return False 
  @property
  def last_finished_job(self):
    return self._last_finished_job
