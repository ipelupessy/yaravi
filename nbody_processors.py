import cPickle
import processors

# naieve function to find biggest factor smaller than root
def fsmrt(N):
  n=int(N**0.5)
  while N%n:
    n-=1
  return n

ceil=lambda x,y: (x/y+(x%y>0))

class pickled_argument(object):
  def __init__(self,func,pickled=(True,True)):
    self.func=func
    self.pickled=pickled
  def __call__(self,arg1,arg2,*args,**kwargs):
    return self.func(cPickle.loads(arg1) if self.pickled[0] else arg1,
                     cPickle.loads(arg2) if self.pickled[1] else arg2, 
                     *args,**kwargs)

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
  def __init__(self,nslices=1,nblocks=None,pre_pickle=False):
    self.nslices=nslices
    if nblocks is None: 
      nblocks=nslices
    self.pre_pickle=pre_pickle
    n=fsmrt(nblocks)
    ijblocks=(n,nblocks/n)
    self.nblocks=nblocks
    self.ijblocks=ijblocks
  def evaluate(self,func, iparts,jparts,*args, **kwargs ):
    nbunch=kwargs.get('nbunch', ceil(len(iparts),self.nslices))
    if self.pre_pickle:
      jparts=cPickle.dumps(jparts,-1)
      func=pickled_argument(func,(False,True))
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
      func=pickled_argument(func)
    
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
    processors.MultiProcessor.__init__(self,preamble=preamble,nproc=nproc)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks,pre_pickle=pre_pickle)

class AmuseProcessor(processors.AmuseProcessor,NbodyProcessor):
  def __init__(self,hosts=[],preamble="pass", pre_pickle=True,channel_type="mpi",verbose=False,
                      nslices=None,nblocks=None):
    processors.AmuseProcessor.__init__(self,hosts=hosts,preamble=preamble,
                                         channel_type=channel_type,verbose=verbose)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks,pre_pickle=pre_pickle)
    
class pp_Processor(processors.pp_Processor,NbodyProcessor):
  def __init__(self,preamble="pass",pre_pickle=True,ppservers=(),nslices=None,nblocks=None):    
    processors.pp_Processor.__init__(self,preamble=preamble, ppservers=ppservers)
    if nslices is None: nslices=self.nproc
    if nblocks is None: nblocks=nslices
    NbodyProcessor.__init__(self,nslices=nslices,nblocks=nblocks,pre_pickle=pre_pickle)

class Local_Processor(object):
  def exec_(self,arg):
    pass
  def evaluate(self,func, iparts,*args, **kwargs):
    return func(iparts,*args)
  evaluate2=evaluate
