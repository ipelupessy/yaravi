from mpmath import mp
import mp_integrator
from mp_integrator import Local_Processor, MultiProcessor,AmuseProcessor,pp_Processor
import operator

from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.rfi.core import PythonCodeInterface

mp.dps=20

class mpNbodyImplementation(object):
    def __init__(self):
        self.particles=dict()
        self.index_counter=0
        self._time=mp.mpf(0.)
        self.processor="Local_Processor()"
        self.timestep_parameter=mp.mpf('1.')

    def set_timestep_parameter(self, eta):
        self.timestep_parameter=eta
        return 0
    
    def get_timestep_parameter(self,eta):
        eta.value=self.timestep_parameter
        return 0
      
    def set_processor(self, processor):
        self.processor=processor
        return 0
    
    def get_processor(self,processor):
        processor.value=self.processor
        return 0

    def get_time(self,time):
        time.value=float(self._time)
        return 0
    
    def initialize_code(self):
        mp_integrator.pproc=eval(self.processor)
        self.integrator=mp_integrator.floating_point_exact_BS(dt_param=self.timestep_parameter)
        self.integrator.time=self._time
        return 0  

    def commit_parameters(self):
        return 0    
      
    def new_particle(self, index, mass,x, y, z, vx, vy, vz, radius):
        id_ = self.index_counter
        self.index_counter+=1
        self.particles[id_]=mp_integrator.particle(mass,x,y,z,vx,vy,vz)
        index.value=id_
        return 0

    def set_state(self,index,mass,x,y,z,vx,vy,vz,radius):
        try:
          p=self.particles[index]
          p.m=mp.mpf(mass)
          p.x=mp.mpf(x)
          p.y=mp.mpf(y)
          p.z=mp.mpf(z)
          p.vx=mp.mpf(vx)
          p.vy=mp.mpf(vy)
          p.vz=mp.mpf(vz)
          return 0
        except:
          return -1

    def get_state(self,index,mass,x,y,z,vx,vy,vz,radius):
        try:
          p=self.particles[index]
          mass.value=float(p.m)
          x.value=float(p.x)
          y.value=float(p.y)
          z.value=float(p.z)
          vx.value=float(p.vx)
          vy.value=float(p.vy)
          vz.value=float(p.vz)
          radius.value=0.
          return 0
        except:
          return -1
          
    def get_position(self,index,x,y,z):
        try:
          p=self.particles[index]
          x.value=float(p.x)
          y.value=float(p.y)
          z.value=float(p.z)
          return 0
        except:
          return -1

    def get_velocity(self,index,vx,vy,vz):
        try:
          p=self.particles[index]
          vx.value=float(p.vx)
          vy.value=float(p.vy)
          vz.value=float(p.vz)
          return 0
        except:
          return -1

    def get_mass(self,index,mass):
        try:
          p=self.particles[index]
          mass.value=float(p.m)
          return 0
        except:
          return -1

    def get_radius(self,index,radius):
        try:
          p=self.particles[index]
          radius.value=0.
          return 0
        except:
          return -1

    def remove_particle(self,index):
        try:
          p=self.particles.pop(index)
          return 0
        except:
          return -1

    def commit_particles(self):
        self.integrator.particles=[]
        for key in sorted(self.particles.keys()):
            self.integrator.particles.append(self.particles[key])
        self.integrator.commit_particles()
        return 0

    def recommit_particles(self):
        self.commit_particles()
        return 0
        
    def evolve_model(self,tend):
        try:
          self.integrator.evolve( mp.mpf(tend))
          for i,key in enumerate(sorted(self.particles.keys())):
            self.particles[key]=self.integrator.particles[i]
          self._time=self.integrator.time
          return 0
        except:
          return -1

    def set_begin_time(self,t):
        self.integrator.time=mp.mpf(t)
        return 0
        
    def synchronize_model(self):
        return 0    

    def get_potential_energy(self,e):
        e.value=float(mp_integrator.potential_energy(self.integrator.particles))
        return 0

    def get_kinetic_energy(self,e):
        e.value=float(mp_integrator.kinetic_energy(self.integrator.particles))
        return 0


class mpNbodyInterface(PythonCodeInterface,
                     GravitationalDynamicsInterface, LiteratureReferencesMixIn):
    """
    mpNbody is a code to solve the astrophysical N-body problem guaranteed to
    floating precision. It uses an python implementation of the Brutus solver
    (Bulirsch-Stoer with adaptive timestepping). It has crude support for 
    parallelization. Note that (re)commit_particles calls generate a checkpoint
    in the code: the solution found is guaranteed to floating point precision from 
    the last check point.

    .. [#] based on Boekholt et al. 

    """

    def __init__(self, **options):
        processor=options.pop("processor")
        
        PythonCodeInterface.__init__(
            self,
            mpNbodyImplementation,
            'mpnbody_worker',
            **options)
        LiteratureReferencesMixIn.__init__(self)
        
        if processor:
          self.set_processor(processor)

    @legacy_function
    def set_processor():
        function = LegacyFunctionSpecification()  
        function.addParameter('processor', dtype='string', direction=function.IN,
            description = "string to generate processor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_processor():
        function = LegacyFunctionSpecification()  
        function.addParameter('processor', dtype='string', direction=function.OUT,
            description = "string to generate processor")
        function.result_type = 'int32'
        return function


class mpNbody(GravitationalDynamics):

    def __init__(self, convert_nbody=None, **options):
        if not options.has_key("processor"):
          options["processor"]=None
        
        nbody_interface = mpNbodyInterface(**options)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )
