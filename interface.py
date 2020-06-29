import mp
import mp_integrator
from nbody_processors import Local_Processor, MultiProcessor,AmuseProcessor,pp_Processor
import operator

from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.rfi.core import PythonCodeInterface

from amuse.units import units,nbody_system

class YaraviImplementation(object):
    def __init__(self):
        self.particles=dict()
        self._index_counter=0
        self._time=mp.mpf(0.)
        self._begin_time=mp.mpf(0.)
        self.processor="Local_Processor()"
        self.timestep_parameter=mp.mpf('1.')
        self.epsilon_squared=mp.mpf('0.')
        self.initial_target_error=mp.mpf("1.e-16")
        self.factor=mp.mpf("1.e10")
        self.integrator_method="floating_point_exact"
        self.initial_dps=15

    def set_timestep_parameter(self, eta):
        self.timestep_parameter=mp.mpf(eta)
        return 0
    
    def get_timestep_parameter(self,eta):
        eta.value=float(self.timestep_parameter)
        return 0

    def set_eps2(self, eps2):
        self.epsilon_squared=mp.mpf(eps2)
        return 0
    
    def get_eps2(self,eps2):
        eps2.value=float(self.epsilon_squared)
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
        return 0  

    def commit_parameters(self):
        if self.integrator_method=="floating_point_exact":
          self.integrator=mp_integrator.floating_point_exact_BS(
                                         target_error=self.initial_target_error,
                                         factor=self.factor,
                                         dt_param=self.timestep_parameter,
                                         dps=self.initial_dps)
        elif self.integrator_method=="BS":
          self.integrator=mp_integrator.default_BS(
                                         target_error=self.initial_target_error,
                                         dt_param=self.timestep_parameter,
                                         dps=self.initial_dps)
        self._time=self._begin_time
        mp_integrator.pproc=eval(self.processor)
        return 0

    def recommit_parameters(self):
        return -2
      
    def new_particle(self, index, mass,x, y, z, vx, vy, vz, radius):
        id_ = self._index_counter
        self._index_counter+=1
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

    def recommit_particles(self):
        self.integrator.time=self._time
        self.integrator.particles=[]
        for key in sorted(self.particles.keys()):
            self.integrator.particles.append(self.particles[key])
        self.integrator.commit_particles()
        return 0

    def commit_particles(self):
        self.recommit_particles()
        return 0
        
    def evolve_model(self,tend):
        try:
          self.integrator.evolve( mp.mpf(tend))
          for i,key in enumerate(sorted(self.particles.keys())):
            self.particles[key]=self.integrator.particles[i]
          self._time=self.integrator.time
          return 0
        except Exception as ex:
          print ex
          return -1

    def set_begin_time(self,t):
        self._begin_time=mp.mpf(t)
        return 0

    def get_begin_time(self,t):
        t.value=float(self._begin_time)
        return 0

    def set_initial_target_error(self,t):
        self.initial_target_error=mp.mpf(t)
        return 0

    def get_initial_target_error(self,t):
        t.value=float(self.initial_target_error)
        return 0

    def set_factor(self,t):
        self.factor=mp.mpf(t)
        return 0

    def get_factor(self,t):
        t.value=float(self.factor)
        return 0

    def set_initial_dps(self,t):
        self.initial_dps=t
        return 0

    def get_initial_dps(self,t):
        t.value=self.initial_dps
        return 0

    def set_integrator(self,t):
        self.integrator_method=t
        return 0

    def get_integrator(self,t):
        t.value=self.integrator_method
        return 0


    def synchronize_model(self):
        return 0    

    def get_potential_energy(self,e):
        e.value=float(mp_integrator.potential_energy(self.integrator.particles))
        return 0

    def get_kinetic_energy(self,e):
        e.value=float(mp_integrator.kinetic_energy(self.integrator.particles))
        return 0

    def get_total_energy(self,e):
        e.value=float(mp_integrator.kinetic_energy(self.integrator.particles)+
                      mp_integrator.potential_energy(self.integrator.particles))
        return 0

    def get_current_error(self,err):
        err.value=float(self.integrator.current_error())
        return 0

class YaraviInterface(PythonCodeInterface,
                     GravitationalDynamicsInterface, LiteratureReferencesMixIn):
    """
    Yaravi
    
    Yet AnotheR Arithmetic-precision Varying Integrator
    
    Yaravi is a code to solve the astrophysical N-body problem 
    guaranteed to python float( 64-bit floating point) precision. It 
    uses an python implementation of a Bulirsch-Stoer integrator with 
    adaptive timestepping). It has crude support for parallelization. 
    Note that (re)commit_particles calls generate a checkpoint in the 
    code: the solution found is guaranteed to floating point precision 
    from the last check point.

    .. [#] integrator: Boekholt et al., python implementation F.I. Pelupessy

    """

    def __init__(self, **options):
        processor=options.pop("processor")
        
        PythonCodeInterface.__init__(
            self,
            YaraviImplementation,
            'yaravi_worker',
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

    @legacy_function      
    def get_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep_parameter', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep_parameter', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.IN,
          unit=nbody_system.length**2)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_eps2():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.OUT,
          unit=nbody_system.length**2)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_initial_target_error():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_target_error', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function      
    def get_initial_target_error():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_target_error', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('factor', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function      
    def get_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('factor', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_initial_dps():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_dps', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function      
    def get_initial_dps():
        function = LegacyFunctionSpecification()
        function.addParameter('initial_dps', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_integrator():
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='s', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function      
    def get_integrator():
        function = LegacyFunctionSpecification()
        function.addParameter('integrator', dtype='s', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_current_error():
        function = LegacyFunctionSpecification()
        function.addParameter('error_estimate', dtype='d', direction=function.OUT,unit=units.none)
        function.result_type = 'i'
        return function
        
class Yaravi(GravitationalDynamics):

    def __init__(self, convert_nbody=None, **options):
        if not options.has_key("processor"):
          options["processor"]=None
        
        nbody_interface = YaraviInterface(**options)

        GravitationalDynamics.__init__(
            self,
            nbody_interface,
            convert_nbody,
            **options
        )

    def define_properties(self, object):
        GravitationalDynamics.define_properties(self,object)
        object.add_property("get_total_energy")

    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations (not working)", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        
        object.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter", 
            "timestep_parameter", 
            "timestep parameter for gravity calculations", 
            default_value = 1.
        )

        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

        object.add_method_parameter(
            "get_processor", 
            "set_processor",
            "processor", 
            "startup command determining the processor of the integrator (Local_Processor)", 
            default_value = "Local_Processor()"
        )

        object.add_method_parameter(
            "get_initial_target_error",
            "set_initial_target_error", 
            "initial_target_error", 
            "initial_target_error for BS integrator", 
            default_value = 1.e-16
        )
        
        object.add_method_parameter(
            "get_factor",
            "set_factor", 
            "factor", 
            "factor for decrease in target error", 
            default_value = 1.e10
        )        
        
        object.add_method_parameter(
            "get_initial_dps",
            "set_initial_dps", 
            "initial_dps", 
            "initial dps", 
            default_value = 15
        )        

        object.add_method_parameter(
            "get_integrator",
            "set_integrator", 
            "integrator", 
            "integrator (floating_point_exact or BS)", 
            default_value = "floating_point_exact"
        )
        
