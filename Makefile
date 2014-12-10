# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

MPMATH_AVAILABLE := $(shell $(PYTHON) -c "import mpmath"  1>&2 2> /dev/null && echo "yes" || echo "no")
GMPY2_AVAILABLE := $(shell $(PYTHON) -c "import gmpy2"  1>&2 2> /dev/null && echo "yes" || echo "no")

all: test yaravi_worker

yaravi_worker: interface.py mp_integrator.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x amuse.community.yaravi.interface YaraviInterface YaraviImplementation -o $@
	
test:
	@echo
	@echo "Testing import of mpmath:"
ifeq ($(MPMATH_AVAILABLE),no)
	$(error "Python import mpmath not available - install this first")
endif
	@echo "Testing import of gmpy2:"
ifeq ($(GMPY2_AVAILABLE),no)
	$(error "Python import gmpy2 not available - install this first")
endif
	@echo "Tests successful!"
	@echo
	
	
clean:
	$(RM) -f *.bck *.pyc
	$(RM) -f yaravi_worker
