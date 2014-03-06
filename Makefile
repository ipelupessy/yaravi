# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:yaravi_worker

yaravi_worker: interface.py mp_integrator.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x amuse.community.yaravi.interface YaraviInterface YaraviImplementation -o $@
	
clean:
	$(RM) -f *.bck *.pyc
	$(RM) -f yaravi_worker
