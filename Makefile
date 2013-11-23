# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:mpnbody_worker

mpnbody_worker: interface.py mp_nbody.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x amuse.community.mp_nbody.interface mpNbodyInterface mpNbodyImplementation -o $@
	
clean:
	$(RM) -f *.bck *.pyc
	$(RM) -f mpnbody_worker
