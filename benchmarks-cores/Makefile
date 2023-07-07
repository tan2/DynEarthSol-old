
CASE := test-3d-tiny.cfg
#CASE := test-tiny-long.cfg
#CASE := test-big.cfg
#CASE := test-topo.cfg
FRAME := 4
OMP := 1
TFLAG = # time -v
DEBUG := 0
NSYS := 0
VTK := 0
NDIMS := 3
EXE := ../dynearthsol$(NDIMS)d
FUNC := compute_mass

GITID = $(shell git rev-parse HEAD | cut -c1-7)

ifeq ($(DEBUG), 1)
	DEBUGFLAG = valgrind --tool=helgrind --log-file=valgrind-out.txt -v
endif
ifeq ($(NSYS), 1)
	DEBUGFLAG = nsys profile --trace=cuda,nvtx,osrt,openmp --output=benchmark.${GITID}.${CASE}.${OMP}
#	NSYSSUM = nsys stats -o benchmark.${GITID}.${CASE}.${OMP}.csv benchmark.${GITID}.${CASE}.${OMP}.qdrep 
endif
ifeq ($(VTK), 1)
	VTKFLAG = ${MAKE} vtk
endif


# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f benchmark.*; OMP_NUM_THREADS=${OMP} ${TFLAG} ${DEBUGFLAG} ${EXE} ${CASE}
RESTART_TEST = OMP_NUM_THREADS=${OMP} ${TFLAG} ${DEBUGFLAG} ${EXE} ${CASE}

ORIG = orig-${CASE}
CURR = ~/data/jobs/${CASE}

all: cmp

auto:
	${MAKE} CASE=test-big.cfg OMP=4
	${MAKE} CASE=test-topo.cfg OMP=4
	${MAKE} CASE=test-model.cfg OMP=4

set:
	# run test, store test output to old/
	${RUN_TEST}
	${MAKE} store

cmp:
	${RUN_TEST}
	${NSYSSUM} > /dev/null
	${VTKFLAG}
	python3 compare.py ${ORIG} ${FRAME}
	cp benchmark.info benchmark.${GITID}.${CASE}.${OMP}.info

restart:
	${RESTART_TEST}
	${NSYSSUM} > /dev/null
	${VTKFLAG}
	python3 compare.py ${ORIG} ${FRAME}
	cp benchmark.info benchmark.${GITID}.${CASE}.${OMP}.info

store:
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv benchmark.* ${ORIG}
	cp ../snapshot.diff ${ORIG}

setjob:
	# st job to compare
	rm -rf orig-job
	ln -s ~/data/jobs/${CASE} orig-job

cmpjob:
	python3 compare.py orig-job ${CURR} ${FRAME}

vtk:
	python ../2vtk.py -m  benchmark

.PHONY: all set cmp store
