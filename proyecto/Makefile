APP=strassen_mpi_omp
ARGS=matrix_1000.txt matrix_1000.txt
NUM_THREADS=4
PROC=2

CC_FLAGS=-O2 -m64 -lm -lc -fopenmp
RUN_FLAGS=-np ${PROC} -mca btl ^openib
CLANG=c
MPICC=mpicc
#MPICC=/usr/lib64/openmpi/1.4-gcc/bin/mpicc
MPIRUN=mpirun
#MPIRUN=/usr/lib64/openmpi/1.4-gcc/bin/mpirun

.PHONY: run clean
run: ${APP}.o
	bash -c "export OMP_NUM_THREADS=${NUM_THREADS} && ${MPIRUN} ${RUN_FLAGS} ${APP}.o ${ARGS} --bind-to-socket"

${APP}.o: ${APP}.${CLANG}
	${MPICC} -o ${APP}.o ${APP}.${CLANG} ${CC_FLAGS}
	@echo "Compilado."

clean:
	@rm *.o