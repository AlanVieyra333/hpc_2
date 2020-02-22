APP=pi_mpi
PROC=2
CTYPE=c

run: ${APP}.o
	@mpirun -np ${PROC} -mca btl ^openib ${APP}.o

${APP}.o: ${APP}.${CTYPE}
	@mpicc -o ${APP}.o ${APP}.${CTYPE} -O2 -m64 -lm -lc
	@echo "Compilado."
