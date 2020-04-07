APP=strassen

.PHONY: run clear

build: ${APP}.o generateMatrix.o

run: ${APP}.o
	./${APP}.o matrix_a.txt matrix_b.txt

${APP}.o: ${APP}.cpp
	g++ -std=c++11 -o ${APP}.o ${APP}.cpp -O2

generateMatrix.o: generateMatrix.c
	gcc -o generateMatrix.o generateMatrix.c -O2

clear:
	rm -f -R *.o
