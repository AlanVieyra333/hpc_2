APP=strassen

.PHONY: run clear

build: ${APP}.o generateMatrix.o

run: build
	./${APP}.o matrix_a.txt matrix_b.txt

${APP}.o: ${APP}.cpp
	g++ -std=c++11 -o ${APP}.o ${APP}.cpp

generateMatrix.o: generateMatrix.c
	gcc -o generateMatrix.o generateMatrix.c

clear:
	rm -f -R *.o
