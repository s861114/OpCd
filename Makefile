OBJECT2 = ~/gssl/lib/libgsl.a ~/gssl/lib/libgslcblas.a
OUTFLAG = -fopenmp  
CFLAG = -O3 -lm 
CFLAGI = -O3 -I/usr/local/include 

CC= g++

all: main.out

main.out: main.o wksp_initial_setting.o wksp_h.o wksp.o
	$(CC) $(OUTFLAG) $(CFLAG) -g -o main.out main.o wksp_initial_setting.o wksp_h.o wksp.o /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a

main.o: main.cpp wksp.h
	$(CC) $(CFLAGI) -c main.cpp

wksp_initial_setting.o: wksp_initial_setting.cpp
	$(CC) $(CFLAGI) -c wksp_initial_setting.cpp

wksp_h.o: wksp_h.cpp
	$(CC) $(CFLAGI) -c wksp_h.cpp

wksp.o: wksp.cpp
	$(CC) $(CFLAGI) -c wksp.cpp


clean: 
	rm *.o *.out *.png *.data *.txt
