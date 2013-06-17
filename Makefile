OBJECT2 = ~/gssl/lib/libgsl.a ~/gssl/lib/libgslcblas.a
OUTFLAG = -openmp  
CFLAG = -fast -I/home/linux/yun9854/gssl/include 

CC= /opt/intel/bin/icc

all: main.out

main.out: main.cpp wksp.h wksp_initial_setting.cpp wksp_h.cpp wksp.cpp
	$(CC) $(OUTFLAG) $(CFLAG) -o main.out main.cpp wksp_initial_setting.cpp wksp_h.cpp wksp.cpp $(OBJECT2)

clean: 
	rm *.o *.out *.png *.data *.txt
