bpnet: main.o libbpfind.o 
	g++  -static main.o libbpfind.o -o bpnet.linux -lgfortran -lquadmath
libbpfind.o: libbpfind.f
	gfortran -c -static libbpfind.f 
main.o: main.cpp
	g++ -c -static main.cpp -I./include
	
.PHONY : clean
clean :
	-rm -f *.o 

