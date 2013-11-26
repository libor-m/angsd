CC=gcc
CXX=g++

FLAGS=-O3 


all: supersim sfstools calcStat bgid emOptim2 msToGlf

bfgs.o: bfgs.c bfgs.h
	$(CC) $(FLAGS) -c bfgs.c

optimSFS: testfolded18.c bfgs.o
	$(CC) $(FLAGS) testfolded18.c -std=c99 bfgs.o -o optimSFS -lm -lpthread

simnextgen: simnextgen10.c
	$(CC) $(FLAGS) simnextgen10.c -o simnextgen -lm -lz

supersim: supersim.cpp rbeta.cpp
	$(CXX) $(FLAGS) supersim.cpp -lz -o supersim

calcStat: calcStat24.cpp
	$(CXX) $(FLAGS) calcStat24.cpp -lz -o calcStat

sfstools: sfstools.cpp
	$(CXX) $(FLAGS) sfstools.cpp -o sfstools -lz

emOptim: emOptim.cpp
	$(CXX) $(FLAGS) emOptim.cpp -o emOptim -lz -lpthread

emOptim2: emOptim2.cpp
	$(CXX) $(FLAGS) emOptim2.cpp bfgs.cpp -o emOptim2 -lz -lpthread


bgid:	bgid6.cpp stats.cpp ../kstring.o ../bgzf.o
	$(CXX) $(FLAGS) bgid6.cpp -I../  ../kstring.o ../bgzf.o ../knetfile.o -lz -o bgid



msToGlf: msToGlf6.c
	$(CC) $(FLAGS) msToGlf6.c -O3 -o msToGlf -std=gnu99 -lz -lm



clean:
	rm  -f *.o optimSFS simnextgen supersim emOptim calcStat sfstools bgid emOptim2 *~
