B1;2cCC=gcc
CXX=g++

FLAGS=-O3 



all: optimSFS supersim sfstools emOptim calcStat bgid

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

bgid:	bgid5.cpp stats.cpp ../kstring.o ../bgzf.o
	$(CXX) $(FLAGS) bgid5.cpp -I../  ../kstring.o ../bgzf.o ../knetfile.o -lz -o bgid


clean:
	rm  -f *.o optimSFS simnextgen supersim emOptim calcStat sfstools bgid *~
