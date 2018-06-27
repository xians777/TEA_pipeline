BAMTOOLS=/home/jc624/eclipse/workspaces/genomics-new/mybamtools/src/
TEA_INCLUDE=/home/jc624/eclipse/workspaces/genomics-new/include/

BAMTOOLS_LIB=/home/jc624/eclipse/workspaces/genomics-new/mybamtools/lib
TEA_LD=/home/jc624/eclipse/workspaces/genomics-new/lib/

LIBS=-lboost_filesystem -lboost_system -lboost_regex -lboost_thread -lbamtools-utils -lbamtools -ljsoncpp -lm -lrt -lz

CC = g++
CFLAGS = -O3 -Wall  -static -std=c++11 -pthread -fopenmp -I$(BAMTOOLS) -I$(TEA_INCLUDE)  -L$(BAMTOOLS_LIB) -L$(TEA_LD)
GCC = gcc
GCFLAGS = -O3 -Wall -static  -pthread -fopenmp -I$(BAMTOOLS) -I$(TEA_INCLUDE) -L$(BAMTOOLS_LIB) -L$(TEA_LD)

all: tea

xfaidx.o : ../src/third/xfaidx.cpp ../src/third/xfaidx.h
	$(CC) $(CFLAGS) -c ../src/third/xfaidx.cpp
faidx.o : ../src/third/faidx.c ../src/third/faidx.h
	$(GCC) $(GCFLAGS) -c ../src/third/faidx.c
razf.o : ../src/third/razf.c ../src/third/razf.h
	$(GCC) $(GCFLAGS) -c ../src/third/razf.c
gzstream.o : ../src/third/gzstream.cpp ../src/third/gzstream.h
	$(CC) $(CFLAGS) -c ../src/third/gzstream.cpp
ResultCompartor.o : ../src/tea/ResultCompartor.cpp ../src/tea/ResultCompartor.hpp
	$(CC) $(CFLAGS) -c ../src/tea/ResultCompartor.cpp
TEA.o : ../src/tea/TEA.cpp ../src/tea/TEA.hpp
	$(CC) $(CFLAGS) -c ../src/tea/TEA.cpp
TEAOptionParser.o : ../src/tea/TEAOptionParser.cpp ../src/tea/TEAOptionParser.hpp
	$(CC) $(CFLAGS) -c ../src/tea/TEAOptionParser.cpp
ReadGroup.o : ../src/meerkat/preprocess/ReadGroup.cpp ../src/meerkat/preprocess/ReadGroup.h
	$(CC) $(CFLAGS) -c ../src/meerkat/preprocess/ReadGroup.cpp
BWACaller.o : ../src/meerkat/BWACaller.cpp ../src/meerkat/BWACaller.hpp
	$(CC) $(CFLAGS) -c ../src/meerkat/BWACaller.cpp
BinaryEncoder.o : ../src/castle/BinaryEncoder.cpp ../src/castle/BinaryEncoder.hpp
	$(CC) $(CFLAGS) -c ../src/castle/BinaryEncoder.cpp
IOUtils.o : ../src/castle/IOUtils.cpp ../src/castle/IOUtils.hpp
	$(CC) $(CFLAGS) -c ../src/castle/IOUtils.cpp
MathUtils.o : ../src/castle/MathUtils.cpp ../src/castle/MathUtils.hpp
	$(CC) $(CFLAGS) -c ../src/castle/MathUtils.cpp
MurmurHash3.o : ../src/castle/MurmurHash3.cpp ../src/castle/MurmurHash3.hpp
	$(CC) $(CFLAGS) -c ../src/castle/MurmurHash3.cpp
OptionParser.o : ../src/castle/OptionParser.cpp ../src/castle/OptionParser.hpp
	$(CC) $(CFLAGS) -c ../src/castle/OptionParser.cpp
StringUtils.o : ../src/castle/StringUtils.cpp ../src/castle/StringUtils.hpp
	$(CC) $(CFLAGS) -c ../src/castle/StringUtils.cpp
ParallelRunner.o : ../src/castle/ParallelRunner.cpp ../src/castle/ParallelRunner.hpp
	$(CC) $(CFLAGS) -c ../src/castle/ParallelRunner.cpp
TimeChecker.o : ../src/castle/TimeChecker.cpp ../src/castle/TimeChecker.hpp
	$(CC) $(CFLAGS) -c ../src/castle/TimeChecker.cpp
tea.o : ../src/tea.cpp
	$(CC) $(CFLAGS) -c ../src/tea.cpp
tea : xfaidx.o faidx.o razf.o gzstream.o ResultCompartor.o TEA.o TEAOptionParser.o ReadGroup.o BWACaller.o BinaryEncoder.o IOUtils.o MathUtils.o MurmurHash3.o OptionParser.o ParallelRunner.o  StringUtils.o TimeChecker.o tea.o
	$(CC) $(CFLAGS) -o TEA_mem1 xfaidx.o faidx.o razf.o gzstream.o ResultCompartor.o TEA.o TEAOptionParser.o ReadGroup.o BWACaller.o BinaryEncoder.o IOUtils.o MathUtils.o MurmurHash3.o OptionParser.o StringUtils.o ParallelRunner.o TimeChecker.o tea.o $(LIBS)
