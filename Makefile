BOOST = /usr/local/Cellar/boost/1.66.0  # for unit tests only

CXX      = /usr/local/opt/llvm/bin/clang++
CXXFLAGS = -Wall -pedantic -O3 -fopenmp -I ./include
LDFLAGS  = -L /usr/local/opt/llvm/lib

SRCS := example.cpp $(wildcard include/*.hpp include/*/*.hpp)

run_csa : $(SRCS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) example.cpp -o $@

.PHONY : clean
clean :
	rm -f run_csa
	rm -rf run_csa.dSYM/
