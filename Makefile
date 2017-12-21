CXXFLAGS=-g -O3 -Wall -std=c++17
LDFLAGS=`icu-config --ldflags` -licuio -lboost_program_options -lgmp -lgmpxx
#CXX=clang++

ufag: ufag.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f *.o ufag
