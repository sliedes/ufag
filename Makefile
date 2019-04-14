CXXFLAGS=`pkg-config --cflags icu-io` -g -O3 -Wall -std=c++17
LDFLAGS=`pkg-config --libs icu-io` -licuio -lboost_program_options
#CXX=clang++

ufag: ufag.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f *.o ufag
