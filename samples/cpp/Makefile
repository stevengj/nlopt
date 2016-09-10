CXX = g++

# you may need to set PKG_CONFIG_PATH to the directory containing nlopt.pc
CXXFLAGS = -Wall $(shell pkg-config --cflags nlopt)
LDFLAGS = $(shell pkg-config --libs nlopt)

all: example_cpp

# you may need to set LD_LIBRARY_PATH to the directory containing the
# NLopt shared library when running the example
example_cpp: example.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f *.o *.exe example_cpp
