CC = gcc

# you may need to set PKG_CONFIG_PATH to the directory containing nlopt.pc
CFLAGS = -Wall $(shell pkg-config --cflags nlopt)
LDFLAGS = $(shell pkg-config --libs nlopt)

all: example_c

# you may need to set LD_LIBRARY_PATH to the directory containing the
# NLopt shared library when running the example
example_c: example.c
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *.o *.exe example_c
