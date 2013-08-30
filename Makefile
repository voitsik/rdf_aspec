PROG = rdf_aspec

CC   = gcc
#CC   = icc
CFLAGS = -O3 -march=native -Wall -W -ffast-math
#CFLAGS = -O3 -xAVX -Wall -W
LDFLAGS = -lfftw3f

OBJS = rdf_aspec.o

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(PROG) *.o

.PHONY: all clean
