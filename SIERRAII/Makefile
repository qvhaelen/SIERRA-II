CC=g++
CFLAGS= -c  -O3
LDFLAGS= 

objects := $(patsubst %.cc,%.o,$(wildcard *.cc))


all: sierraII

sierraII: $(objects)
	$(CC)  -o sierraII $(objects)

%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o 
