CC=g++

#FLAGS= -ggdb -pg -O3 -static
FLAGS=-O3 -static -fopenmp

SOURCE = align refseq main param reads utilities 
OBJS= $(patsubst %,%.o,$(SOURCE))

all: fr-hit

%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@
fr-hit: $(OBJS)
	$(CC) $(FLAGS) $^ -o $@ 

clean:
	rm -f *.o fr-hit

