CC = g++ $(CFLAGS)
CFLAGS = -g -pg -O3 -funroll-all-loops -ffast-math

task: main.o spec_mesh.o spec_mesh.h defs.h
	$(CC) -o task main.o spec_mesh.o

main.o: main.cpp spec_mesh.h defs.h
	$(CC) -c main.cpp

spec_mesh.o: spec_mesh.cpp spec_mesh.h
	$(CC) -c spec_mesh.cpp

clean:
	rm -f *.o *~ task
