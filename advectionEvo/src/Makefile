obj = main.o evolve.o init.o grid.o cell.o

all: $(obj)
	nvcc $(obj) -o app

%.o: %.cu
	nvcc -x cu -I. -dc $< -o $@

clean:
	rm -f *.o app
