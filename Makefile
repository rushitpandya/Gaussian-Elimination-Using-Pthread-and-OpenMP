ll: openmp pthread

openmp:
	gcc gauss_openmp.c -std=gnu99 -fopenmp -lm -o openmp

pthread:
	gcc gauss_pthread.c -std=gnu99 -lpthread -lm -o pthread

clean:
	rm openmp pthread
