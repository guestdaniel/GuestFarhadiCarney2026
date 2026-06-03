gcc -c -Wall demo.c 
gcc -c -Wall complex.c 
gcc -shared -o demo.so demo.o complex.o
