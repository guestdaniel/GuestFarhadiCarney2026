# This script includes the core compilation steps for the Zilany, Bruce, and Carney (2014) model
# In practice, these are executed via BinaryBuilder in the build_tarballs.jl script, but it could be informative to see them as a shell script so, here they are:

# Step 1: Compile model_IHC, model_Syanpse, complex, and test with warnings (-Wall), position independent code (-fPIC), and optimization level three (-O3)
gcc -c -Wall -fPIC -O3 model_IHC_debug.c
gcc -c -Wall -fPIC -O3 model_Synapse_debug.c
gcc -c -Wall -fPIC -O3 complex.c 

# Step 2: Compile everything together into a shared library object
gcc -shared -o libzbc2014debug.so model_IHC_debug.o model_Synapse_debug.o complex.o
