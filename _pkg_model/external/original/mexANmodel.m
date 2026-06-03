% Compile source code into MEX functions.  Requires C compiler.
% Run "mex -setup" first.
mex -v model_IHC.c complex.c  
mex -v model_Synapse.c complex.c 
mex -v model_IHC_BEZ2018.c complex.c  
mex -v model_Synapse_BEZ2018.c complex.c 
