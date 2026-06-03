clear all;

% mex -v model_Synapse_first_actual_imp.c complex.c 

 mex model_IHC_sample_actual_imp.c complex.c 

 mex model_IHC_sample_lsr_actual_imp.c complex.c 
 mex model_IHC_first.c complex.c
% mex -v model_IHC_sample.c complex.c  
% % clear all;
% mex -v model_IHC_sample_lsr.c complex.c  
% mex -v model_Synapse_first.c complex.c