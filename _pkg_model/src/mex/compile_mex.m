%% compile_mex.m
% Author: Daniel R. Guest
% Date: 7/26/2023

% Delete existing intermediate files and compilex Mex file
delete *.obj;
delete *.o;
delete *.mex*;

% Branch based on platform
if ispc
	% Compile individual C files
	mex -c complex.c sfie.c adaptation.c model.c;  
	
	% Compile Mex wrapper
	mex sim_efferent_model_mex.c complex.obj sfie.obj adaptation.obj model.obj;
elseif ismac
	% Compile individual C files
	mex -c complex.c sfie.c adaptation.c model.c;  
	
	% Compile Mex wrapper
	mex sim_efferent_model_mex.c complex.o sfie.o adaptation.o model.o;
elseif isunix
	% Compile individual C files
	mex -c complex.c sfie.c adaptation.c model.c;  
	
	% Compile Mex wrapper
	mex sim_efferent_model_mex.c complex.o sfie.o adaptation.o model.o;
end