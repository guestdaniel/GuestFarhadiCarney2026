%% measure_memory_leak.m
% This script is intended to repeatedly run the efferent model, either by
% single for loop or by parfor, and to look for evidence of a memory leak.
% Current memory usage will be recorded and logged during runs so that 
% memory usage can be reconstructed afterward even if a crash occurs.

% Set up logger
fn = fullfile("logs", string(datetime("today")) + ".txt");

% Loop over series of empty efferent model runs
mem = zeros(500, 1);
for idx = 1:500
	sim_efferent_model(zeros(50000, 1), 1000.0);
	temp = memory;
	mem(idx) = temp.MemUsedMATLAB;
end
plot(mem);