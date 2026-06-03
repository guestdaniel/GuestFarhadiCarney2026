%% Simple version
n_step = 250;
mems = zeros(n_step, 1);
for ii = 1:n_step
	[~, hsr, lsr, ic, gain] = sim_efferent_model(...
		zeros(50000, 1),...
		[1000.0]...
	);
	temp = memory;
	mems(ii) = temp.MemUsedMATLAB;
end
plot(mems*1e-6);
hold on;
xlabel('Sample');
ylabel('MATLAB memory usage (MB)');

mems = zeros(n_step, 1);
for ii = 1:n_step
	sim_ihc_zbc2014(zeros(50000, 1), 1000.0);
	temp = memory;
	mems(ii) = temp.MemUsedMATLAB;
end

plot(mems*1e-6);
hold off;
legend(["Efferent", "Old"]);

%% Parfor version
n_step = 5000;
mems = zeros(n_step, 1);
parfor ii = 1:n_step
	[~, hsr, lsr, ic, gain] = sim_efferent_model(...
		zeros(50000, 1),...
		[1000.0]...
	);
	temp = memory;
	mems(ii) = temp.MemUsedMATLAB;
end
plot(mems*1e-6);
hold on;
xlabel('Sample');
ylabel('MATLAB memory usage (MB)');

mems = zeros(n_step, 1);
parfor ii = 1:n_step
	sim_ihc_zbc2014(zeros(50000, 1), 1000.0);
	temp = memory;
	mems(ii) = temp.MemUsedMATLAB;
end

plot(mems*1e-6);
hold off;
legend(["Efferent", "Old"]);