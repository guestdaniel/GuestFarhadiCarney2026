%% Demo 1: Same input, no fGn, comparison between different PLA modes
% 0 -> old approximation
% 1 -> original
% 2 -> new approximation
figure;
fs = 100e3;                                      % sample rate (Hz)
dur = 0.5;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
x = 20e-6 * 10^(50.0/20.0) * sin(2*pi * 1000.0 * t)*sqrt(2);
ihc = model_IHC(x, 1000.0, 1, 1/100e3, length(x)/100e3, 1.0, 1.0, 2);
tiledlayout(1, 3);
for ii = 0:2
	nexttile;
	[r, ~, ~] = model_Synapse_lightspeed(ihc, 1000.0, 1, 1/100e3, 3, 0, ii, 100e3);
	plot(r)
end

%% Demo 2: Performance at various fs as function of stimulus duration, powerlaw mode
figure;
tiledlayout(1, 3);
durs = 0.01 .* 2 .^ (0:0.5:7);
fss = [10e3, 20e3, 100e3];
for idx_fs = 1:length(fss)
	nexttile;
	for idx_mode = 1:3
		data{idx_mode} = zeros(length(durs), 1);
		for idx_dur = 1:length(durs)
			x = zeros(1, round(durs(idx_dur)*100e3));
			tic;
			[r, ~, ~] = model_Synapse_lightspeed(x, 1000.0, 1, 1/100e3, 3, 0, idx_mode-1, fss(idx_fs));
			data{idx_mode}(idx_dur) = toc;
		end
	end
	hold on;
	for ii = 1:3
		plot(durs, data{ii});
	end
	set(gca, 'Xscale', 'log');
	set(gca, 'yscale', 'log');
	legend(["Old approx", "Orig", "New approx"]);
	hold off;
	title(num2str(fss(idx_fs)));
end


%% Demo 3: Performance at various durations as function of fs, powerlaw mode
figure;
tiledlayout(1, 4);
durs = [0.1, 0.2, 0.4, 0.8];
fss = [10e3, 20e3, 50e3, 100e3];
for idx_dur = 1:length(durs)
	nexttile;
	for idx_mode = 1:3
		data{idx_mode} = zeros(length(fss), 1);
		for idx_fs = 1:length(fss)
			x = zeros(1, round(durs(idx_dur)*100e3));
			tic;
			[r, ~, ~] = model_Synapse_lightspeed(x, 1000.0, 1, 1/100e3, 3, 0, idx_mode-1, fss(idx_fs));
			data{idx_mode}(idx_fs) = toc;
		end
	end
	hold on;
	for ii = 1:3
		plot(fss, data{ii});
	end
	set(gca, 'yscale', 'log');
	if idx_dur == length(durs)
		legend(["Old approx", "Orig", "New approx"]);
	end
	xlabel("Sampling rate (Hz)");
	ylabel("Runtime (s)");
	hold off;
	title(num2str(durs(idx_dur)) + " s stimulus");
end