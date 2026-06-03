% Use brute-force to find CF that yields maximum response to 50-dB 1-kHz
% pure tone
f = 250.0:25.0:2000.0;
rates = zeros(size(f));
for ii = 1:length(f)
	rates(ii) = meanrate_at_cf(f(ii), -1);
end
[~, idx_max] = max(rates);

fprintf('Best CF via brute force, noiseType==-1, is %f\n', f(idx_max));

% Use fmincon to find the CF that yields the maximum response to 50-dB
% 1-kHz pure tone
[x, ~, ~, ~, ~, ~, ~] = ...
	fmincon(@(x) -meanrate_at_cf(x, -1), [500.0], [], [], [], [], [250.0], [2000.0]);
fprintf('Best CF via fmincon, noiseType==-1, is %f\n', x);
[x, ~, ~, ~, ~, ~, ~] = ...
	fmincon(@(x) -meanrate_at_cf(x, 0), [500.0], [], [], [], [], [250.0], [2000.0]);
fprintf('Best CF via fmincon, noiseType==0, is %f\n', x);
for idx_attempt = 1:5
	[x, ~, ~, ~, ~, ~, ~] = ...
		fmincon(@(x) -meanrate_at_cf(x, 1), [500.0], [], [], [], [], [250.0], [2000.0]);
	fprintf('Best CF via fmincon, noiseType==1 attempt #%d, is %f\n', idx_attempt, x);
end

% Function to compute mean HSR response to 1-kHz pure tone at CF of x 
function out = meanrate_at_cf(x, noiseType)
	% Create stimulus
	fs = 100e3;                                      % sample rate (Hz)
	dur = 0.1;                                       % duration (seconds)
	t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
	stim = 20e-6 * 10^(50.0/20.0) * sin(2*pi * 1000.0 * t)*sqrt(2);

	% Run model
	[~, hsr, ~, ~, ~] = sim_efferent_model(stim, [x], noiseType=noiseType);
	out = mean(hsr);
end