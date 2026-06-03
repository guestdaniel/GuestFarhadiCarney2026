% Compare implementations
y1 = ffGn_2014(5000, 1/100e3, 0.9, 0, 10.0, 20.0);
y2 = ffGn_urear_2020b(5000, 1/100e3, 0.9, 0, 10.0, 20.0);
sum(y1 - y2);

% Compare driven/undriven rates
figure;
fs = 100e3;                                      % sample rate (Hz)
dur = 0.1;                                       % duration (seconds)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)
x = 20e-6 * 10^(50.0/20.0) * sin(2*pi * 1000.0 * t)*sqrt(2);
n_rep = 500;
spont_2014 = zeros(n_rep, 1);
spont_2018 = zeros(n_rep, 1);
drive_2014 = zeros(n_rep, 1);
drive_2018 = zeros(n_rep, 1);
parfor ii = 1:n_rep
	fprintf('Running rep %d\n', ii);
	[~, hsr, ~, ~, ~] = sim_efferent_model_ffGn_2014(zeros(size(x)), [1000.0]);
	spont_2014(ii) = mean(hsr);
	[~, hsr, ~, ~, ~] = sim_efferent_model_ffGn_urear_2020b(zeros(size(x)), [1000.0]);
	spont_2018(ii) = mean(hsr);
	[~, hsr, ~, ~, ~] = sim_efferent_model_ffGn_2014(x, [1000.0]);
	drive_2014(ii) = mean(hsr);
	[~, hsr, ~, ~, ~] = sim_efferent_model_ffGn_urear_2020b(x, [1000.0]);
	drive_2018(ii) = mean(hsr);
end
figure;
hold on;
ksdensity(spont_2014);
ksdensity(spont_2018);
ksdensity(drive_2014);
ksdensity(drive_2018);
hold off;