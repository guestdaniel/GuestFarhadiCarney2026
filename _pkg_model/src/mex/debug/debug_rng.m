% Start by completely re-initializing RNG
rng('default');

% Collect "old" numbers by forcibly setting seed with deprecated `randn`
% call
randn('seed', 37);  % set seed
old = randn(25, 1);
rng('default');     % reset seed
matched = randn(RandStream('v4', 'Seed', 37), 25, 1);
new = randn(RandStream('mt19937ar', 'Seed', 37), 25, 1);

% Compare frozen and fresh via streams
frozen = {};
fresh = {};
for repeat = 1:5
	frozen{repeat} = randn(RandStream('v4', 'Seed', 37), 25, 1);
	fresh{repeat} = randn(RandStream.getGlobalStream(), 25, 1);
end

% Examine "fresh" vs "frozen" waveforms
figure;
subplot(1, 2, 1);
hold on;
for ii = 1:10
	plot(ffGn(100e3, 1/args.fs, 0.9, 0, 1.0, 100.0));
end
hold off;
subplot(1, 2, 2);
hold on;
for ii = 1:10
	plot(ffGn(100e3, 1/args.fs, 0.9, 1, 1.0, 100.0));
end
hold off;