%% demo.m
% This script is a series of examples that aim to demonstrate how to use
% the Carney-lab efferent model code. Each example includes a brief
% description and key ideas to take away.

%% Example #1: Complex tone, 6-10F0, 50 dB SPL per component, sine phase
% Here, we synthesize a five-component complex tone with frequencies at
% 6F0-10F0 for an F0 of 200 Hz at 50 dB SPL per component and then simulate 
% and visualize the efferent-model response to this stimulus. 
%
% **Key ideas:** scaling acoustic waveforms, running the model

% Set parameters
fs = 100e3;                                      % sample rate (Hz)
fs_down = 10e3;                                  % sample rate for plot (Hz)
dur = 0.5;                                       % duration (seconds)
f0 = 200.0;                                      % F0 (Hz)
t = 0.0:(1/fs):(dur - 1/fs);                     % sample times (s)  
n_cf = 41;                                       % number of channels (#)
cf = exp(linspace(log(0.5e3), log(3e3), n_cf));  % CFs (Hz)

% Construct stimulus
x = zeros(size(t));
for harm_no = 6:10
	% Here, we synthesize one harmonic, scale it to 50 dB SPL, and add it
	% to the stimulus
	x = x + 20e-6 * 10^(50.0/20.0) * sin(2*pi* harm_no*f0 * t)/0.7071;
end
x = x .* tukeywin(length(x), 0.05)';

% Simulate efferent-model response
[~, hsr, lsr, gain] = sim_efferent_model(x, cf, display_info=true, noiseType=-1);
responses = {hsr, lsr, gain};  % store all output matrices in cell array

% Resample responses down to lower sampling rate for plotting
t_resampled = 0.0:(1/fs_down):(dur - 1/fs_down);
for ii = 1:length(responses)
	responses{ii} = resample(responses{ii}, fs_down, fs, Dimension=2);
end

% Define labels/limits for plot
labels = ["HSR", "LSR", "Gain"];
limits = {[0.0, 600.0], [0.0, 200.0], [0.0, 1.0]};

% Plot as colorplots/neurograms
figure;
tiledlayout(length(responses), 1);
for ii = 1:length(responses)
	nexttile;
	imagesc(t_resampled, 1:length(cf), responses{ii});
	set(gca, 'ydir', 'normal');
	caxis(limits{ii});
	yticks(1:10:n_cf);
	yticklabels(round(cf(1:10:n_cf)));
	xlabel('Time (s)');
	ylabel('CF (Hz)');
	title(labels(ii));
	xlim([0.1 0.25]);  % zoom into only 100-250 ms range
	colorbar;
end
set(gcf, "Units", "normalized");
set(gcf, "Position", [0.2 0.3 0.5 0.6])

%% Example #2: Variability in  spontaneous rate
% In the Zilany/Bruce/Carney et al. models, model versions from 2009 and
% onward include some amount of fractional Gaussian noise added to the
% synapse output inside the power-law adaptation loop. This noise has
% several effects, one of which is to cause the spontaneous rate of the
% model to "drift" over time with a long-range temporal dependence.
%
% In some cases, this behavior is undesireable, such as when parameters are
% being optimized using a numerical optimization procedure that assumes a
% deterministic function. Thus, two options are provided for control.
% First, by passing `noiseType=0`, the "fresh" noise is replaced with a
% frozen sample of noise that is reused on subsequent calls to the model
% --- this makes the model completely deterministic. Second, by passing
% `noiseType=-1`, the noise is eliminated entirely from the model. In other
% words, the "sample" of noise used is a sample containing zeros at time
% points. This can have unintended side effects, so this is not recommended
% unless you know what you are doing (but is shown here as a learning
% opportunity).
% 
% Below, example spontaneous responses are shown for each noiseType
% setting. For clarity, different waveforms are slightly separated by
% vertical offets so that identical responses do not lie on top of each
% other. To simulate spontaneous responses, we simply simulate responses to
% silence for some 500 ms. 
%
% **Key ideas:** Spontaneous rate, fractional Gaussian noise (fGN)

% Create figure
figure;
noiseTypes = [-1, 0, 1];
labels = [
	"No fGn (noiseType == -1)", ...
	"Frozen fGn (noiseType = 0)", ...
	"Fresh fGn (noiseType = 1)" ...
];

% Loop through possibilities, generate responses and plot 
for ii = 1:3
	% Set up plot and titles
	subplot(1, 3, ii);
	subtitle(labels(ii));
	xlabel('Time (s)');
	ylabel('Firing rate (sp/s)');
	ylim([0.0, 750.0]);

	% Run simulations 10 times
	hold on;
	for jj = 1:10
		% Run model
		[~, hsr, ~, ~] = sim_efferent_model( ...
			zeros(1, 50000), ...
			[1000.0], ...
			noiseType=noiseTypes(ii) ...
		);

		% Plot result
		t = 0.0:(1/100e3):(0.5-1/100e3);
		plot(t, hsr + jj);
	end
	hold off;
end
set(gcf, "Units", "normalized");
set(gcf, "Position", [0.2 0.3 0.5 0.2])

%% Example #3: Comparison of true vs approximate power law for many pure tones
% In the AN model, computation of power-law adaptation is one of the most
% computationally demanding stages of the model. One way to get around this
% is to replace "true" power-law adaptation with an approximation composed
% of many parallel exponential adaptation processes with time constants
% varying over a wide range (to capture the behavior of the power law at
% different time scales). These exponential adaptation processes can be
% efficiently implemented as IIR filters that depend only on the previous
% sample, rather than all previous samples (as for the "true" power-law
% adaptation).
% 
% This demo simulates responses to pure tones under the "true" power-law
% adaptation or under its approximate implementation and compares them
% side-by-side. Note that the demo relies on `parfor` to speed things up
% (in line 252). This can be replaced with `for` if you do not have the
% parallel toolbox or are running into issues. Note that, because the
% approximation is highly accurate, these responses will look nearly
% identical (slight discrepancies can be noted following pure-tone
% offsets).
%
% **Key ideas:* Power-law adaptation, power-law approximation

% Set parameters
fs = 100e3;                                      % sample rate (Hz)
dur = 0.2;                                       % duration (s)
dur_post = 0.1;                                  % duration of post-stimulus simulation time (s)
t = 0.0:(1/fs):(dur+dur_post - 1/fs);            % sample times (s)
freqs = [500.0, 1000.0, 2000.0, 4000.0, 8000.0]; % test frequencies (Hz)
time_windows = {[0.05, 0.06], [0.1, 0.3]};  % time windows for analysis (s)

% Pre-allocate storage
true = zeros(length(t), length(freqs));
approx = zeros(length(t), length(freqs));

% Loop over stimuli and do calculations (in parallel using parfor)
parfor idx_freq = 1:length(freqs)
	% Synthesize stimulus for this frequency (50 dB SPL pure tone)
	stim = 20e-6 * 10^(50.0/20.0) * ...
		sin(2*pi * freqs(idx_freq) * (0.0:(1/fs):(dur - 1/fs))) * sqrt(2);
	stim = [stim zeros(1, round(dur_post*fs))];

	% Run efferent model with true and approximate power-law adaptation
	[~, true(:, idx_freq), ~, ~] = ...
		sim_efferent_model(stim, freqs(idx_freq), powerlaw_mode=1, noiseType=0);
	[~, approx(:, idx_freq), ~, ~] = ...
		sim_efferent_model(stim, freqs(idx_freq), powerlaw_mode=2, noiseType=0);	
end

% Create figure
figure;

% Loop through frequencies
for idx_freq = 1:length(freqs)
	% Loop through different time scales
	for idx_tw = 1:length(time_windows)
		% Create subplot
		subplot( ...
			length(freqs), ...
			length(time_windows), ...
			length(time_windows)*(idx_freq-1) + idx_tw ...
		);

		% Plot true and approximate results on top of each other
		plot(t, true(:, idx_freq)); hold on;
		plot(t, approx(:, idx_freq)); hold off;

		% Add legend to upper left
		if idx_freq == 1 && idx_tw == 1
			legend(["True power-law adaptation", "Approximate power-law adaptation"]);
		end

		% Set limits, titles, and labels
		xlim(time_windows{idx_tw});
		xlabel("Time (s)");
		ylabel("Firing rate (sp/s)");
		ylim([0.0, 1000.0]);
		title(sprintf( ...
			"Freq = %4.0f Hz, time window = [%4.3f, %4.3f] s", ...
			freqs(idx_freq), ...
			time_windows{idx_tw}(1), ...
			time_windows{idx_tw}(2) ...
		));
	end
end
set(gcf, "Units", "normalized");
set(gcf, "Position", [0.2 0.3 0.5 0.5])

%% Example #4: Performance benefits of approximate power-law implementation
% See above for more details about power-law adaptation approximation. 
%
% Here, we test the approximation and compare its performance to the true
% power-law adaptation implementation. Feel free to adjust the durations
% below to determine the performance gains you can expect for different
% stimuli, although beware that durations beyond ~1s start to require
% prohibitively long compute times for the true power law (so it may take 
% quite some time to generate). The first graph simply depicts compute
% time, while the other graph shows the "performance gain" (expressed as a 
% ratio between runtimes, with larger values indictaing greater
% improvements in speed) from switching to the approximate power-law 
% adaptatation scheme.
%
% Note that time estimates for shorter stimuli can be somewhat unreliable
% --- if the graphs look nonmonotonic or have outliers, do not trust those
% points, since they can be unduly influenced by brief delays or
% interruptions in compute due to system processes or other factors. Also
% since we estimate compute time multiple times to eliminate some of this
% randomness, the plot can take quite some time to generate. 
%
% **Key ideas:* Model performance, power-law approximation, runtime

% Set parameters
dur_min = 5e-2;
dur_max = 5e-1;
durs = exp(linspace(log(dur_min), log(dur_max), 10));  % stimulus durations (s)
n_rep = 5;                                             % how many repeats to do

% Pre-allocate storage
durs_true = zeros(length(durs), n_rep);
durs_approx = zeros(length(durs), n_rep);

% Time compute time for each stimulus duration
for ii = 1:length(durs)
	% Run model with real power-law adaptation
	for jj = 1:n_rep
		tic;
		sim_efferent_model(zeros(1, round(durs(ii)*100e3)), [1000.0], powerlaw_mode=1);
		durs_true(ii, jj) = toc;
	end

	% Run model with approximate power-law adaptation
	for jj = 1:n_rep
		tic;
		sim_efferent_model(zeros(1, round(durs(ii)*100e3)), [1000.0], powerlaw_mode=2);
		durs_approx(ii, jj) = toc;
	end
end

% Plot compute time
figure;
subplot(1, 2, 1);
errorbar(durs, mean(durs_true, 2), 1.96 * std(durs_true, 0, 2)/sqrt(n_rep), 'k'); hold on;
grid on;
errorbar(durs, mean(durs_approx, 2), 1.96 * std(durs_true, 0, 2)/sqrt(n_rep), 'r');
plot(10 .^ (-3.0:0.01:3.0), 10 .^ (-3.0:0.01:3.0), 'k'); hold off;  % plot unity line
set(gca, "xscale", "log");
set(gca, "yscale", "log");
ylim([1e-2, 1e2]);
xlim([1e-2, 1e1]);
legend(["True power law adaptation", "Approximate power law adaptation"]);
xlabel("Stimulus duration (s)");
ylabel("Compute time (s)");

% Plot speed up
subplot(1, 2, 2);
plot(durs,  mean(durs_true, 2) ./  mean(durs_approx, 2), 'k'); grid on;
set(gca, "xscale", "log");
set(gca, "yscale", "log");
ylim([1e0, 1e3]);
xlim([2e-2, 2e0]);
xlabel("Stimulus duration (s)");
ylabel("Performance gain (true/approx compute time)");

set(gcf, "Units", "normalized");
set(gcf, "Position", [0.2 0.3 0.5 0.5])

%% Example #5: COHC and CIHC
% Every channel in the simulation has an associated COHC and CIHC value
% that governs, respectively, the gain provided by OHCs and the amplitude
% of IHC responses. Altering these values can simulate different types of
% hearing loss. Here, we show a simulation of different CFs with different
% COHCs values in response to a pure tone presented at 4 kHz. Guidance is
% available in the literature for how to configure these parameters
% appropriately based on audiograms. (note: noiseType is fixed at -1 to
% make differences between NH and HI easier to see in a single simulation).
%
% **Key ideas:* Hearing loss, COHC, CIHC

% Choose parameters
freq = 4e3;                    % frequency (Hz)
level = 45.0;                  % level (dB SPL)
dur = 0.1;                     % duration (s)
cfs = exp(linspace(log(1e3), log(8e3), 31));    % CFs for each channel
cohcs_nh = ones(size(cfs));                  % COHC values for each channel, normal-hearing
cohcs_hi = linspace(0.7, 0.0, length(cfs));  % COHC values for each channel, hearing-impaired
									         % HL assumed to be gently sloping to 8 kHz
% Synthesize stimulus
t = 0.0:(1/fs):(dur-1/fs);
x = sin(2*pi*t * freq);
x = 10^(level/20)*20e-6 * x/rms(x);
x = x .* tukeywin(length(x), 0.1)';

% Simulate responses
[~, hsr_nh, ~, ~] = sim_efferent_model(x, cfs, cohc=cohcs_nh, noiseType=-1);
[~, hsr_hi, ~, ~] = sim_efferent_model(x, cfs, cohc=cohcs_hi, noiseType=-1);

% Plot average rate at each CF
figure;
plot(cfs, mean(hsr_nh, 2)); hold on;
set(gca, "Xscale", "log");
plot(cfs, mean(hsr_hi, 2)); hold off;
ylim([0.0, 500.0]);
xticks([1e3, 2e3, 4e3, 8e3]);
xlim([1e3/sqrt(2), 8e3*sqrt(2)]);
xlabel("CF (Hz)");
ylabel("Firing rate (sp/s)");
legend(["Normal hearing", "Hearing impaired"]);

%% Example #6: Single-channel versus multi-channel simulations of gain control
% The present model of MOC efferent gain control has a parameter that
% controls how "wide" the WDR-driven efferent projections are,
% `moc_width`. When set to a value of zero, each channel's WDR-driven 
% gain factor is determined only by WDR-pathway responses in that channel. 
% When set to a value greater than zero, WDR-driven gain factors for a
% range of channels are geometrically averaged to produce a final
% WDR-driven gain factor. Note that this is a
% highly experimental and novel feature of the model and and such is
% constantly being updated. 
% 
% The simulations below explore the consequences of this width parameter on
% responses to (and resulting gain control evoked by) a filtered noise
% stimulus at various sound levels.  To make the differences easier to
% observe, we disable IC-driven gain control and set the strength of
% WDR-driven gain control to be quite high. The expectation is that
% single-channel simulations and multi-channel simulations with narrow
% widths should resemble one another, whereas once the width is wider than
% the response area of the narrowish-band stimulus, gain control should
% lessen (because the gain factor for the plotted channel will be an
% average of channels with strong gain control signals [those excited by
% the noise] and channels with weaker gain control signals [those not
% excited by the noise]).
%
% **Key ideas:* IC vs WDR gain control, MOC width

% Set parameters
fs = 100e3;
cf = 1e3;                                                    % middle CF 
cfs = exp(linspace(log(1e3 * 2^(-1)), log(1e3 * 2^1), 21));  % range of CFs around middle CF
dur = 0.5;                                   % duration (dur)
Nos = [-10.0:10.0:20.0];                     % spectrum levels of noise (dB SPL)
[b, a] = butter(4, [1e3*2^(-1/6), 1e3*2^(1/6)] / (fs/2));  % bandpass filter for noise
t = 0.0:(1/fs):(dur-1/fs);                   % time axis (s)

% Create figure
figure;
tiledlayout(1, 4);

% First, we plot results for single-channel simulations
nexttile; hold(gca, "on");
for idx_level = 1:length(Nos)
	level = Nos(idx_level) + 10*log10(fs/2);
	x = randn(1, round(dur*fs));
	x = 10^(level/20)*20e-6 * x/rms(x);
	x = filter(b, a, x);
	[~, hsr, ~, gain] = sim_efferent_model(x, cf, moc_weight=10.0, noiseType=-1);
	plot(t, gain);
end
hold(gca, "off");
ylim([0.0, 1.1]);
xlabel("Time (s)");
ylabel("Gain factor");
legend(string(Nos));
title("Single-channel")

% Next, we plot results for multi-channel simulations with various
% moc_width values, looking only at time-varying gain factor from the
% CF in the middle of the simulated range
for width = [0 0.5 1.0]
	nexttile; hold(gca, "on");
	for idx_level = 1:length(Nos)
		level = Nos(idx_level) + 10*log10(fs/2);
		x = randn(1, round(dur*fs));
		x = 10^(level/20)*20e-6 * x/rms(x);
		x = filter(b, a, x);
		[~, hsr, ~, gain] = sim_efferent_model(x, cfs, noiseType=-1, moc_weight=10.0, moc_width=width);
		plot(t, gain(11, :));  % target CF is 11th channel
	end
	hold(gca, "off");
	ylim([0.0, 1.1]);
	xlabel("Time (s)");
	ylabel("Gain factor");
	legend(string(Nos));
	title(sprintf("Multichannel (width = %0.1f oct)", width))
end

%% Example 7: Rate-level function
% [[fill in description of RLF]]
% Set parameters
fs = 100e3;
cf = 10e3;
levels = 0.0:5.0:90.0;
dur = 0.03;

% Simulate response to each RLF stimulus
rlf = zeros(size(levels));
for idx_level = 1:length(levels)
	% Do stim
	stim = 20e-6 * 10^(levels(idx_level)/20.0) * sin(2*pi * cf * (0.0:(1/fs):(dur - 1/fs))) * sqrt(2);
	stim = stim .* tukeywin(length(stim), 0.1)';

	% Do sim
	[~, hsr, ~, ~] = sim_efferent_model(stim, cf, display_info=true, noiseType=-1);
	rlf(idx_level) = mean(hsr);
end

% Plot
figure;
plot(levels, rlf);
xlim([0 90]);
ylim([0 500]);