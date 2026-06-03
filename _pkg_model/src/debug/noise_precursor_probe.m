function stim = noise_precursor_probe(args)
	arguments
		% Set overall params
		args.cf = 8e3;            % CF (Hz)
		args.bw_pre = 1/3;        % bandwidth (octaves)
		args.level_pre = 20.0;    % overall level of precursor (dB SPL)
		args.level_probe = 20.0;  % overall level of probe (dB SPL)
		args.dur_pre = 0.3;       % duration of precursor (s)
		args.dur_ppi = 0.1;       % duration precursor-probe interval (s)
		args.dur_probe = 0.01;    % duration probe (s)
		args.fs = 100e3;          % sampling rate (Hz)
	end
	len_pre = round(args.dur_pre*args.fs);
	
	% Create precursor
	precursor = randn(len_pre, 1);
	[b, a] = butter(4, [args.cf * 2 ^ (-args.bw_pre/2), args.cf * 2 ^ (args.bw_pre/2)] ./ (args.fs/2));
	precursor = filter(b, a, precursor);
	precursor = raised_cosine_ramp(precursor', 0.001, args.fs)';
	precursor = 20e-6 * 10^(args.level_pre/20.0) * precursor/rms(precursor);
	
	% Create probe
	t = 0.0:(1/args.fs):(args.dur_probe - 1/args.fs);
	probe = sin(2 * pi * args.cf * t');
	probe = raised_cosine_ramp(probe', 0.001, args.fs)';
	probe = 20e-6 * 10^(args.level_probe/20.0) * probe/rms(probe);
	
	% Create total stimulus
	stim = [precursor; zeros(round(args.dur_ppi*args.fs), 1); probe];
end