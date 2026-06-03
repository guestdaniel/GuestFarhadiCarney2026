function quicklook_efferent(stim, cf, fs, args)
	arguments
		stim
		cf
		fs=100e3
		args.moc_beta=0.02
		args.moc_offset=2.0
		args.powerlaw_mode=2
	end

	% Pass stimulus through efferent model
	[ihc, hsr, lsr, gain] = sim_efferent_model( ...
		[stim; zeros(5000, 1)], ...
		cf, ...
		noiseType=-1, ...
		moc_offset=args.moc_offset, ...
		powerlaw_mode=args.powerlaw_mode ...
	);

	% Plot everything
	figure;
	tiledlayout(4, 1, "TileSpacing", "compact", "Padding", "compact");
	responses = {ihc, hsr, lsr, gain};
	labels = ["IHC", "HSR", "LSR", "Gain"];
	for idx_resp = 1:length(responses)
		resp = responses{idx_resp};
		t = 0.0:(1/fs):(length(resp)/fs - 1/fs);
		nexttile;
		plot(t, resp);
		title(labels(idx_resp));
		if idx_resp < length(responses)
			xticklabels([]);
		end
		if idx_resp == length(responses)
			ylim([0.0, 1.1]);
		else
			ylim([0.0, max(resp) * 1.2])
		end
	end
end

