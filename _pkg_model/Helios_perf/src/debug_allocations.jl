using Helios

len_total = 10_000
n_chan = 11
controlout = [zeros(len_total) for _ in 1:n_chan]
c1out = [zeros(len_total) for _ in 1:n_chan]
c2out = [zeros(len_total) for _ in 1:n_chan]
ihcout = [zeros(len_total) for _ in 1:n_chan]
expout_hsr = [zeros(len_total) for _ in 1:n_chan]
sout1_hsr = [zeros(len_total) for _ in 1:n_chan]
sout2_hsr = [zeros(len_total) for _ in 1:n_chan]
synout_hsr = [zeros(len_total) for _ in 1:n_chan]
expout_lsr = [zeros(len_total) for _ in 1:n_chan]
sout1_lsr = [zeros(len_total) for _ in 1:n_chan]
sout2_lsr = [zeros(len_total) for _ in 1:n_chan]
synout_lsr = [zeros(len_total) for _ in 1:n_chan]
hsrout = [zeros(len_total) for _ in 1:n_chan]
lsrout = [zeros(len_total) for _ in 1:n_chan]
cnout = [zeros(len_total) for _ in 1:n_chan]
icout = [zeros(len_total) for _ in 1:n_chan]
mocwdr = [zeros(len_total) for _ in 1:n_chan]
mocic = [zeros(len_total) for _ in 1:n_chan]
gain = [ones(len_total) for _ in 1:n_chan]
gainpostmix = [ones(len_total) for _ in 1:n_chan]
ffGn_hsr, ffGn_lsr = Helios.prepare_ffGn(len_total, false, n_chan)
y = randn(len_total)
cfs = collect(200.0:200.0:2200.0)

@allocated sim_gfc2023!(ffGn_hsr, ffGn_lsr, controlout, c1out, c2out, ihcout, expout_hsr, sout1_hsr, sout2_hsr, synout_hsr, expout_lsr, sout1_lsr, sout2_lsr, synout_lsr, hsrout, lsrout, cnout, icout, mocwdr, mocic, gain, gainpostmix, y, cfs; 
    dur_pad_left=0.0,
    dur_pad_right=0.0,
)
