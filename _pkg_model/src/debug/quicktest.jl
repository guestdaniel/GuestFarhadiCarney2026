using Helios
using UnicodePlots
using AuditorySignalUtils

stim = scale_dbspl(pure_tone(1000.0, 0.0, 0.1, 100e3), 50.0)
resp = sim_gfc2023_dict(stim, 1000.0; dur_pad_left=0.0, clip_left=false)

lineplot(resp["hsr"])