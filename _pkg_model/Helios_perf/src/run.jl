using Helios_perf

funcs = [
    p_performance, 
    p_performance_timevec, 
    p_performance_prealloc, 
    () -> p_channel_scaling(20; moc_width=0.0), 
    () -> p_channel_scaling(20; moc_width=0.8), 
    p_performance_prealloc_benefit,
]
for func in funcs[[4, 5]]
    println("Running $func...")
    func()
end
