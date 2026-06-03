export Simulation, Response, simulate, analyze, modelobj, samprate, stimobj, status, getcf, hascache, axis, spont

"""
    Simulation <: Component

Abstract type for encoding how to combine Components to achieve a result 

Simulations are composed of a handful of Components or iterable containers of Components.
Simulations must provide an implementation of `simulate(::Simulation)`, which describes how
to combine Components to produce a result or output (most typically, a numeric value or some
Array of numeric values, such as a vector of firing rates).

# Methods
- `simulate(::Simulation)`: produce the output of the simulation
- `analyze(::Simulation, args...)`: analyze a single simulation result
"""
abstract type Simulation <: Component end

function simulate(::Simulation) "Hitting this when you're not supposed to!" end
simulate(sim::Simulation, ::Config) = simulate(sim)

function analyze(::Simulation) end
function modelobj(x::Simulation) x.model end
function stimobj(x::Simulation) x.stim end
function samprate(x::Simulation) samprate(modelobj(x)) end
function getcf(x::Simulation) getcf(modelobj(x)) end
function hascache(x::Simulation, config::Config) isfile(config, simulate, x) end
function axis(x::Simulation) end

"""
    Response{S, M}

Simulation consisting of generating a response from a model::M for a stimulus::S
"""
@with_kw struct Response{S, M} <: Simulation where {S <: AbstractStimulus, M <: Model}
    stimulus::S
    model::M
end

function simulate(r::Response)
    compute(r.model, r.stimulus)
end

timeaxis(r::Response) = collect(0.0:(1/r.model.fs):(r.stimulus.dur - 1/r.model.fs))
cfaxis(r::Response) = r.model.cf

"""
    spont(model::Model)

Return simulated estimate of spont rate based on 200 ms of silence
"""
function spont(model::Model; config=Config())
    @memo config mean(compute(model, zeros(20000)))
end