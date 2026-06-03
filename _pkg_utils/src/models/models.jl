export Model, Null, compute, extract, samprate, getstage, getcoi, issinglechannel, ismultichannel, islogout, getcf, getcentercf, singlechannel

"""
    Model

Abstract type for a mapping from input to ouptut

Below is a description of the informal method interface developed around Model:
- `compute(::Model, ::AbstractStimulus)` should compute the output of a model given a stimulus
- `compute(::Model, ::Vector)` should compute the output of a model given a sound-pressure
  waveform
"""
abstract type Model <: Component end
@with_kw struct Null <: Model 
    cf::Vector{Float64}=[1000.0]
end

compute(m::Model, s::AbstractStimulus) = compute(m, synthesize(s))
samprate(m::Model) = m.fs
getcf(m::Model) = m.n_chan == 1 ? m.cf[m.coi[1]] : m.cf
getcentercf(m::Model) = ismultichannel(m) ? getcf(m) : middle(getcf(m))
getstage(m::Model) = m.stage
getcoi(m::Model) = m.coi
issinglechannel(m::Model) = length(getcoi(m)) == 1
ismultichannel(m::Model) = length(getcoi(m)) > 1
islogout(m::Model) = false
function copy_and_alter(s, field::Symbol, value)
    T = typeof(s)
    fields = fieldnames(T)
    
    # Create a tuple of values to construct the new struct
    new_values = map(f -> f === field ? value : getproperty(s, f), fields)
    
    # Use the splat operator `...` to pass the values to the constructor
    return T(new_values...)
end
singlechannel(m::Model) = copy_and_alter(m, :coi, [middle(1:length(m.cf))])

# Handle logic of extract response from model output
function extract(m::Model, r) 
# Branch based on opossible cases
    if m.n_chan == 1
        # If we have a single-channel model, return single-channel response as Vector{Float64}
        return r[1]
    elseif length(m.coi) == 1
        # If we have a multichannel model but only one channel of interest, return primary 
        # channel response at Vector{Float64}
        return r[m.coi[1]]
    else
        # Otherwise return all channels as vector of vectors
        r[m.coi]
    end
end