export standardfields, standardclims, standardylims, standardcmap, standardlabel, standardavg, standardfreqs

"""
    standardfields(x)

Returns typical fields worth looking at when generating automated labels, names, etc.
"""
standardfields(::Type{AuditoryNerveZBC2014}) = [:τ_e_ic, :τ_i_ic, :d_i_ic, :S_ic, :A_ic]
standardfields(::Type{InferiorColliculusSFIEBE}) = [:τ_e_ic, :τ_i_ic, :d_i_ic, :S_ic, :A_ic]
standardfields(::Type{InferiorColliculusSFIEBS}) = [:species, :fiber_type, :audiogram]
function standardfields(x::String)
    if x == "MOC"
        return [:moc_weight_wdr, :moc_weight_ic]
    end
end

"""
    standardclims

Returns tuple containing suitable default color limits for type
"""
standardclims(::Type{AuditoryNerveZBC2014}) = (0.0, 1000.0)
standardclims(::Type{InferiorColliculusSFIEBE}) = (0.0, 200.0)
standardclims(::Type{InferiorColliculusSFIEBS}) = (0.0, 500.0)

"""
    standardylims

Returns tuple containing suitable default ylims for type
"""
standardylims(::Type{AuditoryNerveZBC2014}) = (0.0, 500.0)
standardylims(::Type{InferiorColliculusSFIEBE}) = (0.0, 100.0)
standardylims(::Type{InferiorColliculusSFIEBS}) = (0.0, 500.0)
function standardylims(x::String)
    @match x begin
        "hsr" => (0.0, 1000.0)
        "lsr" => (0.0, 150.0)
        "mocwdr" => (0.0, 70.0)
        "ic" => (0.0, 200.0)
        "mocic" => (0.05, 70.0)
        "gain" => (0.05, 1.1)
        "gainpostmix" => (0.01, 1.1)
        "spectrum" => (20.0, 70.0)
        "ihc" => (-1e-2, 1e-2)
    end
end

"""
    standardcmap

Returns colormap or identifier that can be used in place of colormaps in Makie
"""
standardcmap(x) = :acton

"""
    standardlabel

Returns suitable ylabel for type
"""
standardlabel(x) = "Firing rate (sp/s)"

"""
    standardavg

Returns a function that computes a reasonable average given input type
"""
function standardavg(x::String)
    @match x begin
        "ihc" => rms
        "hsr" => mean
        "lsr" => mean
        "ic" => mean
    end
end

"""
    standardfreqs

Returns standard frequencies for different analyses (e.g., audiometric)
"""
standardfreqs() = 0.125 .* 2 .^ (0:7)