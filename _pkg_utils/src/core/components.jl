export Component
export id, cachepath, validate, hashparts

"""
    Component

Fundamental ingredient of auditory computational modeling (e.g., Model, Stimulus)

Component is an abstract type that provides various methods useful for computational 
modeling work, such as generating unique ID values for unique combinations of field values.

# Methods:
- `id(x)`: Returns string that uniquely maps to `x`
- `id(x, y)` Returns a string that uniquely maps to the combination of `x` and `y`
- `cachepath(x...)`: Returns path of the cache file for this combination of components
- `compute(x)` or `compute(x, args...)`: runs the core computations of the Component to
  yield an output
"""
abstract type Component end

"""
    hash(comp::Component)

Compute hash for a component
"""
function Base.hash(comp::C) where C <: Component
    # Start hashid with type
    hashids = [hashcomptype(comp)]

    # Get hashes for field names
    names = hashfieldnames(comp)

    # Get hashes for field values
    values = hashfieldvalues(comp)

    # Combine with hash for component type
    hashids = vcat(hashids, (names .+ values)...)

    # Reduce all hashes together via summation
    reduce(+, hashids)
end

function Base.hash(comps::Vector{C}) where C <: Component
    reduce(+, Base.hash.(comps))
end

function hashcomptype(::C) where C <: Component
    hash(string(C))
end

function hashfieldnames(::C) where C <: Component
    hash.(fieldnames(C))
end

function hashfieldvalues(comp::C) where C <: Component
    # Map over the fields of comp 
    map(fieldnames(C)) do fn
        # Get value of field fn
        val = getfield(comp, fn)

        # Branch based on type of fn
        if typeof(val) <: Function
            # If the val is a subtype of Function, we merely hash its name
            hash(string(val))
        else
            # Otherwise, we hash it normally
            hash(val)
        end
    end
end

function hashparts(comp::C) where C <: Component
    names = hashfieldnames(comp)
    values = hashfieldvalues(comp)
    println("typehash: $(string(hashcomptype(comp); base=16))")
    map(zip(fieldnames(C), names, values)) do (name, namehash, valuehash)
        println("$name: $(string(namehash; base=16)) + $(string(valuehash; base=16)) = $(string(namehash + valuehash; base=16))")
    end
    println("Final hash: $(string(Base.hash(comp); base=16))")
end

Base.isequal(x::Component, y::Component) = (Base.hash(x) == Base.hash(y))
function Base.:(==)(x::Component, y::Component)
    Base.isequal(x, y)
end

"""
    id(comp::Component[; accesses=nothing, connector="_"])

Return string uniquely mapping to field values of `comp`

Uses `DrWatson.savename` to generate a string that uniquely maps to `comp`, up to certain
limits, such as only working for fields of a limited set of types. Generates IDs for the
type:

    a=1_b=2_c=xyz

where a, b, c are field names and the values following the equals signs are the
corresponding field values converted to strings.

# Arguments:
- `comp::Component` Input component
- `accesses=nothing`: If not nothing, selects which fields of `comp` are included in
  generating the ID
- `connector="_"`: String used to connect between field names and values
"""
function id(comp::Component; accesses=nothing, connector="_", kwargs...)
    savename(
        string(typeof(comp)),
        comp; 
        accesses=accesses === nothing ? fieldnames(typeof(comp)) : accesses,
        allowedtypes=(
            Real, 
            String, 
            Symbol, 
            Function,
            Component,
            Audiogram,
            Vector{<:Component},
        ), 
        connector=connector,
        kwargs...
    )
end

"""
    id(comps...[; kwargs...])

Return string uniquely mapping to combination of Components in `comps`
"""
function id(comps::Vararg{Component, N}; connector_super="_", connector="_", kwargs...) where N
    join(map(x -> id(x; connector=connector, kwargs...), comps), connector_super)
end

# We override DrWatson.access to give us recursive ID-generating superpowers
function DrWatson.access(comp::Component, key)
    if typeof(getproperty(comp, key)) <: Component
        id(getproperty(comp, key))
    else
        getproperty(comp, key)
    end
end

# Add a generic function for validating components 
validate(::Component) = print("Valid!")

# Add a custom display function for components
# function Base.display(comp::C) where {C <: Component}
#     # Print title
#     printstyled("$C <: $(supertype(C))\n"; bold=true, color=:blue)

#     # Print contents
#     names = fieldnames(C)
#     len_max = maximum(length.(string.(names)))
#     for name in names
#         printstyled("   $name: "; italic=true, bold=true)
#         display(getfield(comp, name))
#         print("\n")
#     end
# end