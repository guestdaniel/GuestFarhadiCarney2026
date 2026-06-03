# Utilities.jl
Toolbox that binds together a lot of code I frequently reuse in Carney lab Julia code and provides an interface to our various experimental models.

# Design sketch
A few ubiquitous design features:

## Nesting plotting functions
Plotting functions are as `cplot(stuff..., fig=Figure(), ax=Axis(fig[1, 1]))` almost wherever possible, and then return `return fig, ax` wherever possible.
This style makes it easy to embed larger plots as subplots in a larger layout or to adjust features of a figure or axis object *post hoc*. 

## Components
Utilities provides a `Component` abstract type that other types that are designed to implement model features/components/analyses/etc are encouraged to extend.
A single `Component` could encapsulate an entire elaborate analysis, or merely a single model component.
Then a loose functional interface relates different components together, as in calls like
```julia
m = Model()
s = Stimulus()

ep = run(m, s)
```

```julia
stim = Stimulus()
model = Model()
summary = ExcitationPattern()

ep = run(model, stim)
```

### Interface
The following functional interfaces are provided for `Component`s:

- `id(x)` or `id(x1, x2, ...)` returns strings that uniquely identify the combination of `Component`s in question (up to the limit imposed by the fact that some fields of a `Component` may not be convertible into a string, such as vectors or matrices of data)