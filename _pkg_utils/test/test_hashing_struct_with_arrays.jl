using Test
using Parameters
using Utilities

# Start with a simple Dummy component that just includes a float and a vector field
@with_kw struct DummyVector <: Component
    x::Float64=1.0
    y::Vector{Float64}=rand(50)
end

# When we force vector to be the same, ensure hash is the same
x = DummyVector(; x=1.0, y=collect(1:50))
y = DummyVector(; x=1.0, y=collect(1:50))
@test x == y
@test hash(x) == hash(y)

# # When we don't, ensure they are different
# x = Dummy()
# y = Dummy()
# @test x != y
# @test hash(x) != hash(y)

# # Create thousands and thousands of dummies, ensure none hash to same value
# hashes = [hash(Dummy()) for ]
