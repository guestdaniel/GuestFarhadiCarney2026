using Parameters
using Test
using Utilities

# # ==========================================================================================
# # Super simple dummy component, test basics
# # ==========================================================================================
@with_kw struct Dummy <: Component
    x::Float64=1.0
    y::Float64=2.0
    z::String="abc"
end

x = Dummy()
@test id(x) == "Dummy_x=1.0_y=2.0_z=abc"

x = Dummy()
y = Dummy(; x=9.2)
@test id(x, y) == "Dummy_x=1.0_y=2.0_z=abc_Dummy_x=9.2_y=2.0_z=abc"

# ==========================================================================================
# Recursive Dummy component, test recursion rules
# ==========================================================================================
@with_kw struct RDummy <: Component
    x::Float64=1.0
    y::Float64=2.0
    z::String="abc"
    d::Dummy=Dummy(5.0, 10.0, "qqq")
end

x = RDummy()
@test id(x) == "RDummy_d=Dummy_x=5.0_y=10.0_z=qqq_x=1.0_y=2.0_z=abc"