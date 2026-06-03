x = [zeros(5, 5), zeros(5, 5)]
ccall(
    (:test_3d_load, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\test.so"),
    Cvoid,
    (Ptr{Ptr{Cdouble}}, Int64, Int64, Int64),
    pointer.(x),
    5,
    5,
    5,
)