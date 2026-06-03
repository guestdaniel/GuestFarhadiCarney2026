x = zeros(5000)

ccall(
    (:accept_vector, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\demo\\demo.so"),
    Cvoid,
    (
        Ptr{Cdouble},
    ),
    x,
)

x = [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]]
ccall(
    (:accept_matrix, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\demo\\demo.so"),
    Cvoid,
    (
        Ptr{Ptr{Cdouble}},
    ),
    x,
)

ccall(
    (:allocate_matrix_complex2, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\demo\\demo.so"),
    Cvoid,
    (),
)

ccall(
    (:allocate_3d_array, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\demo\\demo.so"),
    Cvoid,
    (),
)

ccall(
    (:modify_vector_of_matrix_inplace, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\demo\\demo.so"),
    Cvoid,
    (),
)