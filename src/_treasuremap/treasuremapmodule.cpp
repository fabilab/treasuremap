#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <igraph_constructors.h>
#include "treasuremap_layout.c"

namespace py = pybind11;

py::array_t<double> layout_treasuremap(
        py::array_t<double> seed,
        py::array_t<long> edges,
        py::array_t<double> distances,
        py::array_t<long> is_fixed,
        double min_dist,
        double sampling_probability,
        int epochs,
        int ndim) {
    
    long nedges, nvertices;
    py::buffer_info buf;
    double *ptr_double;
    long *ptr_long;
    igraph_vector_int_t igraph_edges;
    igraph_vector_t igraph_distances;
    igraph_vector_int_t igraph_is_fixed;
    igraph_t igraph_graph;
    igraph_matrix_t igraph_res;

    igraph_real_t igraph_sampling_prob = (igraph_real_t) sampling_probability;
    igraph_real_t igraph_min_dist = (igraph_real_t) min_dist;
    igraph_integer_t igraph_epochs = (igraph_integer_t) epochs;
    igraph_integer_t igraph_ndim = (igraph_integer_t) ndim;

    // TODO: copy data into igraph data structures...

    // edge list
    if ((edges.ndim() != 2) || (edges.shape(1) != 2)) {
        throw std::runtime_error("edges must be an Mx2 matrix of integers");
    }
    nedges = edges.shape(0);
    if (igraph_vector_int_init(&igraph_edges, 2 * nedges))
        throw std::runtime_error("Could not initialize edge list of graph");
    buf = edges.request();
    ptr_long = static_cast<long *>(buf.ptr);
    for (long i = 0; i < buf.shape[0]; i++)
        for (long j = 0; j < buf.shape[1]; j++)
            VECTOR(igraph_edges)[2 * i + j] = ptr_long[2 * i + j];

    // distances
    if ((distances.ndim() != 1) || (distances.size() != edges.shape(0))) {
        throw std::runtime_error("distances must be an vector of size equal to edges");
    }
    if(igraph_vector_init(&igraph_distances, nedges))
        throw std::runtime_error("Could not initialize distance list");
    buf = distances.request();
    ptr_double = static_cast<double *>(buf.ptr);
    for (long i = 0; i < buf.shape[0]; i++)
        VECTOR(igraph_distances)[i] = ptr_double[i];

    // is_fixed
    if (is_fixed.ndim() != 1) {
        throw std::runtime_error("is_fixed must be an vector of size equal to vertices");
    }
    nvertices = is_fixed.size();
    if (igraph_vector_int_init(&igraph_is_fixed, nvertices))
        throw std::runtime_error("Could not initialize fixed vertices list");
    buf = is_fixed.request();
    ptr_long = static_cast<long *>(buf.ptr);
    for (long i = 0; i < buf.shape[0]; i++)
        VECTOR(igraph_is_fixed)[i] = ptr_long[i];

    // output
    if ((seed.ndim() != 2) || (seed.shape(0) != nvertices) || (seed.shape(1) != ndim))
        throw std::runtime_error("res must be an N x ndim matrix of length equal to vertices");
    if (igraph_matrix_init(&igraph_res, nvertices, ndim))
        throw std::runtime_error("Could not initialize embedding matrix");

    // Initialize graph
    if (igraph_create(&igraph_graph, &igraph_edges, nvertices, false))
        throw std::runtime_error("Could not initialize graph");


    // Call C fuction
    igraph_error_t return_code = igraph_layout_treasuremap(
            &igraph_graph,
            &igraph_res,
            true,
            &igraph_distances,
            igraph_min_dist,
            igraph_epochs,
            igraph_sampling_prob,
            igraph_ndim,
            &igraph_is_fixed);

    // Destroy intermediate data structures
    igraph_destroy(&igraph_graph); // FIXME: requires more deallocation?
    igraph_vector_int_destroy(&igraph_is_fixed);
    igraph_vector_destroy(&igraph_distances);
    igraph_vector_int_destroy(&igraph_edges);

    if (return_code) {
        igraph_matrix_destroy(&igraph_res);
        throw std::runtime_error("Low-level C call failed");
    }

    // Convert output back to Python
    auto result = py::array_t<double>(nvertices * 2);
    buf = result.request();
    ptr_double = static_cast<double *>(buf.ptr);
    for (long i = 0; i < nvertices; i++)
        for (long j = 0; j < 2; j++)
            ptr_long[i * 2 + j] = MATRIX(igraph_res, i, j);
    // NOTE: C++ STL takes care of destructor when going out of scope
    std::vector<int> new_shape{(int)nvertices, 2};
    result = result.reshape(new_shape);

    // Destroy low-level object
    igraph_matrix_destroy(&igraph_res);

    return result;
}



PYBIND11_MODULE(_treasuremap, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: _treasuremap
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    m.def("_layout_treasuremap", &layout_treasuremap, R"pbdoc(
        TREASUREMAP algorithm.
    )pbdoc",
    py::arg("res"), py::arg("edges"), py::arg("distances"), py::arg("is_fixed"),
    py::arg("min_dist") = 0.01, py::arg("sampling_probability") = 0.3,
    py::arg("epochs") = 100, py::arg("ndim") = 2);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
