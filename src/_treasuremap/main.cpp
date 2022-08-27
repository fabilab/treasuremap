#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <igraph_constructors.h>
#include "treasuremap.c"

namespace py = pybind11;

int layout_treasuremap(
        py::array_t<double> res,
        py::array_t<int> edges,
        py::array_t<double> distances,
        py::array_t<int> is_fixed,
        double min_dist,
        double sampling_probability,
        int epochs,
        int ndim) {
    
    long nedges, nvertices;
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

    // distances
    if ((distances.ndim() != 1) || (distances.size() != edges.shape(0))) {
        throw std::runtime_error("distances must be an vector of size equal to edges");
    }
    if(igraph_vector_init(&igraph_distances, nedges))
        throw std::runtime_error("Could not initialize distance list");

    // is_fixed
    if (is_fixed.ndim() != 1) {
        throw std::runtime_error("is_fixed must be an vector of size equal to vertices");
    }
    nvertices = is_fixed.size();
    if (igraph_vector_int_init(&igraph_is_fixed, nvertices))
        throw std::runtime_error("Could not initialize fixed vertices list");

    // output
    if ((res.ndim() != 2) || (res.shape(0) != nvertices) || (res.shape(1) != ndim))
        throw std::runtime_error("res must be an N x ndim matrix of length equal to vertices");
    if (igraph_matrix_init(&igraph_res, nvertices, ndim))
        throw std::runtime_error("Could not initialize embedding matrix");
    //
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

    // Convert output back to Python
    // TODO

    // Destroy intermediate data structures
    igraph_matrix_destroy(&igraph_res);
    igraph_destroy(&igraph_graph); // FIXME: requires more deallocation?
    igraph_vector_int_destroy(&igraph_is_fixed);
    igraph_vector_destroy(&igraph_distances);
    igraph_vector_int_destroy(&igraph_edges);

    if (return_code) {
        // Raise Python exception
        return 1;
    }
    return 0;
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

    m.def("fit", &layout_treasuremap, R"pbdoc(
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
