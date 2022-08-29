#include <Python.h>
#include <igraph_constructors.h>
#include <igraph_version.h>
#include "treasuremap_layout.c"


static PyObject* treasuremapmodule_treasuremap(PyObject *self, PyObject *args, PyObject * kwds) {

    static char *kwlist[] = { "nvertices", "nedges", "edges", "dist", "seed", "is_fixed", "min_dist", "sampling_prob", "epochs", "dim", NULL };
    long nedges, nvertices;
    PyObject *edges_o, *dist_o, *seed_o = Py_None, *is_fixed_o;
    double min_dist, sampling_probability;
    long epochs, ndim;

    igraph_vector_int_t igraph_edges;
    igraph_vector_t igraph_distances;
    igraph_vector_int_t igraph_is_fixed;
    igraph_t igraph_graph;
    igraph_matrix_t igraph_res;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "llOOOOddli", kwlist,
                                    &nvertices, &nedges,
                                    &edges_o, &dist_o, &seed_o, &is_fixed_o,
                                    &min_dist, &sampling_probability,
                                    &epochs, &ndim)) {
        return NULL;
    }


    return PyLong_FromLong(2);

}


static PyMethodDef treasuremapmodule_methods[] = {
    {"layout_treasuremap", (PyCFunction) treasuremapmodule_treasuremap, METH_VARARGS | METH_KEYWORDS, "Run the C-level treasuremap function"},
    {NULL, NULL, 0, NULL}
};


/**
 * Module definition table
 */
static struct PyModuleDef treasuremapmodule = {
    PyModuleDef_HEAD_INIT,
    "treasuremap._treasuremap",
    "CPython extension for treasuremap",
    -1,
    treasuremapmodule_methods,
};


/* Module initialization function */
PyMODINIT_FUNC
PyInit__treasuremap(void)
{
    PyObject* m;
   
    m = PyModule_Create(&treasuremapmodule);

    /* More useful constants */
    {
      const char* version;
      igraph_version(&version, 0, 0, 0);
      PyModule_AddStringConstant(m, "__igraph_version__", version);
    }
    PyModule_AddStringConstant(m, "__build_date__", __DATE__);

    return m;
}

//    // edge list
//    if (edges.size() > 0) {
//        if ((edges.ndim() != 2) || (edges.shape(1) != 2)) {
//            throw std::runtime_error("edges must be an Mx2 matrix of integers");
//        }
//    }
//    nedges = edges.size() / 2;
//    if (igraph_vector_int_init(&igraph_edges, 2 * nedges))
//        throw std::runtime_error("Could not initialize edge list of graph");
//    buf = edges.request();
//    ptr_long = static_cast<long *>(buf.ptr);
//    for (long i = 0; i < buf.shape[0]; i++)
//        for (long j = 0; j < buf.shape[1]; j++)
//            VECTOR(igraph_edges)[2 * i + j] = ptr_long[2 * i + j];
//
//    // distances
//    if (nedges > 0) {
//        if ((distances.ndim() != 1) || (distances.size() != edges.shape(0))) {
//            igraph_vector_int_destroy(&igraph_edges);
//            throw std::runtime_error("distances must be an vector of size equal to edges");
//        }
//    } else if (distances.size() != 0) {
//        igraph_vector_int_destroy(&igraph_edges);
//        throw std::runtime_error("inconsistent input: zero edges but nonzero distances");
//    }
//    if(igraph_vector_init(&igraph_distances, nedges))
//        throw std::runtime_error("Could not initialize distance list");
//    buf = distances.request();
//    ptr_double = static_cast<double *>(buf.ptr);
//    for (long i = 0; i < buf.shape[0]; i++)
//        VECTOR(igraph_distances)[i] = ptr_double[i];
//
//    // is_fixed
//    if (is_fixed.ndim() != 1) {
//        igraph_vector_int_destroy(&igraph_edges);
//        igraph_vector_destroy(&igraph_distances);
//        throw std::runtime_error("is_fixed must be an vector of size equal to vertices");
//    }
//    nvertices = is_fixed.size();
//    if (igraph_vector_int_init(&igraph_is_fixed, nvertices)){
//        igraph_vector_int_destroy(&igraph_edges);
//        igraph_vector_destroy(&igraph_distances);
//        throw std::runtime_error("Could not initialize fixed vertices list");
//    }
//    buf = is_fixed.request();
//    ptr_long = static_cast<long *>(buf.ptr);
//    for (long i = 0; i < buf.shape[0]; i++)
//        VECTOR(igraph_is_fixed)[i] = ptr_long[i];
//
//    // output
//    if (use_seed) {
//        if ((seed.ndim() != 2) || (seed.shape(0) != nvertices) || (seed.shape(1) != ndim)) {
//            igraph_vector_int_destroy(&igraph_is_fixed);
//            igraph_vector_int_destroy(&igraph_edges);
//            igraph_vector_destroy(&igraph_distances);
//            throw std::runtime_error("res must be an N x ndim matrix of length equal to vertices");
//        }
//    }
//    if (igraph_matrix_init(&igraph_res, nvertices, ndim)) {
//        igraph_vector_int_destroy(&igraph_is_fixed);
//        igraph_vector_int_destroy(&igraph_edges);
//        igraph_vector_destroy(&igraph_distances);
//        throw std::runtime_error("Could not initialize embedding matrix");
//    }
//    // If a seed is used (main operation), then it is taken as a pointer in the C function
//    // and altered in place into a result
//    if (use_seed) {
//        buf = seed.request();
//        ptr_double = static_cast<double *>(buf.ptr);
//        for (long i = 0; i < buf.shape[0]; i++)
//            for (long j = 0; j < buf.shape[1]; j++)
//                MATRIX(igraph_res, i, j) = ptr_long[2 * i + j];
//    }
//
//    // Initialize graph
//    if (igraph_create(&igraph_graph, &igraph_edges, nvertices, false)) {
//        igraph_matrix_destroy(&igraph_res);
//        igraph_vector_int_destroy(&igraph_is_fixed);
//        igraph_vector_int_destroy(&igraph_edges);
//        igraph_vector_destroy(&igraph_distances);
//        throw std::runtime_error("Could not initialize graph");
//    }
//
//    // Call C fuction
//    igraph_error_t return_code = igraph_layout_treasuremap(
//            &igraph_graph,
//            &igraph_res,
//            use_seed,
//            &igraph_distances,
//            igraph_min_dist,
//            igraph_epochs,
//            igraph_sampling_prob,
//            igraph_ndim,
//            &igraph_is_fixed);
//
//    // Destroy intermediate data structures
//    igraph_destroy(&igraph_graph);
//    igraph_vector_int_destroy(&igraph_is_fixed);
//    igraph_vector_destroy(&igraph_distances);
//    igraph_vector_int_destroy(&igraph_edges);
//
//    if (return_code != IGRAPH_SUCCESS) {
//        igraph_matrix_destroy(&igraph_res);
//        throw std::runtime_error("Low-level C call failed");
//    }
//
//    // Convert output back to Python
//    auto result = py::array_t<double>(nvertices * 2);
//    //buf = result.request();
//    //ptr_double = static_cast<double *>(buf.ptr);
//    //for (long i = 0; i < nvertices; i++)
//    //    for (long j = 0; j < 2; j++)
//    //        ptr_long[i * 2 + j] = MATRIX(igraph_res, i, j);
//
//    // Destroy low-level object
//    igraph_matrix_destroy(&igraph_res);
//
//    return result;
//}
//
//
//
//PYBIND11_MODULE(_treasuremap, m) {
//    m.doc() = R"pbdoc(
//        Pybind11 example plugin
//        -----------------------
//        .. currentmodule:: _treasuremap
//        .. autosummary::
//           :toctree: _generate
//           add
//           subtract
//    )pbdoc";
//
//    m.def("_layout_treasuremap", &layout_treasuremap, R"pbdoc(
//        TREASUREMAP algorithm.
//    )pbdoc",
//    py::arg("res"), py::arg("edges"), py::arg("distances"), py::arg("use_seed"), py::arg("is_fixed"),
//    py::arg("min_dist") = 0.01, py::arg("sampling_probability") = 0.3,
//    py::arg("epochs") = 100, py::arg("ndim") = 2);
//
//#ifdef VERSION_INFO
//    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
//#else
//    m.attr("__version__") = "dev";
//#endif
//}
