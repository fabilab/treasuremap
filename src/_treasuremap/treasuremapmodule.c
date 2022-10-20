#include <Python.h>
#include <igraph_constructors.h>
#include <igraph_version.h>

#include "convert.h"
#include "treasuremap_layout.h"


static PyObject* treasuremapmodule_fit_ab(PyObject *self, PyObject *args, PyObject *kwds) {

    static char *kwlist[] = { "min_dist", NULL };
    double min_dist, a, b;
    igraph_error_t return_code;
    PyObject *list, *item;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist,&min_dist)) {
        return NULL;
    }

    return_code = fit_ab(min_dist, &a, &b);
    if (return_code != IGRAPH_SUCCESS) {
        return NULL;
    }

    // Convert to Python list
    list = PyList_New(2);
    if (!list) {
      return NULL;
    }
    item = PyFloat_FromDouble(a);
    if (!item) {
        Py_DECREF(list);
        return NULL;
    }
    PyList_SetItem(list, 0, item);  /* will not fail */
    item = PyFloat_FromDouble(b);
    if (!item) {
        Py_DECREF(list);
        return NULL;
    }
    PyList_SetItem(list, 1, item);  /* will not fail */

    return list;
}


static PyObject* treasuremapmodule_treasuremap(PyObject *self, PyObject *args, PyObject * kwds) {

    static char *kwlist[] = { "nvertices", "nedges", "edges", "dist", "seed", "is_fixed", "min_dist", "sampling_prob", "epochs", "dim", "negative_sampling_rate", "a", "b", "distances_are_connectivities", NULL };
    long nedges, nvertices;
    PyObject *edges_o, *dist_o, *seed_o = Py_None, *is_fixed_o;
    double min_dist, sampling_probability;
    long epochs, ndim, negative_sampling_rate = 5;
    int distances_are_connectivities;
    PyObject *output;
    bool use_seed = false;
    double a = -1, b = -1;

    igraph_error_t return_code = 0;
    igraph_vector_int_t igraph_edges;
    igraph_vector_t igraph_distances;
    igraph_vector_bool_t igraph_is_fixed;
    igraph_t igraph_graph;
    igraph_matrix_t igraph_res;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "llOOOOddll|lddp", kwlist,
                                    &nvertices, &nedges,
                                    &edges_o, &dist_o, &seed_o, &is_fixed_o,
                                    &min_dist, &sampling_probability,
                                    &epochs, &ndim, &negative_sampling_rate, &a, &b,
                                    &distances_are_connectivities
                                    ))
        return NULL;

    if (igraphmodule_PyObject_to_vector_bool_t(is_fixed_o, &igraph_is_fixed))
        return NULL;

    if (igraphmodule_PyObject_float_to_vector_t(dist_o, &igraph_distances)) {
        igraph_vector_bool_destroy(&igraph_is_fixed);
        return NULL;
    }

    // With a seed, inject it into the output
    use_seed = (seed_o != Py_None);
    if (use_seed) {
        if (igraphmodule_PyList_to_matrix_t(seed_o, &igraph_res)) {
            igraph_vector_bool_destroy(&igraph_is_fixed);
            igraph_vector_destroy(&igraph_distances);
            return NULL;
        }
    // Otherwise, just make a zeroes matrix, it will be randomized
    // by the C function later on
    } else {
        if (igraph_matrix_init(&igraph_res, nvertices, ndim)) {
            igraph_vector_bool_destroy(&igraph_is_fixed);
            igraph_vector_destroy(&igraph_distances);
            return NULL;
        }
    }

    // edgelist (already flattened)
    if (igraphmodule_PyObject_to_edgelist(edges_o, &igraph_edges)) {
        igraph_matrix_destroy(&igraph_res);
        igraph_vector_bool_destroy(&igraph_is_fixed);
        igraph_vector_destroy(&igraph_distances);
        return NULL;
    }

    // Initialize graph
    if (igraph_create(&igraph_graph, &igraph_edges, nvertices, false)) {
        igraph_matrix_destroy(&igraph_res);
        igraph_vector_bool_destroy(&igraph_is_fixed);
        igraph_vector_int_destroy(&igraph_edges);
        igraph_vector_destroy(&igraph_distances);
        return NULL;
    }
    igraph_vector_int_destroy(&igraph_edges);

    /* Fit a and b parameter to find smooth approximation to
     * probability distribution in embedding space */
    if ((a < 0) || (b < 0)) {
        return_code = fit_ab(min_dist, &a, &b);
        if (return_code != IGRAPH_SUCCESS) {
            igraph_destroy(&igraph_graph);
            igraph_vector_bool_destroy(&igraph_is_fixed);
            igraph_vector_destroy(&igraph_distances);
            return NULL;
        }
    }

    // Call C fuction
    return_code = igraph_layout_treasuremap(
            &igraph_graph,
            &igraph_res,
            use_seed,
            &igraph_distances,
            min_dist,
            epochs,
            sampling_probability,
            (int) ndim,
            &igraph_is_fixed,
            a,
            b,
            negative_sampling_rate,
            distances_are_connectivities
            );

    // Destroy intermediate data structures
    igraph_destroy(&igraph_graph);
    igraph_vector_bool_destroy(&igraph_is_fixed);
    igraph_vector_destroy(&igraph_distances);

    if (return_code != IGRAPH_SUCCESS) {
        igraph_matrix_destroy(&igraph_res);
        return NULL;
    }

    // Convert output back to Python
    output = igraphmodule_matrix_t_to_PyList(&igraph_res);

    // Destroy low-level object
    igraph_matrix_destroy(&igraph_res);

    return output;
}


/* Connect interruption signals from inside igraph: these do not work directly
 * with interuptions in our C code */
static igraph_error_t treasuremapmodule_igraph_interrupt_hook(void* data) {
  if (PyErr_CheckSignals()) {
    IGRAPH_FINALLY_FREE();
    return IGRAPH_INTERRUPTED;
  }
  return IGRAPH_SUCCESS;
}


static PyMethodDef treasuremapmodule_methods[] = {
    {"layout_treasuremap", (PyCFunction) treasuremapmodule_treasuremap, METH_VARARGS | METH_KEYWORDS, "Run the C-level treasuremap function"},
    {"fit_ab", (PyCFunction) treasuremapmodule_fit_ab, METH_VARARGS | METH_KEYWORDS, "Fit a and b parameters from min_dist"},
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

    /* Connect interruptibles signal */
    igraph_set_interruption_handler(treasuremapmodule_igraph_interrupt_hook);

    return m;
}
