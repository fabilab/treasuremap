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


static PyObject* treasuremapmodule_compute_crossentropy(
        PyObject *self, PyObject *args, PyObject *kwds) {

    static char *kwlist[] = { "edges", "weights", "layout", "a", "b", NULL };
    PyObject *edges_o, *layout_o, *weights_o;
    double a, b;

    igraph_error_t return_code = 0;
    igraph_vector_int_t igraph_edges;
    igraph_t igraph_graph;
    igraph_matrix_t igraph_layout;
    igraph_vector_t igraph_weights;
    igraph_real_t cross_entropy;
    igraph_integer_t nvertices;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdd", kwlist,
                                    &edges_o, &weights_o, &layout_o, &a, &b))
        return NULL;

    // weights
    if (igraphmodule_PyObject_float_to_vector_t(weights_o, &igraph_weights)) {
        return NULL;
    }

    // layout
    if (igraphmodule_PyList_to_matrix_t(layout_o, &igraph_layout)) {
        igraph_vector_destroy(&igraph_weights);
        return NULL;
    }
    nvertices = igraph_matrix_size(&igraph_layout) / 2;

    // edgelist (already flattened)
    if (igraphmodule_PyObject_to_edgelist(edges_o, &igraph_edges)) {
        igraph_vector_destroy(&igraph_weights);
        igraph_matrix_destroy(&igraph_layout);
        return NULL;
    }

    // Initialize graph
    if (igraph_create(&igraph_graph, &igraph_edges, nvertices, false)) {
        igraph_vector_destroy(&igraph_weights);
        igraph_matrix_destroy(&igraph_layout);
        igraph_vector_int_destroy(&igraph_edges);
        return NULL;
    }
    igraph_vector_int_destroy(&igraph_edges);

    igraph_umap_compute_cross_entropy(
        &igraph_graph,
        &igraph_weights,
        &igraph_layout,
        (igraph_real_t) a,
        (igraph_real_t) b,
        &cross_entropy
    );

    // Destroy intermediate data structures
    igraph_destroy(&igraph_graph);
    igraph_matrix_destroy(&igraph_layout);
    igraph_vector_destroy(&igraph_weights);

    if (return_code != IGRAPH_SUCCESS) {
        return NULL;
    }

    return PyFloat_FromDouble((double) cross_entropy);

}

static PyObject* treasuremapmodule_compute_weights(PyObject *self, PyObject *args, PyObject *kwds) {

    static char *kwlist[] = {
        "nvertices",
        "nedges",
        "edges",
        "dist",
        NULL };
    long nedges, nvertices;
    PyObject *edges_o, *dist_o;
    PyObject *output;

    igraph_error_t return_code = 0;
    igraph_vector_int_t edges;
    igraph_vector_t distances;
    igraph_vector_t *distancesp = NULL;
    igraph_t graph;
    igraph_vector_t res;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "llOO", kwlist,
                                    &nvertices, &nedges,
                                    &edges_o, &dist_o))
        return NULL;

    /* Get distances */
    if (dist_o != Py_None) {
        if (igraphmodule_PyObject_float_to_vector_t(dist_o, &distances)) {
            return NULL;
        }
        distancesp = &distances;
    }

    // Prepare the output weights
    if (igraph_vector_init(&res, 0)) {
        igraph_vector_destroy(&distances);
        return NULL;
    }

    // edgelist (already flattened)
    if (igraphmodule_PyObject_to_edgelist(edges_o, &edges)) {
        igraph_vector_destroy(&res);
        if (dist_o != Py_None) igraph_vector_destroy(&distances);
        return NULL;
    }

    // Initialize graph: it is directed and typically asymmetric
    // The compute_weights function will symmetrize and assign
    // zero weight to redundant edges
    if (igraph_create(&graph, &edges, nvertices, true)) {
        igraph_vector_destroy(&res);
        igraph_vector_int_destroy(&edges);
        if (dist_o != Py_None) igraph_vector_destroy(&distances);
        return NULL;
    }
    igraph_vector_int_destroy(&edges);

    // Compute weights: this is to be understood as a symmetric
    // graph in which the redundant edges (distances that were
    // present in both direction) are now one edge with the joint
    // weight and one edge with negative weight.
    return_code = treasuremap_compute_weights(
            &graph,
            distancesp,
            &res
            );

    // Destroy intermediate data structures
    igraph_destroy(&graph);
    if (dist_o != Py_None) igraph_vector_destroy(&distances);

    if (return_code != IGRAPH_SUCCESS) {
        igraph_vector_destroy(&res);
        return NULL;
    }

    // Convert output back to Python
    output = igraphmodule_vector_t_to_PyList(&res);

    // Destroy low-level object
    igraph_vector_destroy(&res);

    return output;

}


static PyObject* treasuremapmodule_treasuremap(PyObject *self, PyObject *args, PyObject * kwds) {

    static char *kwlist[] = {
        "nvertices",
        "nedges",
        "edges",
        "dist",
        "seed",
        "is_fixed",
        "min_dist",
        "epochs",
        "dim",
        "negative_sampling_rate",
        "a", "b",
        "distances_are_weights",
        NULL };
    long nedges, nvertices;
    PyObject *edges_o, *dist_o = Py_None, *seed_o = Py_None, *is_fixed_o = Py_None;
    double min_dist;
    long epochs, ndim, negative_sampling_rate = 5;
    int distances_are_weights;
    PyObject *output;
    double a = -1, b = -1;

    igraph_error_t return_code = 0;
    igraph_vector_int_t edges;
    igraph_vector_t distances, weights;
    igraph_vector_bool_t is_fixed;
    igraph_vector_bool_t *is_fixedp = NULL;
    igraph_vector_t *distancesp = NULL, *weightsp = NULL;
    igraph_t graph;
    igraph_matrix_t res;

    /* Parse arguments */
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "llOOOOdll|lddp", kwlist,
                                    &nvertices, &nedges,
                                    &edges_o, &dist_o, &seed_o, &is_fixed_o,
                                    &min_dist, &epochs, &ndim,
                                    &negative_sampling_rate, &a, &b,
                                    &distances_are_weights
                                    ))
        return NULL;

    if ((ndim != 2) && (ndim != 3)) {
        PyErr_SetString(PyExc_ValueError, "Number of dimensions must be 2 or 3");
        return NULL;
    }

    if (min_dist < 0) {
        PyErr_SetString(PyExc_ValueError, "Minimum distance must be nonnegative");
        return NULL;
    }

    if (epochs < 0) {
        PyErr_SetString(PyExc_ValueError, "Epochs must be nonnegative");
        return NULL;
    }

    if (negative_sampling_rate < 0) {
        PyErr_SetString(PyExc_ValueError, "Negative sampling rate must be nonnegative");
        return NULL;
    }

    if (is_fixed_o != Py_None) {
        if (igraphmodule_PyObject_to_vector_bool_t(is_fixed_o, &is_fixed))
            return NULL;
        is_fixedp = &is_fixed;
    }

    if (dist_o != Py_None) {
        if (igraphmodule_PyObject_float_to_vector_t(dist_o, &distances)) {
            if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
            return NULL;
        }
        distancesp = &distances;
    }

    // With a seed, inject it into the output
    if (igraphmodule_PyList_to_matrix_t(seed_o, &res)) {
        if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
        if (dist_o != Py_None) igraph_vector_destroy(&distances);
        return NULL;
    }

    // edgelist (already flattened)
    if (igraphmodule_PyObject_to_edgelist(edges_o, &edges)) {
        igraph_matrix_destroy(&res);
        if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
        if (dist_o != Py_None) igraph_vector_destroy(&distances);
        return NULL;
    }

    // Initialize graph
    if (igraph_create(&graph, &edges, nvertices, false)) {
        igraph_vector_int_destroy(&edges);
        igraph_matrix_destroy(&res);
        if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
        if (dist_o != Py_None) igraph_vector_destroy(&distances);
        return NULL;
    }
    igraph_vector_int_destroy(&edges);

    // Compute weights, unless distances are already weights
    if (distances_are_weights) {
        weightsp = distancesp;
    } else {
        if (igraph_vector_init(&weights, 0)) {
            igraph_destroy(&graph);
            if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
            if (dist_o != Py_None) igraph_vector_destroy(&distances);
            return NULL;
        }
        if (treasuremap_compute_weights(&graph, distancesp, &weights)) {
            igraph_destroy(&graph);
            if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
            if (dist_o != Py_None) igraph_vector_destroy(&distances);
            igraph_vector_destroy(&weights);
            return NULL;
        }
        weightsp = &weights;
    }

    /* Fit a and b parameter to find smooth approximation to
     * probability distribution in embedding space */
    if ((a < 0) || (b < 0)) {
        return_code = fit_ab(min_dist, &a, &b);
        if (return_code != IGRAPH_SUCCESS) {
            igraph_destroy(&graph);
            if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
            if (dist_o != Py_None) igraph_vector_destroy(&distances);
            if (!distances_are_weights) igraph_vector_destroy(&weights);
            return NULL;
        }
    }

    // Compute layout
    return_code = igraph_layout_treasuremap(
            &graph,
            &res,
            weightsp,
            min_dist,
            epochs,
            (int) ndim,
            is_fixedp,
            a,
            b,
            negative_sampling_rate
            );

    // Destroy intermediate data structures
    igraph_destroy(&graph);
    if (is_fixed_o != Py_None) igraph_vector_bool_destroy(&is_fixed);
    if (dist_o != Py_None) igraph_vector_destroy(&distances);
    if (!distances_are_weights) igraph_vector_destroy(&weights);

    if (return_code != IGRAPH_SUCCESS) {
        igraph_matrix_destroy(&res);
        return NULL;
    }

    // Convert output back to Python
    output = igraphmodule_matrix_t_to_PyList(&res);

    // Destroy low-level object
    igraph_matrix_destroy(&res);

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
    {"compute_weights", (PyCFunction) treasuremapmodule_compute_weights, METH_VARARGS | METH_KEYWORDS, "Compute weights from distances"},
    {"compute_crossentropy", (PyCFunction) treasuremapmodule_compute_crossentropy, METH_VARARGS | METH_KEYWORDS, "Compute cross-entropy"},
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
