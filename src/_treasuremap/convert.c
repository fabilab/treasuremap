#include <Python.h>
#include <igraph_constructors.h>

/**
 * \brief Converts a Python object to an igraph \c igraph_real_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * \param object the Python object to be converted; \c NULL is accepted but
 *        will keep the input value of v
 * \param v the result is returned here
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_real_t(PyObject *object, igraph_real_t *v) {
  igraph_real_t value;

  if (object == NULL) {
    return 0;
  } else if (PyLong_Check(object)) {
    value = PyLong_AsDouble(object);
  } else if (PyFloat_Check(object)) {
    value = PyFloat_AsDouble(object);
  } else if (PyNumber_Check(object)) {
    PyObject *i = PyNumber_Float(object);
    if (i == NULL) {
      return 1;
    }
    value = PyFloat_AsDouble(i);
    Py_DECREF(i);
  } else {
    PyErr_BadArgument();
    return 1;
  }

  if (PyErr_Occurred()) {
    return 1;
  } else {
    *v = value;
    return 0;
  }
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of floats to an igraph \c igraph_vector_t
 * The incoming \c igraph_vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_float_to_vector_t(PyObject *list, igraph_vector_t *v) {
  PyObject *item, *it;
  Py_ssize_t size_hint;
  int ok;
  igraph_real_t number;

  if (PyUnicode_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with numbers */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing numbers");
    return 1;
  }

  /* if the list is a sequence, we can pre-allocate the vector to its length */
  if (PySequence_Check(list)) {
    size_hint = PySequence_Size(list);
    if (size_hint < 0) {
      /* should not happen but let's try to recover */
      size_hint = 0;
    }
  } else {
    size_hint = 0;
  }

  /* initialize the result vector */
  if (igraph_vector_init(v, 0)) {
    //igraphmodule_handle_igraph_error();
    return 1;
  }

  /* if we have a size hint, allocate the required space */
  if (size_hint > 0) {
    if (igraph_vector_reserve(v, size_hint)) {
      //igraphmodule_handle_igraph_error();
      igraph_vector_destroy(v);
      return 1;
    }
  }

  /* try to use an iterator first */
  it = PyObject_GetIter(list);
  if (it) {
    while ((item = PyIter_Next(it)) != 0) {
      ok = 1;

      if (igraphmodule_PyObject_to_real_t(item, &number)) {
        PyErr_SetString(PyExc_ValueError, "iterable must yield numbers");
        ok=0;
      }

      Py_DECREF(item);

      if (!ok) {
        igraph_vector_destroy(v);
        Py_DECREF(it);
        return 1;
      }

      if (igraph_vector_push_back(v, number)) {
        //igraphmodule_handle_igraph_error();
        igraph_vector_destroy(v);
        Py_DECREF(it);
        return 1;
      }
    }
    Py_DECREF(it);
  } else {
    /* list is not iterable; maybe it's a single number? */
    PyErr_Clear();
    if (igraphmodule_PyObject_to_real_t(list, &number)) {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      igraph_vector_destroy(v);
      return 1;
    } else {
      igraph_vector_push_back(v, number);
    }
  }

  return 0;
}


/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of objects to an igraph \c igraph_vector_bool_t
 * The incoming \c igraph_vector_bool_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_bool_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_bool_t(PyObject *list,
    igraph_vector_bool_t *v) {
  PyObject *item;
  Py_ssize_t i, j;

  if (PyUnicode_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    PyObject *it = PyObject_GetIter(list);
    if (it) {
      PyObject *item;
      igraph_vector_bool_init(v, 0);
      while ((item = PyIter_Next(it)) != 0) {
        if (igraph_vector_bool_push_back(v, PyObject_IsTrue(item))) {
          //igraphmodule_handle_igraph_error();
          igraph_vector_bool_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }
        Py_DECREF(item);
      }
      Py_DECREF(it);
      return 0;
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }
    return 0;
  }

  j=PySequence_Size(list);
  igraph_vector_bool_init(v, j);
  for (i=0; i<j; i++) {
    item=PySequence_GetItem(list, i);
    if (item) {
      VECTOR(*v)[i]=PyObject_IsTrue(item);
      Py_DECREF(item);
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PySequence_GetItem
       * at this point */
      igraph_vector_bool_destroy(v);
      return 1;
    }
  }
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of lists to an \c igraph_matrix_t
 *
 * \param o the Python object representing the list of lists
 * \param m the address of an uninitialized \c igraph_matrix_t
 * \return 0 if everything was OK, 1 otherwise. Sets appropriate exceptions.
 */
int igraphmodule_PyList_to_matrix_t(PyObject* o, igraph_matrix_t *m) {

  Py_ssize_t nr, nc, n, i, j;
  PyObject *row, *item;
  int was_warned = 0;

  /* calculate the matrix dimensions */
  if (!PySequence_Check(o) || PyUnicode_Check(o)) {
    PyErr_SetString(PyExc_TypeError, "matrix expected (list of sequences)");
    return 1;
  }

  nr = PySequence_Size(o);
  nc = 0;
  for (i = 0; i < nr; i++) {
    row = PySequence_GetItem(o, i);
    if (!PySequence_Check(row)) {
      Py_DECREF(row);
      PyErr_SetString(PyExc_TypeError, "matrix expected (list of sequences)");
      return 1;
    }
    n = PySequence_Size(row);
    Py_DECREF(row);
    if (n > nc) {
      nc = n;
    }
  }

  igraph_matrix_init(m, nr, nc);
  for (i = 0; i < nr; i++) {
    row = PySequence_GetItem(o, i);
    n = PySequence_Size(row);
    for (j = 0; j < n; j++) {
      item = PySequence_GetItem(row, j);
      if (PyLong_Check(item)) {
        MATRIX(*m, i, j) = PyLong_AsLong(item);
      } else if (PyLong_Check(item)) {
        MATRIX(*m, i, j) = PyLong_AsLong(item);
      } else if (PyFloat_Check(item)) {
        MATRIX(*m, i, j) = PyFloat_AsDouble(item);
      } else if (!was_warned) {
        //PY_IGRAPH_WARN("non-numeric value in matrix ignored");
        was_warned=1;
      }
      Py_DECREF(item);
    }
    Py_DECREF(row);
  }

  return 0;
}


