#include <Python.h>
#include <igraph_constructors.h>

#include "convert.h"

/**
 * \brief Converts a PyLong to an igraph \c igraph_integer_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * This function differs from the next one because it is less generic,
 * i.e. the Python object has to be a PyLong
 *
 * \param object the PyLong to be converted
 * \param v the result is stored here
 * \return 0 if everything was OK, 1 otherwise
 */
int PyLong_to_integer_t(PyObject* obj, igraph_integer_t* v) {
  if (IGRAPH_INTEGER_SIZE == 64) {
    /* here the assumption is that sizeof(long long) == 64 bits; anyhow, this
     * is the widest integer type that we can convert a PyLong to so we cannot
     * do any better than this */
    long long int dummy = PyLong_AsLongLong(obj);
    if (PyErr_Occurred()) {
      return 1;
    }
    *v = dummy;
  } else {
    /* this is either 32-bit igraph, or some weird, officially not-yet-supported
     * igraph flavour. Let's try to be on the safe side and assume 32-bit. long
     * ints are at least 32 bits so we will fit, otherwise Python will raise
     * an OverflowError on its own */
    long int dummy = PyLong_AsLong(obj);
    if (PyErr_Occurred()) {
      return 1;
    }
    *v = dummy;
  }
  return 0;
}

/**
 * \brief Converts a Python object to an igraph \c igraph_integer_t
 *
 * Raises suitable Python exceptions when needed.
 *
 * \param object the Python object to be converted
 * \param v the result is stored here
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_integer_t(PyObject *object, igraph_integer_t *v) {
  int retval;
  igraph_integer_t num;

  if (object == NULL) {
  } else if (PyLong_Check(object)) {
    retval = PyLong_to_integer_t(object, &num);
    if (retval)
      return retval;
    *v = num;
    return 0;
  } else if (PyNumber_Check(object)) {
    /* try to recast as PyLong */
    PyObject *i = PyNumber_Long(object);
    if (i == NULL)
      return 1;
    /* as above, plus decrement the reference for the temp variable */
    retval = PyLong_to_integer_t(i, &num);
    Py_DECREF(i);
    if (retval)
      return retval;
    *v = num;
    return 0;
  }
  PyErr_BadArgument();
  return 1;
}

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
 * \brief Converts a Python list of ints to an igraph \c igraph_vector_int_t
 * The incoming \c igraph_vector_int_t should be uninitialized.
 * Raises suitable Python exceptions when needed.
 *
 * This function is almost identical to
 * \ref igraphmodule_PyObject_to_vector_t . Make sure you fix bugs
 * in both cases (if any).
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_int_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vector_int_t(PyObject *list, igraph_vector_int_t *v) {
  PyObject *it = 0, *item;
  igraph_integer_t value = 0;
  Py_ssize_t i, j, k;
  int ok;

  if (PyUnicode_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integers");
    return 1;
  }

  if (!PySequence_Check(list)) {
    /* try to use an iterator */
    it = PyObject_GetIter(list);
    if (it) {
      if (igraph_vector_int_init(v, 0)) {
        //igraphmodule_handle_igraph_error();
        Py_DECREF(it);
        return 1;
      }

      while ((item = PyIter_Next(it)) != 0) {
        if (!PyNumber_Check(item)) {
          PyErr_SetString(PyExc_TypeError, "iterable must return numbers");
          ok = 0;
        } else {
          ok = (igraphmodule_PyObject_to_integer_t(item, &value) == 0);
        }

        if (ok == 0) {
          igraph_vector_int_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }

        if (igraph_vector_int_push_back(v, value)) {
          //igraphmodule_handle_igraph_error();
          igraph_vector_int_destroy(v);
          Py_DECREF(item);
          Py_DECREF(it);
          return 1;
        }

        Py_DECREF(item);
      }

      Py_DECREF(it);
    } else {
      PyErr_SetString(PyExc_TypeError, "sequence or iterable expected");
      return 1;
    }

    return 0;
  }

  j = PySequence_Size(list);

  if (igraph_vector_int_init(v, j)) {
    //igraphmodule_handle_igraph_error();
    Py_XDECREF(it);
    return 1;
  }

  for (i = 0, k = 0; i < j; i++) {
    item = PySequence_GetItem(list, i);
    if (item) {
      if (!PyNumber_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "sequence elements must be integers");
        ok = 0;
      } else {
        ok = (igraphmodule_PyObject_to_integer_t(item, &value) == 0);
      }

      Py_XDECREF(item);

      if (!ok) {
        igraph_vector_int_destroy(v);
        return 1;
      }

      VECTOR(*v)[k]=value;
      k++;
    } else {
      /* this should not happen, but we return anyway.
       * an IndexError exception was set by PyList_GetItem
       * at this point */
      igraph_vector_int_destroy(v);
      return 1;
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

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python integer list
 *
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList(const igraph_vector_t *v) {
  PyObject *list, *item;
  Py_ssize_t n, i;

  n = igraph_vector_size(v);
  if (n < 0) {
    return NULL;
  }

  list = PyList_New(n);
  if (!list) {
    return NULL;
  }

  for (i = 0; i < n; i++) {
    item = PyFloat_FromDouble(VECTOR(*v)[i]);

    if (!item) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, item);  /* will not fail */
  }

  return list;
}


/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_matrix_t to a Python list of lists
 *
 * \param m the \c igraph_matrix_t containing the matrix to be converted
 * \return the Python list of lists as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_matrix_t_to_PyList(const igraph_matrix_t *m) {
  PyObject *list, *row, *item;
  Py_ssize_t nr, nc, i, j;

  nr = igraph_matrix_nrow(m);
  nc = igraph_matrix_ncol(m);
  if (nc < 0 || nc < 0) {
    return NULL;
  }

  // create a new Python list
  list = PyList_New(nr);
  if (!list) {
    return NULL;
  }

  // populate the list with data
  for (i = 0; i < nr; i++) {
    row = PyList_New(nc);
    if (!row) {
      Py_DECREF(list);
      return NULL;
    }
  
    for (j = 0; j < nc; j++) {
      item = PyFloat_FromDouble(MATRIX(*m, i, j));
      if (!item) {
        Py_DECREF(row);
        Py_DECREF(list);
        return NULL;
      }

      PyList_SetItem(row, j, item);  /* will not fail */
    }

    PyList_SetItem(list, i, row);  /* will not fail */   
  }

  // return the list
  return list;
}


/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python iterable of non-negative integer pairs (i.e. an
 * edge list) to an igraph \c igraph_vector_t
 *
 * The incoming \c igraph_vector_int_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 *
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_int_t containing the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_edgelist(
    PyObject *list, igraph_vector_int_t *v
) {
  PyObject *item, *i1, *i2, *it;
  int ok;
  igraph_integer_t idx1=0, idx2=0;

  if (PyUnicode_Check(list)) {
    /* It is highly unlikely that a string (although it is a sequence) will
     * provide us with integers or integer pairs */
    PyErr_SetString(PyExc_TypeError, "expected a sequence or an iterable containing integer or string pairs");
    return 1;
  }

  it = PyObject_GetIter(list);
  if (!it)
    return 1;

  igraph_vector_int_init(v, 0);

  while ((item = PyIter_Next(it)) != 0) {
    ok = 1;
    if (!PySequence_Check(item) || PySequence_Size(item) != 2) {
      PyErr_SetString(PyExc_TypeError, "iterable must return pairs of integers or strings");
      ok = 0;
    } else {
      i1 = PySequence_GetItem(item, 0);
      i2 = i1 ? PySequence_GetItem(item, 1) : 0;
      ok = (i1 != 0 && i2 != 0);
      ok = ok && !igraphmodule_PyObject_to_integer_t(i1, &idx1);
      ok = ok && !igraphmodule_PyObject_to_integer_t(i2, &idx2);
      Py_XDECREF(i1); Py_XDECREF(i2); /* PySequence_ITEM returned new ref */
    }

    Py_DECREF(item);

    if (ok) {
      if (igraph_vector_int_push_back(v, idx1)) {
        //igraphmodule_handle_igraph_error();
        ok = 0;
      }
      if (ok && igraph_vector_int_push_back(v, idx2)) {
        //igraphmodule_handle_igraph_error();
        ok = 0;
      }
    }

    if (!ok) {
      igraph_vector_int_destroy(v);
      Py_DECREF(it);
      return 1;
    }
  }

  Py_DECREF(it);
  return 0;
}
