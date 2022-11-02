/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Tamas Nepusz <ntamas@gmail.com>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

/************************ Miscellaneous functions *************************/

/** \defgroup python_interface_conversion Converting between Python and igraph data types
 * \ingroup python_interface */

#ifndef TREASUREMAP_CONVERT_H
#define TREASUREMAP_CONVERT_H

#include <Python.h>
#include <igraph_constructors.h>

/* Conversion from PyObject to igraph types */

int igraphmodule_PyObject_to_integer_t(PyObject *object, igraph_integer_t *v);
int igraphmodule_PyObject_to_real_t(PyObject *object, igraph_real_t *v);

int igraphmodule_PyObject_to_vector_t(PyObject *list, igraph_vector_t *v,
		igraph_bool_t need_non_negative);
int igraphmodule_PyObject_float_to_vector_t(PyObject *list, igraph_vector_t *v);
int igraphmodule_PyObject_to_vector_int_t(PyObject *list, igraph_vector_int_t *v);
int igraphmodule_PyObject_to_vector_bool_t(PyObject *list, igraph_vector_bool_t *v);

int igraphmodule_PyList_to_matrix_t(PyObject *o, igraph_matrix_t *m);

PyObject* igraphmodule_vector_t_to_PyList(const igraph_vector_t *v);
PyObject* igraphmodule_matrix_t_to_PyList(const igraph_matrix_t *m);

int igraphmodule_PyObject_to_edgelist(PyObject *list, igraph_vector_int_t *v);
#endif
