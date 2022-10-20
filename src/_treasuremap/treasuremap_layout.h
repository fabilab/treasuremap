/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
   */

#ifndef TREASUREMAP_LAYOUT
#define TREASUREMAP_LAYOUT

#include <igraph_layout.h>
#include <igraph_interface.h>
#include <igraph_lapack.h>
#include <igraph_matrix.h>
#include <igraph_random.h>
#include <igraph_nongraph.h>
#include <igraph_interrupt.h>

#include <math.h>

// FIXME
#define UMAP_DEBUG

/* NOTE: support also a spread? */
igraph_error_t fit_ab(igraph_real_t min_dist, igraph_real_t *a_p, igraph_real_t *b_p);

/* This is the main function that works for any dimensionality of the embedding
 * (currently hard-constrained to 2 or 3 ONLY in the initialization). */
igraph_error_t igraph_layout_treasuremap(
        const igraph_t *graph,
        igraph_matrix_t *res,
        igraph_bool_t use_seed,
        const igraph_vector_t *distances,
        igraph_real_t min_dist,
        igraph_integer_t epochs,
        igraph_real_t sampling_prob,
        igraph_integer_t ndim,
        const igraph_vector_bool_t *is_fixed,
        igraph_real_t a,
        igraph_real_t b,
        igraph_integer_t negative_sampling_rate,
        int distances_are_connectivities
        );

#endif
