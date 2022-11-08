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

#include "treasuremap_layout.h"

/* rho is just the size of the distance from each vertex and its closest neighbor */
/* sigma is the the decay from each vertex, depends on its rho and the rest of its
 * neighbor distances */

/* Find sigma for this vertex by binary search */
static igraph_error_t igraph_i_umap_find_sigma(const igraph_vector_t *distances,
        const igraph_vector_int_t *eids,
        igraph_real_t rho, igraph_real_t *sigma_p, igraph_real_t target) {

    igraph_real_t sigma = 1;
    igraph_real_t sum;
    igraph_real_t tol = 0.01;
    igraph_integer_t maxiter = 100;
    igraph_integer_t no_of_neis = igraph_vector_int_size(eids);
    igraph_integer_t eid;
    igraph_real_t step = sigma;
    igraph_integer_t seen_max = 0;

    /* Binary search */
    for (igraph_integer_t iter = 0; iter < maxiter; iter++) {
        sum = 0;
        for (igraph_integer_t j = 0; j < no_of_neis; j++) {
            eid = VECTOR(*eids)[j];
            sum += exp(-(VECTOR(*distances)[eid] - rho) / sigma);
        }

#ifdef SIGMA_DEBUG
        fprintf(stderr, "SIGMA function (no_of_neis = %" IGRAPH_PRId ")- sum: %g, "
               "target: %g, rho: %g, sigma: %g\n", no_of_neis, sum, target, rho, sigma);
#endif

        if (sum < target) {
            /* going back up after having seen an upper bound */
            if (seen_max == 1) {
                step /= 2;
            /* we need to go up but have not seen an upper bound yet
             * first iteration we want to increase by sigma, else we must come from
             * below, so we are sitting at 2 * step, we want to move to 4 * step */
            } else if (iter > 0) {
                step *= 2;
            }
            sigma += step;
        /* overshooting, we have definitely seen the max */
        } else {
            seen_max = 1;
            step /= 2;
            sigma -= step;
        }

        /* Check for convergence */
        if (fabs(sum - target) < tol) {
            break;
        }
    }

    *sigma_p = sigma;

    return IGRAPH_SUCCESS;
}


/* Convert the graph with distances into a probability graph with exponential decay.
 * The original UMAP implementation calls these probabilities "membership strengths"
 * and computes them in a separate function. Here, we have a single function computing
 * the sigmas (scale factors), rhos (minimal distance or connectivity correction)
 * and the probability graph in one go. */
igraph_error_t treasuremap_compute_weights(
        const igraph_t *graph,
        const igraph_vector_t *distances,
        igraph_vector_t *weights) {

    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_neis, eid, i, j, k, l;
    igraph_vector_int_t eids;
    igraph_real_t rho, dist_max, dist, sigma, weight, weight_inv, sigma_target;
    igraph_vector_int_list_t neighbors_seen;
    igraph_vector_list_t weights_seen;

    /* Check input distances */
    if (distances != NULL) {
        if (igraph_vector_size(distances) != no_of_edges) {
            IGRAPH_ERROR("Weights must be the same number as the edges in the graph.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t distance_min = igraph_vector_min(distances);
            if (distance_min < 0) {
                IGRAPH_ERROR("Distance vector must be nonnegative.", IGRAPH_EINVAL);
            } else if (isnan(distance_min)) {
                IGRAPH_ERROR("Distance vector must not contain NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    /* Resize weights */
    IGRAPH_CHECK(igraph_vector_resize(weights, no_of_edges));

    /* Trivial case */
    if (no_of_edges == 0) {
        return IGRAPH_SUCCESS;
    }

    /* Initialize vectors and matrices */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&neighbors_seen, no_of_vertices);
    IGRAPH_VECTOR_LIST_INIT_FINALLY(&weights_seen, no_of_vertices);

    /* Iterate over vertices x, like in the paper */
    for (i = 0; i < no_of_vertices; i++) {
        /* Edges into this vertex */
        IGRAPH_CHECK(igraph_incident(graph, &eids, i, IGRAPH_OUT));
        no_of_neis = igraph_vector_int_size(&eids);

        /* Vertex has no neighbors */
        if (no_of_neis == 0) {
            continue;
        }

        /* Find rho for this vertex, i.e. the minimal non-self distance */
        if (distances != NULL) {
            rho = VECTOR(*distances)[VECTOR(eids)[0]];
            dist_max = rho;
            for (j = 1; j < no_of_neis; j++) {
                dist = VECTOR(*distances)[VECTOR(eids)[j]];
                rho = fmin(rho, dist);
                dist_max = fmax(dist_max, dist);
            }
        } else {
            rho = dist_max = 0;
        }

        /* If the maximal distance is rho, all neighbors are nearest. */
        if (dist_max == rho) {
            /* This is a special flag for later on: it means weight = 1 */
            sigma = -1;

        /* Else, find sigma for this vertex, from its rho plus binary search */
        } else {
            sigma_target = log2(no_of_neis);
            IGRAPH_CHECK(igraph_i_umap_find_sigma(distances,
                        &eids, rho, &sigma,
                        sigma_target));
        }

        /* Convert to weight */
        for (j = 0; j < no_of_neis; j++) {
            eid = VECTOR(eids)[j];

            /* Basically, nodes closer than rho have probability 1, but nothing disappears */
            weight = sigma < 0 ? 1 : exp(-(VECTOR(*distances)[eid] - rho) / sigma);

            #ifdef UMAP_CONNECTIVITYDEBUG
            fprintf(stderr, "distance: %g\n", VECTOR(*distances)[eid]);
            fprintf(stderr, "weight: %g\n", weight);
            #endif

            /* Tag the edge as seen, for symmetrization */
            k = IGRAPH_OTHER(graph, eid, i);
            if (k == i) {
                IGRAPH_ERROR("Input graph must contain no loops.", IGRAPH_EINVAL);
            }

            igraph_vector_int_push_back(&(VECTOR(neighbors_seen)[i]), k);
            igraph_vector_push_back(&(VECTOR(weights_seen)[i]), weight);
        }

    }

    /* Symmetrize the adjacency matrix without casting the whole dense thing.
     * This takes the union of the fuzzy sets. In other words, compute the probability
     * (weight) of this edge being real as seen by either one vertex or the other:
     *
     *             P_{-> | <-} = P_{->} + P_{<-} - P_{->} * P_{<-}
     *
     * The original implementation calls this "set operations": union, i.e. what we do
     * here, is the default. One could also conceive doing the intersection (joint
     * probability) or a mix of union and intersection.
     *
     * NOTE: Symmetrization causes a loss of information, hence we need fewer edges
     * to represent connectivities than we have in the distance matrix which is
     * typically not symmetric. To avoid generating another graph altogether, we
     * encode the redundancy as a negative weight and skip that edge later on in
     * the gradient descent. */
    for (eid=0; eid < no_of_edges; eid++) {
        i = IGRAPH_FROM(graph, eid);
        k = IGRAPH_TO(graph, eid);

        /* Direct weight, if found */
        weight = 0;
        no_of_neis = igraph_vector_int_size(&(VECTOR(neighbors_seen)[i]));
        for (l=0; l < no_of_neis; l++) {
            if (VECTOR(VECTOR(neighbors_seen)[i])[l] == k) {
                weight = VECTOR(VECTOR(weights_seen)[i])[l];
                /* Tag this weight so we can ignore it later on if the opposite
                 * directed edge is found. It's ok to retag */
                VECTOR(VECTOR(weights_seen)[i])[l] = -1;
                break;
            }
        }

        if (weight < 0) {
            VECTOR(*weights)[eid] = 0;
            continue;
        }

        /* Weight of the opposite edge, if found */
        weight_inv = 0;
        no_of_neis = igraph_vector_int_size(&(VECTOR(neighbors_seen)[k]));
        for (l=0; l < no_of_neis; l++) {
            if (VECTOR(VECTOR(neighbors_seen)[k])[l] == i) {
                weight_inv = VECTOR(VECTOR(weights_seen)[k])[l];
                /* Tag this weight so we can ignore it later on if the opposite
                 * directed edge is found. It's ok to retag */
                VECTOR(VECTOR(weights_seen)[k])[l] = -1;
                break;
            }
        }

        if (weight_inv < 0) {
            VECTOR(*weights)[eid] = 0;
            continue;
        }

        VECTOR(*weights)[eid] = weight + weight_inv - weight * weight_inv;
    }

    igraph_vector_list_destroy(&weights_seen);
    igraph_vector_int_list_destroy(&neighbors_seen);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/* cross-entropy */
igraph_error_t igraph_umap_compute_cross_entropy(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        const igraph_matrix_t *layout,
        igraph_real_t a, igraph_real_t b,
        igraph_real_t *cross_entropy
        ) {

    igraph_real_t mu, nu, xd, yd, sqd;
    igraph_integer_t from, to;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_matrix_t edge_seen;

    IGRAPH_MATRIX_INIT_FINALLY(&edge_seen, no_of_vertices, no_of_vertices);

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. - sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * NOTE: the sum goes over the whole adjacency matrix, i.e. all potential edges,
     * not just the actual edges. That is because otherwise there's no benefit from
     * repelling unconnected edges.
     * */
    *cross_entropy = 0;
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        mu = VECTOR(*weights)[eid];
        if (mu < 0) {
            continue;
        }

        /* Find vertices */
        from = IGRAPH_FROM(graph, eid);
        to = IGRAPH_TO(graph, eid);
        /* Find distance in layout space */
        xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
        yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
        sqd = xd * xd + yd * yd;

        /* Find probability associated with distance using fitted Phi */
        nu = 1.0 / (1 + a * pow(sqd, b));

        /* DEBUG PRINT */
        if ((nu < 0.001) || (1 - nu < 0.001)) {
            if ((mu > 0) && (mu < 1)) {
                //fprintf(stderr, "CE, nu for eid %ld: %g, sqd: %g, mu: %g\n", eid, nu, sqd, mu);
            }
        }

        /* NOTE: we are in charge of corner cases where the log is trumped by the linear term */
        /* Term 1*/
        if (mu > 0)
            *cross_entropy -= mu * log(nu);

        /* Term 2*/
        if (mu < 1)
            *cross_entropy -= (1 - mu) * log(1 - nu);

        MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
    }
    /* Add the entropy from the missing edges */
    for (igraph_integer_t from = 0; from < no_of_vertices; from++) {
        for (igraph_integer_t to = 0; to < from; to++) {
            if (MATRIX(edge_seen, from, to) > 0) {
                continue;
            }

            /* Find distance in layout space */
            xd = (MATRIX(*layout, from, 0) - MATRIX(*layout, to, 0));
            yd = (MATRIX(*layout, from, 1) - MATRIX(*layout, to, 1));
            sqd = xd * xd + yd * yd;

            /* Find probability associated with distance using fitted Phi */
            nu = 1.0 / (1 + a * pow(sqd, b));

            /* Term 2*/
            *cross_entropy -= log(1 - nu);

            MATRIX(edge_seen, from, to) = MATRIX(edge_seen, to, from) = 1;
        }
    }

    igraph_matrix_destroy(&edge_seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* clip forces to avoid too rapid shifts */
static igraph_real_t igraph_i_umap_clip_force(igraph_real_t force, igraph_real_t limit) {
    if (force < -limit)
        return -limit;
    if (force > limit)
        return limit;
    return force;
}

static igraph_real_t igraph_i_umap_attract(
        igraph_real_t dsq,
        igraph_real_t a,
        igraph_real_t b)
{
    return - (2 * a * b * pow(dsq, b - 1.)) / (1. + a * pow(dsq, b));
}

static igraph_real_t igraph_i_umap_repel(
        igraph_real_t dsq,
        igraph_real_t a,
        igraph_real_t b)
{
    igraph_real_t dsq_min = CORRECT_DISTANCE_REPULSION * CORRECT_DISTANCE_REPULSION;
    return (2 * b) / (dsq_min + dsq) / (1. + a * pow(dsq, b));
}

static igraph_error_t igraph_i_umap_apply_forces(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_matrix_t *layout,
        igraph_real_t a,
        igraph_real_t b,
        igraph_real_t learning_rate,
        igraph_bool_t avoid_neighbor_repulsion,
        const igraph_vector_bool_t *is_fixed,
        igraph_integer_t negative_sampling_rate,
        igraph_integer_t epoch,
        igraph_integer_t epochs,
        igraph_vector_t *next_epoch_sample_per_edge)
{
    igraph_integer_t no_of_vertices = igraph_matrix_nrow(layout);
    igraph_integer_t ndim = igraph_matrix_ncol(layout);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t from, to, nneis;
    igraph_vector_t from_emb, to_emb, delta;
    igraph_real_t force = 0, dsq, force_d;
    /* The following is only used for small graphs, to avoid repelling your neighbors
     * For large sparse graphs, it's not necessary. For large dense graphs, you should
     * not be doing UMAP.
     * */
    igraph_vector_int_t neis, negative_vertices;
    igraph_integer_t n_negative_vertices = (no_of_vertices - 1 < negative_sampling_rate) ? (no_of_vertices - 1) : negative_sampling_rate;

    /* Initialize vectors */
    IGRAPH_VECTOR_INIT_FINALLY(&from_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&to_emb, ndim);
    IGRAPH_VECTOR_INIT_FINALLY(&delta, ndim);

    if (avoid_neighbor_repulsion) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&negative_vertices, 0);

    /* Iterate over edges. Stronger edges are sampled more often */
    for (igraph_integer_t eid = 0; eid < no_of_edges; eid++) {
        /* Zero or negative weights are a flag for skipping */
        if (VECTOR(*weights)[eid] <= 0) {
            continue;
        }

        /* We sample all and only edges that are supposed to be moved at this time */
        if ((VECTOR(*next_epoch_sample_per_edge)[eid] - epoch) >= 1) {
            continue;
        }

        /* set next epoch at which this edge will be sampled */
        VECTOR(*next_epoch_sample_per_edge)[eid] += 1.0 / VECTOR(*weights)[eid];

        /* we move all vertices on one end of the edges, then we come back for
         * the vertices on the other end. This way we don't move both ends at the
         * same time, which is almost a wasted move since they attract each other */
        int swapflag = (int)(RNG_UNIF01() > 0.5);
        int swapflag_end = swapflag + 2;
        for (; swapflag < swapflag_end; swapflag++) {

            /* half the time, swap the from/to, otherwise some vertices are never moved.
             * This has to do with the graph representation within igraph */
            if (swapflag % 2) {
                from = IGRAPH_FROM(graph, eid);
                to = IGRAPH_TO(graph, eid);
            } else {
                to = IGRAPH_FROM(graph, eid);
                from = IGRAPH_TO(graph, eid);
            }

            /* Skip if the node to move is fixed: it is neither attracted by
             * its neighbor "to", nor repelled by the negative sample of
             * vertices */
            if ((is_fixed != NULL) && (VECTOR(*is_fixed)[from] != 0)) {
                continue;
            }

            /* Current coordinates of both vertices */
            dsq = 0;
            for (igraph_integer_t d = 0; d != ndim; d++) {
                VECTOR(from_emb)[d] = MATRIX(*layout, from, d);
                VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
                VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
                dsq += VECTOR(delta)[d] * VECTOR(delta)[d];
            }

            /* Apply attractive force since they are neighbors */
            /* NOTE: If they are already together, no force needed */
            if (dsq >= MIN_DISTANCE_ATTRACTION * MIN_DISTANCE_ATTRACTION) {
                force = igraph_i_umap_attract(dsq, a, b);
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    force_d = force * VECTOR(delta)[d];
                    /* clip force to avoid too rapid change */
                    force_d = igraph_i_umap_clip_force(force_d, FORCE_LIMIT);

            #ifdef FORCES_DEBUG
                    fprintf(stderr, "force attractive: delta[%ld] = %g, forces[%ld] = %g\n", d, VECTOR(delta)[d], d, force_d);
            #endif

                    MATRIX(*layout, from, d) += learning_rate * force_d;
                }
            }

            /* Random other nodes repel the focal vertex */
            IGRAPH_CHECK(igraph_random_sample(&negative_vertices,
                        0, no_of_vertices - 2, n_negative_vertices));
            for (igraph_integer_t j = 0; j < n_negative_vertices; j++) {
                /* Get random neighbor */
                to = VECTOR(negative_vertices)[j];
                /* obviously you cannot repel yourself */
                if (to >= from) {
                    to++;
                }

                /* Fixed nodes do not repel, otherwise the new data will avoid them.
                 * In practice if you don't have this check you get a few stray
                 * points scattered around their actual cluster */
                if ((is_fixed != NULL) && (VECTOR(*is_fixed)[to] != 0)) {
                    continue;
                }

                /* do not repel neighbors for small graphs, for big graphs this
                 * does not matter as long as the k in knn << number of vertices */
                if (avoid_neighbor_repulsion) {
                    /* NOTE: the efficiency of this step could be improved but it
                     * should be only used for small graphs anyway, so it's fine */
                    igraph_bool_t skip = 0;
                    IGRAPH_CHECK(igraph_incident(graph, &neis, from, IGRAPH_ALL));
                    nneis = igraph_vector_int_size(&neis);
                    for (igraph_integer_t k = 0; k < nneis; k++) {
                        igraph_integer_t eid2 = VECTOR(neis)[k];
                        igraph_integer_t from2, to2;
                        from2 = IGRAPH_FROM(graph, eid2);
                        to2 = IGRAPH_TO(graph, eid2);
                        if (((from2 == from) && (to2 == to)) || ((from2 == to) && (from == to2))) {
                            skip = 1;
                            break;
                        }
                    }
                    if (skip == 1) {
                        continue;
                    }
                }

                /* Get layout of random neighbor and gradient in embedding */
                dsq = 0;
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    VECTOR(to_emb)[d] = MATRIX(*layout, to, d);
                    VECTOR(delta)[d] = MATRIX(*layout, from, d) - MATRIX(*layout, to, d);
                    dsq += VECTOR(delta)[d] * VECTOR(delta)[d];
                }

                /* This repels the other vertex assuming it's a negative example
                 * that is no weight, no edge */
                force = igraph_i_umap_repel(dsq, a, b);
                /* The repulsive force is already *away* from the other (non-neighbor) vertex */
                for (igraph_integer_t d = 0; d != ndim; d++) {
                    force_d = force * VECTOR(delta)[d];

                    /* clip force to avoid too rapid change */
                    force_d = igraph_i_umap_clip_force(force_d, FORCE_LIMIT);

                #ifdef FORCES_DEBUG
                    fprintf(stderr, "force repulsive: delta[%ld] = %g, forces[%ld] = %g\n", d, VECTOR(delta)[d], d, force_d);
                #endif

                    MATRIX(*layout, from, d) += learning_rate * force_d;
                }
            }
        }
    }

    /* Free vectors */
    igraph_vector_int_destroy(&negative_vertices);
    igraph_vector_destroy(&from_emb);
    igraph_vector_destroy(&to_emb);
    igraph_vector_destroy(&delta);
    IGRAPH_FINALLY_CLEAN(4);

    /* Free vector of neighbors if needed */
    if (avoid_neighbor_repulsion) {
        igraph_vector_int_destroy(&neis);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/* Edges with heavier weight/higher probability should be sampled more often. In
 * other words, vertices at each end of those edges should be moved more often. If the
 * edge weight is 1.0, which happens to each nearest neighbor due to the correction via
 * rho, that vertices at the end of that edge are moved each single epoch. Conversely,
 * vertices at the end of weak edges can be moved only once in a while. */
static igraph_error_t igraph_i_umap_optimize_layout_stochastic_gradient(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_real_t a,
        igraph_real_t b,
        igraph_matrix_t *layout,
        igraph_integer_t epochs,
        const igraph_vector_bool_t *is_fixed,
        igraph_integer_t negative_sampling_rate) {

    igraph_real_t learning_rate = 1;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_t next_epoch_sample_per_edge;

#ifdef UMAP_DEBUG
    igraph_real_t cross_entropy, cross_entropy_old;
#endif

    IGRAPH_VECTOR_INIT_FINALLY(&next_epoch_sample_per_edge, no_of_edges);

    /* Explicit avoidance of neighbor repulsion, only useful in small graphs
     * which are never very sparse. This is because negative sampling as implemented
     * relies on an approximation that only works if the graph is sparse, which is never
     * quite true for small graphs (i.e. |V| << |E| << |V|^2 is hard to judge if
     * |V| is small) */
    igraph_bool_t avoid_neighbor_repulsion = 0;
    if (igraph_vcount(graph) < 100) {
        avoid_neighbor_repulsion = 1;
    }

    /* Measure the (variable part of the) cross-entropy terms for debugging:
     * 1. - sum_edge_e mu(e) * log(nu(e))
     * 2. + sum_edge_e (1 - mu(e)) * log(1 - nu(e))
     * The latter is approximated by negative sampling as:
     * 2b. + sum_random_ij 1 * log(1 - nu_ij)
     * whereby the mu = 0 because we assume there's no edge between i and j, and nu_ij
     * is basically their distance in embedding space, lensed through the probability
     * function Phi.
     * */
#ifdef UMAP_DEBUG
    igraph_umap_compute_cross_entropy(
            graph, weights, layout, a, b, &cross_entropy);
#endif

    for (igraph_integer_t e = 0; e < epochs; e++) {
#ifdef UMAP_DEBUG
        fprintf(stderr, "Epoch %ld / %ld\n", e+1, epochs);
#endif

        /* Apply (stochastic) forces */
        igraph_i_umap_apply_forces(
                graph,
                weights,
                layout,
                a, b,
                learning_rate,
                avoid_neighbor_repulsion,
                is_fixed,
                negative_sampling_rate,
                e,
                epochs,
                &next_epoch_sample_per_edge);

#ifdef UMAP_DEBUG
        /* Recompute CE and check how it's going*/
        cross_entropy_old = cross_entropy;
        igraph_umap_compute_cross_entropy(
                graph, weights, layout, a, b, &cross_entropy);

        fprintf(stderr, "Cross-entropy before epoch: %g, after epoch: %g\n", cross_entropy_old, cross_entropy);
#endif

         /* Adjust learning rate */
        learning_rate = 1.0 - (igraph_real_t)(e + 1) / epochs;

        /* Allow interruption here */
        /* NOTE: this ends up being an implicit declaration which I guess is not great, but
         * I'm not sure how to include <Python.h> at this stage... */
        //if (PyErr_CheckSignals()) {
        if(igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) {
            igraph_vector_destroy(&next_epoch_sample_per_edge);
            return IGRAPH_INTERRUPTED;
        }
    }

    igraph_vector_destroy(&next_epoch_sample_per_edge);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* This is the main function that works for any dimensionality of the embedding
 * (currently hard-constrained to 2 or 3 ONLY in the initialization).
 *
 * NOTE: a and b are computed in the caller function by curve fitting (see above)
 *
 * */
igraph_error_t igraph_layout_treasuremap(
        const igraph_t *graph,
        igraph_matrix_t *res,
        const igraph_vector_t *weights,
        igraph_real_t min_dist,
        igraph_integer_t epochs,
        igraph_integer_t ndim,
        const igraph_vector_bool_t *is_fixed,
        igraph_real_t a,
        igraph_real_t b,
        igraph_integer_t negative_sampling_rate
        ) {

    igraph_error_t errorcode;
    igraph_integer_t no_of_vertices = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_fixed = 0;
    igraph_integer_t i;

    /* Distances should be of the right number and nonnegative */
    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weights must be the same number as the edges in the graph.", IGRAPH_EINVAL);
    }

    /* Count the number of fixed vertices, sanity check */
    if (is_fixed != NULL) {
        if (igraph_vector_bool_size(is_fixed) != no_of_vertices) {
            IGRAPH_ERROR("Fixed vertices must be None or a list with length equal to the number of vertices.", IGRAPH_EINVAL);
        }
        for (i = 0; i < no_of_vertices; i++) {
            if (VECTOR(*is_fixed)[i] != 0) {
                no_of_fixed++;
            }
        }
    }

    /* No edges, skip algorithm */
    if (no_of_edges == 0)
        return IGRAPH_SUCCESS;

    /* All vertices are fixed, skip algorithm */
    if (no_of_fixed == no_of_vertices)
        return IGRAPH_SUCCESS;

    RNG_BEGIN();

    /* Minimize cross-entropy between high-d and low-d probability
     * distributions */
    errorcode = igraph_i_umap_optimize_layout_stochastic_gradient(graph,
            weights,
            a, b,
            res,
            epochs,
            is_fixed,
            negative_sampling_rate);
    RNG_END();
    return errorcode;
}
