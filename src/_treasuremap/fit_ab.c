
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

/* Helper function to compute a and b parameters (smoothing probability metric in
 * embedding space). The official UMAP package has a parameter "spread" to set the
 * scaling, but it defaults to 1.0. Here, we assume 1.0 for now, could be extended
 * in the future */
static igraph_error_t igraph_i_umap_get_ab_residuals(igraph_vector_t *residuals,
        igraph_real_t *squared_sum_res, igraph_integer_t nr_points, igraph_real_t a,
        igraph_real_t b, igraph_vector_t *powb, const igraph_vector_t *x, igraph_real_t min_dist)
{
    igraph_real_t tmp;

    *squared_sum_res = 0;
    for (igraph_integer_t i = 0; i < nr_points; i++) {
        /* The ideal probability is:
         *
         *     P(d) = d < min_dist ? 1 : e^{-(d - min_dist)}
         *
         * which is the same as the high-dimensional probability, except
         * min_dist plays the role of rho and sigma is fixed at 1. However,
         * this function has a kink at min_dist (first derivative is not
         * continuous). So we smoothen it with:
         *
         *     Q(d) = ( 1 + a*d^2b )^-1
         *
         * which is quite similar throughout for appropriate a and b. Notice
         * that we do not need to smoothen the high-dimensional probability
         * function because the vertices are not moved in the high-dimensional
         * space, so there is no need for differentiating that function.
         *
         * The residual is of course:
         *
         *    Q(d) - P(d) = ( 1 + a*d^2b )^-1 - [ d < min_dist ? 1 : e^{-(d - min_dist)} ]
         *
         * This function also sets the auxiliary vector powb.
         * */
        VECTOR(*powb)[i] = pow(VECTOR(*x)[i], 2 * b);
        tmp = 1 / (1 + a * VECTOR(*powb)[i]);
        tmp -= VECTOR(*x)[i] <= min_dist ? 1 : exp(-(VECTOR(*x)[i] - min_dist));
        VECTOR(*residuals)[i] = tmp;
        *squared_sum_res += tmp * tmp;
    }
    return IGRAPH_SUCCESS;
}

/* UMAP minimizes the cross-entropy between probability of being a true edge in
 * high and low dimensions. For the low-dimensional computation, it uses a smooth
 * function of the Euclidean distance between two vertices:
 *
 * P(d) = (1 + a*d^2b)^-1
 *
 * where d is the distance and a and b are hyperparameters that basically determine
 * the cutoff distance at which the probability starts to decrease.
 *
 * We fit these two parameters using nonlinear least squares (Gauss-Newton + line search)
 * on a grid of artificial distances. There is only one user-chosen input argument that
 * determines this fit, called min_dist, which is approximately the cutoff distance we
 * are trying to achieve.
 *
 * ADVANCED NOTE:
 * In a way, the whole UMAP layout is invariant upon scaling transformations, of course,
 * so min_dist is basically meaningless. Another way to see this is that for any pair
 * (a,b) that minimize the least squares for dist_min, we can easily find a solution for
 * a new dist_min2 := alpha * dist_min:
 *
 * P(d, a, b) = (1 + a*d^2b)^-1
 *
 * P(alpha * d, a', b') = (1 + a'*(alpha * d)^2b' )^-1
 *
 * that is:
 *
 * a*d^2b = a'*alpha^2b'*d^2b'   for each  d >= 0.
 *
 * So for d = 1        ->  a = a'*alpha^2b'
 * and for d = sqrt(2) ->  a*2^b = a'*alpha^2b'*2^b'
 *
 * which solves as:
 *
 * b' = b
 * a' = a / alpha^2b
 *
 * For instance, if b = 1, a -> 0.01*a moves the fit a decade towards larger min_dist,
 * and a -> 100*a moves the fit a decade towards smaller min_dist.
 * */
igraph_error_t fit_ab(igraph_real_t min_dist, igraph_real_t *a_p, igraph_real_t *b_p)
{
    /* Grid points */
    igraph_vector_t x;
     /* Make a lattice from 0 to this distance */
    igraph_integer_t nr_points = 100;
    igraph_real_t end_point = min_dist*30; /* Need to sample decently around the kink I guess? */
    /* Initial values: a^-2b is the point where f(x) = 0.5, so around min_dist
     * Therefore, initial a is 1/sqrt(min_dist) */
    igraph_real_t b = 1.;
    igraph_real_t a = 1. / sqrt(min_dist);
    /* deltas */
    igraph_real_t da, db;
    /* Residuals */
    igraph_vector_t residuals;
    igraph_real_t squared_sum_res, squared_sum_res_old, squared_sum_res_tmp;
    /* Needed for the Gauss-Newton search */
    igraph_matrix_t jacobian, jTj, jTr;
    igraph_real_t tol = 0.001;
    igraph_real_t maxiter = 100;
    /* Auxiliary vars */
    igraph_real_t tmp;
    igraph_vector_t powb;
    int lapack_info;

    /* Distance lattice */
    IGRAPH_VECTOR_INIT_FINALLY(&x, nr_points);
    /* Residuals */
    IGRAPH_VECTOR_INIT_FINALLY(&residuals, nr_points);
    /* First derivatives, for the fitting (direction) */
    IGRAPH_MATRIX_INIT_FINALLY(&jacobian, nr_points, 2);
    /* Composite matrices/vectors for linear least squares at each iteration */
    IGRAPH_MATRIX_INIT_FINALLY(&jTj, 2, 2);
    IGRAPH_MATRIX_INIT_FINALLY(&jTr, 2, 1);
    /* Auxiliary vars for convenience */
    IGRAPH_VECTOR_INIT_FINALLY(&powb, nr_points);

    /* Distance |x-y| (this is a lattice, there are no actual x and y) */
    for (igraph_integer_t i = 0; i < nr_points; i++) {
        VECTOR(x)[i] = (end_point / nr_points) * i + 0.001; /* added a 0.001 to prevent NaNs */
    }

#ifdef UMAP_DEBUG
    printf("start fit_ab\n");
#endif
    for (igraph_integer_t iter = 0; iter < maxiter; iter++) {
        IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a, b,
                    &powb, &x, min_dist));

        /* break if good fit (conergence to truth) */
        if (squared_sum_res < tol * tol) {
#ifdef UMAP_DEBUG
            printf("convergence to zero (wow!)\n");
#endif
            break;
        }
        /* break if no change (convergence) */
        if ((iter > 0) && fabs(sqrt(squared_sum_res_old) - sqrt(squared_sum_res)) < tol) {
#ifdef UMAP_DEBUG
            printf("no-change absolute convergence\n");
#endif
            break;
        }

        /* Jacobian (first derivatives) of squared residuals at (a, b) */
        for (igraph_integer_t i = 0; i < nr_points; i++) {
            tmp = 1 + a * VECTOR(powb)[i];
            MATRIX(jacobian, i, 0) = - 2 * VECTOR(powb)[i] / tmp / tmp;
            MATRIX(jacobian, i, 1) = MATRIX(jacobian, i, 0) * a * log(VECTOR(x)[i]) * 2;
        }

        /* At each iteration, we want to minimize the linear approximation of the sum of squared
         * residuals:
         *
         * sum_i (Ji @ d(a,b) -r_i)^2
         *
         * Putting the first derivative to zero results in a linear system of 2 equations
         * (for a and b):
         *
         * sum_i J_i^T @ J_i @ d(a,b) = sum_i J_i^T r_i
         * *
         * or more compactly:
         *
         * J^T @ J @ d(a,b) = J^T @ r
         *
         * where J_T is the transpose of the Jacobian. Defining A := J^T @ J, B = J^T @ r:
         *
         * A @ d(a,b) = B
         *
         * This can be solved for d(a,b) using LAPACK within igraph
         * */
        /* Compute A and B, i.e. J^T @ J and J^T @ r */
        MATRIX(jTj, 0, 0) = MATRIX(jTj, 0, 1) = MATRIX(jTj, 1, 0) = MATRIX(jTj, 1, 1) = 0;
        MATRIX(jTr, 0, 0) = MATRIX(jTr, 1, 0) = 0;
        for (igraph_integer_t i = 0; i < nr_points; i++) {
            for (igraph_integer_t j1 = 0; j1 < 2; j1++) {
                for (igraph_integer_t j2 = 0; j2 < 2; j2++) {
                    MATRIX(jTj, j1, j2) += MATRIX(jacobian, i, j1) * MATRIX(jacobian, i, j2);
                }
                MATRIX(jTr, j1, 0) += MATRIX(jacobian, i, j1) * VECTOR(residuals)[i];
            }
        }
        /* LAPACK puts solution into jTr */
        IGRAPH_CHECK(igraph_lapack_dgesv(&jTj, 0, &jTr, &lapack_info));

        /* This might go wrong, in which case we should fail graciously */
        if (lapack_info != 0) {
            igraph_vector_destroy(&x);
            igraph_vector_destroy(&residuals);
            igraph_matrix_destroy(&jacobian);
            igraph_matrix_destroy(&jTj);
            igraph_matrix_destroy(&jTr);
            igraph_vector_destroy(&powb);
            IGRAPH_FINALLY_CLEAN(6);
            IGRAPH_ERROR("Singular matrix in the estimation of a and b for UMAP",  IGRAPH_EINVAL);
        }

        da = -MATRIX(jTr, 0, 0);
        db = -MATRIX(jTr, 1, 0);

        /* Improvement over GN: rough exponential line search for best delta
         * start from largest change, and keep shrinking as long as we are going down
         * */
        squared_sum_res_old = squared_sum_res;
        IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points, a + da,
                    b + db, &powb, &x, min_dist));

#ifdef UMAP_DEBUG
        printf("start line search, SSR before delta: %g, current SSR:, %g\n", squared_sum_res_old,
                squared_sum_res);
#endif
        for (igraph_integer_t k = 0; k < 30; k++) {
            /* Try new parameters */
            da /= 2.0;
            db /= 2.0;
            squared_sum_res_tmp = squared_sum_res;
            IGRAPH_CHECK(igraph_i_umap_get_ab_residuals(&residuals, &squared_sum_res, nr_points,
                        a + da, b + db, &powb, &x, min_dist));

            /* Compare and if we are going back uphill, undo last step and break */
#ifdef UMAP_DEBUG
            printf("during line search, k = %ld, old SSR:, %g, new SSR (half a,b):, %g\n", k,
                    squared_sum_res_tmp, squared_sum_res);
#endif
            if (squared_sum_res > squared_sum_res_tmp - tol) {
                da *= 2;
                db *= 2;
                break;
            }
        }
#ifdef UMAP_DEBUG
        printf("end of line search and iteration, squared_sum_res: %g \n\n", squared_sum_res_tmp);
#endif

        /* assign a, b*/
        a += da;
        b += db;

    }

    /* Free memory and tidy up stack */
    igraph_vector_destroy(&x);
    igraph_vector_destroy(&residuals);
    igraph_matrix_destroy(&jacobian);
    igraph_matrix_destroy(&jTj);
    igraph_matrix_destroy(&jTr);
    igraph_vector_destroy(&powb);
    IGRAPH_FINALLY_CLEAN(6);

#ifdef UMAP_DEBUG
    printf("a, b: %g %g\n", a, b);
#endif

    *a_p = a;
    *b_p = b;

    return IGRAPH_SUCCESS;

}
