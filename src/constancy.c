
#include "constancy.h"
#include "families.h"

// get directions of constancy as sparse matrix, each row a direction
//     of constancy
// pred, group, code, and delta are repred, regroup, recode, and redelta
//     in the asterdata object, all have length nnode
// want_phi is TRUE if we want directions of constancy for unconditional
//     canonical affine submodel and FALSE for conditional
// doit is FALSE if we do not have result vectors allocated yet because
//     we don't know how long they should be, in which case ix, jx, xx
//     are not touched and the length they need to be is returned
// doit is TRUE if we do have result vectors allocated of length nnonzero
//     in which case ix, jx, and xx are modified so that on return they
//     say that for the output sparse matrix the element in row ix[k]
//     and column jx[k] has value xx[k] and all elements not specified
//     this way are zero, and only arguments ix, jx, and xx are modified.
//     Note that ix and jx are 1-origin indexing (for R) not 0-origin as in C.
static int get_elements(int nnode, int *pred, int *group, int *code,
    double *delta, _Bool want_phi, _Bool doit, int nnonzero,
    int *ix, int *jx, double *xx)
{
    int knonzero = 0;
    int nrow = 0;
    _Bool todo[nnode];
    // is_zero[i] means the i-th node is almost surely zero because of
    //     a direction of constancy of the dependence group containing it
    // ancestor_is_zero[i] means the i-th node is almost surely zero
    //     because its predecessor is almost surely zero
    // we need the distinction because directions of constancy for
    //     nodes i such that is_zero[i] are done in the first pass
    //     but those such that ancestor_is_zero[i] need to be done later
    _Bool is_zero[nnode];
    _Bool ancestor_is_zero[nnode];
    for (int i = 0; i < nnode; i++) {
        todo[i] = 1;
        is_zero[i] = 0;
        ancestor_is_zero[i] = 0;
    }

    for (int i = nnode - 1; i >= 0; i--)
        if (todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int myfam = code[i];
            double mydelta[d];
            int mynvec;
            double myvectors[d * d];
            double myrhs[d];
            int myidx[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                mydelta[k] = delta[j];
                todo[j] = 0;
                myidx[k] = j;
            }
            astfam_constancy(&myfam, &d, mydelta, &mynvec, myvectors, myrhs);
            // at this point mynvec rows of myvectors and the corresponding
            // elements of myrhs are the equality constraints on the canonical
            // statistic of this dependence group
            //
            // how we turn this into directions of constancy depend on
            // whether we want conditional or unconditional (indicated
            // by want_phi) and whether the predecessor node is initial
            for (int j = 0; j < mynvec; j++) {
                double myrhsj = myrhs[j];
                for (int k = 0; k < d; k++) {
                    double foo = myvectors[j + k * d];
                    if (foo != 0.0) {
                        if (doit) {
                            if (knonzero >= nnonzero)
                                error("nnonzero too small");
                            ix[knonzero] = nrow + 1;
                            jx[knonzero] = myidx[k] + 1;
                            xx[knonzero] = foo;
                        }
                        knonzero++;

                        if (myrhsj == 0.0) {
                            int nn = 0;
                            for (int kk = 0; kk < d; kk++) {
                                double foofoo = myvectors[j + kk * d];
                                nn += (foofoo != 0.0);
                            }
                            if (nn == 1)
                                is_zero[myidx[k]] = 1;
                        }
                    }
                }
                if (want_phi) {
                    // jpred is 1-origin indexing
                    int jpred = pred[myidx[0]];
                    if (jpred > 0) {
                        // predecessor is not initial
                        double foo = (- myrhsj);
                        if (foo != 0.0) {
                            if (doit) {
                                if (knonzero >= nnonzero)
                                    error("nnonzero too small");
                                ix[knonzero] = nrow + 1;
                                jx[knonzero] = jpred;
                                xx[knonzero] = foo;
                            }
                            knonzero++;
                        }
                    }
                }
                nrow++;
            }
        }

    for (int i = 0; i < nnode; i++)
        if (! is_zero[i]) {
            int predi = pred[i];
            if (predi > 0) {
                // predecessor not initial
                // predi is 1-origin index
                if (is_zero[predi - 1] || ancestor_is_zero[predi - 1]) {
                    ancestor_is_zero[i] = 1;
                    if (doit) {
                        if (knonzero >= nnonzero)
                            error("nnonzero too small");
                        ix[knonzero] = nrow + 1;
                        jx[knonzero] = i + 1;
                        xx[knonzero] = 1.0;
                    }
                    knonzero++;
                    nrow++;
                }
            }
        }

    return knonzero;
}

SEXP aster_constancy(SEXP pred, SEXP group, SEXP code, SEXP delta,
    SEXP isTheta)
{
    if (! isVectorAtomic(pred))
        error("pred must be atomic");
    if (! isInteger(pred))
        error("pred must be integer");
    if (! isVectorAtomic(group))
        error("group must be atomic");
    if (! isInteger(group))
        error("group must be integer");
    if (! isVectorAtomic(code))
        error("code must be atomic");
    if (! isInteger(code))
        error("code must be integer");
    if (! isVectorAtomic(delta))
        error("delta must be atomic");
    if (! isReal(delta))
        error("delta must be double");
    if (! isVectorAtomic(isTheta))
        error("isTheta must be atomic");
    if (! isLogical(isTheta))
        error("isTheta must be logical");

    int nnode = length(pred);
    if (length(group) != nnode)
        error("lengths of pred and group must match");
    if (length(code) != nnode)
        error("lengths of pred and code must match");
    if (length(delta) != nnode)
        error("lengths of pred and delta must match");
    if (length(isTheta) != 1)
        error("isTheta must be scalar");

    int *pred_ptr = INTEGER(pred);
    int *group_ptr = INTEGER(group);
    int *code_ptr = INTEGER(code);
    double *delta_ptr = REAL(delta);
    _Bool want_phi = (! LOGICAL(isTheta)[0]);

    int nnonzero = get_elements(nnode, pred_ptr, group_ptr, code_ptr,
        delta_ptr, want_phi, 0, 0, NULL, NULL, NULL);

    SEXP result, resultnames;
    PROTECT(result = allocVector(VECSXP, 3));
    PROTECT(resultnames = allocVector(STRSXP, 3));
    SET_STRING_ELT(resultnames, 0, mkChar("i"));
    SET_STRING_ELT(resultnames, 1, mkChar("j"));
    SET_STRING_ELT(resultnames, 2, mkChar("x"));
    namesgets(result, resultnames);
    SEXP foo, bar, baz;
    PROTECT(foo = allocVector(INTSXP, nnonzero));
    PROTECT(bar = allocVector(INTSXP, nnonzero));
    PROTECT(baz = allocVector(REALSXP, nnonzero));
    SET_VECTOR_ELT(result, 0, foo);
    SET_VECTOR_ELT(result, 1, bar);
    SET_VECTOR_ELT(result, 2, baz);

    get_elements(nnode, pred_ptr, group_ptr, code_ptr, delta_ptr, want_phi,
        1, nnonzero, INTEGER(foo), INTEGER(bar), REAL(baz));

    UNPROTECT(5);
    return result;
}

