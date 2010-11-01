
#include "constancy.h"
#include "families.h"

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

    _Bool todo[nnode];
    for (int i = 0; i < nnode; i++)
        todo[i] = 1;

    int nvec_total = 0;
    for (int i = nnode - 1; i >= 0; i--)
        if (todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group_ptr[j] - 1)
                d++;
            int myfam = code_ptr[i];
            double mydelta[d];
            int mynvec;
            double myvectors[d * d];
            double myrhs[d];
            for (int j = i, k = d - 1; j >= 0; j = group_ptr[j] - 1, k--) {
                mydelta[k] = delta_ptr[j];
                todo[j] = 0;
            }
            astfam_constancy(&myfam, &d, mydelta, &mynvec, myvectors, myrhs);
            nvec_total += mynvec;
        }

    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, nvec_total, nnode));
    double *result_ptr = REAL(result);
    for (int i = 0; i < nvec_total * nnode; i++)
        result_ptr[i] = 0.0;

    for (int i = 0; i < nnode; i++)
        todo[i] = 1;

    int kvec = 0;
    for (int i = nnode - 1; i >= 0; i--)
        if (todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group_ptr[j] - 1)
                d++;
            int myfam = code_ptr[i];
            double mydelta[d];
            int mynvec;
            double myvectors[d * d];
            double myrhs[d];
            int myidx[d];
            for (int j = i, k = d - 1; j >= 0; j = group_ptr[j] - 1, k--) {
                mydelta[k] = delta_ptr[j];
                todo[j] = 0;
                myidx[k] = j;
            }
            astfam_constancy(&myfam, &d, mydelta, &mynvec, myvectors, myrhs);
            for (int j = 0; j < mynvec; j++) {
                for (int k = 0; k < d; k++) {
                    int result_idx = (nvec_total - 1 - kvec) + j +
                        nvec_total * myidx[k];
                    result_ptr[result_idx] = myvectors[j + d * k];
                }
                if (want_phi) {
                    int jpred = pred_ptr[myidx[0]] - 1;
                    if (jpred >= 0) {
                        int result_idx = (nvec_total - 1 - kvec) + j +
                            nvec_total * jpred;
                        result_ptr[result_idx] = (- myrhs[j]);
                    }
                }
            }
            kvec += mynvec;
        }

    UNPROTECT(1);
    return result;
}

