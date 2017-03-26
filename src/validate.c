
#include <R.h>
#include "aster.h"
#include "families.h"

// validate the data for an aster model, validate everything, not
// just what is not validated in the R function validasterdata, since
// this function might be called from some other C function

void aster_validate(int *nnode, double *resp, int *pred, int *group,
    int *code, double *initial, double *delta)
{
    int n = nnode[0];
    if (n < 1)
        error("graph must have at least one node");

    for (int i = 0; i < n; i++) {
        if (pred[i] < 0)
            error("repred must be nonnegative");
        if (group[i] < 0)
            error("regroup must be nonnegative");
        if (pred[i] >= i + 1)
            error("must have repred[j] < j for all j");
        if (group[i] >= i + 1)
            error("must have regroup[j] < j for all j");
    }

    int nfam;
    astfam_nfam(&nfam);
    if (nfam == 0)
        error("no families set");
    for (int i = 0; i < n; i++)
        if (code[i] <= 0 || code[i] > nfam)
            error("invalid recode, doesn't match number of set family");

    // check all nodes in group have same predecessor and same family
    for (int i = 0; i < n; i++)
        if (group[i] != 0) {
            int j = group[i] - 1; // convert to 0-origin index
            if (pred[i] != pred[j])
                error("all nodes in group must have same predecessor");
            if (code[i] != code[j])
                error("all nodes in group must have same family");
        }

    _Bool is_done[n];
    for (int i = 0; i < n; i++)
        is_done[i] = 0;

    for (int i = n - 1; i >= 0; i--)
        if (! is_done[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int correct_dimension;
            astfam_dimension(&code[i], &correct_dimension);
            if (d != correct_dimension)
                error("size of group does not match dimension of family");
            double ypred = pred[i] > 0 ? resp[pred[i] - 1] : initial[i];
            double my_resp[d];
            double my_delta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                my_resp[k] = resp[j];
                my_delta[k] = delta[j];
                is_done[j] = 1;
            }
            int my_fam = code[i];
            astfam_validate_pred(&my_fam, &ypred);
            astfam_validate_delta(&my_fam, &d, my_delta);
            astfam_validate_resp(&my_fam, &d, my_delta, &ypred, my_resp);
        }
}

// validate conditional canonical parameter vector for an aster model.
// Don't validate everything.  Call aster_validate for that.

void aster_validate_theta(int *nnode, int *pred, int *group, int *code,
    int *want_uam, double *resp, double *delta, double *theta)
{
    int n = nnode[0];
    _Bool *todo = (_Bool *) R_alloc(n, sizeof(_Bool));
    int *zeropred = (int *) R_alloc(n, sizeof(int));

    if (want_uam[0])
        aster_predecessor_zero_unco(nnode, pred, group, code, delta, zeropred);
    else
        aster_predecessor_zero_cond(nnode, pred, resp, zeropred);

    for (int i = 0; i < n; i++)
        todo[i] = 1;

    for (int i = n - 1; i >= 0; i--)
        if ((! zeropred[i]) && todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            double my_theta[d];
            double my_delta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                my_theta[k] = theta[j];
                my_delta[k] = delta[j];
                todo[j] = 0;
            }
            int my_fam = code[i];
            astfam_validate_theta(&my_fam, &d, my_delta, my_theta);
        }
}

// validate conditional mean value parameter vector for an aster model.
// Don't validate everything.  Call aster_validate for that.

void aster_validate_xi(int *nnode, int *pred, int *group, int *code,
    int *want_uam, double *resp, double *delta, double *xi)
{
    int n = nnode[0];
    _Bool *todo = (_Bool *) R_alloc(n, sizeof(_Bool));
    int *zeropred = (int *) R_alloc(n, sizeof(int));

    if (want_uam[0])
        aster_predecessor_zero_unco(nnode, pred, group, code, delta, zeropred);
    else
        aster_predecessor_zero_cond(nnode, pred, resp, zeropred);

    for (int i = 0; i < n; i++)
        todo[i] = 1;

    for (int i = n - 1; i >= 0; i--)
        if ((! zeropred[i]) && todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            double my_xi[d];
            double my_delta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                my_xi[k] = xi[j];
                my_delta[k] = delta[j];
                todo[j] = 0;
            }
            int my_fam = code[i];
            astfam_validate_xi(&my_fam, &d, my_delta, my_xi);
        }
}

