
#include <R.h>
#include "aster.h"
#include "families.h"

void aster_predecessor_zero_cond(int *nnode, int *pred, double *resp,
    int *result)
{
    int n = nnode[0];
    // note: initial cannot be zero (checked in R function
    // validasterdata) so ipred == 0 implies predecessor nonzero
    for (int i = 0; i < n; i++)
        result[i] = pred[i] > 0 && resp[pred[i] - 1] == 0;
}

void aster_predecessor_zero_unco(int *nnode, int *pred, int *group, int *code,
    double *delta, int *result)
{
    int n = nnode[0];
    _Bool *todo = (_Bool *) R_alloc(n, sizeof(_Bool));
    int *revgroup = (int *) R_alloc(n, sizeof(int));

    aster_revlink(nnode, group, revgroup);

    for (int i = 0; i < n; i++) {
        todo[i] = 1;
        result[i] = 0;
    }

    for (int i = 0; i < n; i++)
        if (todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = revgroup[j] - 1) {
                d++;
                todo[j] = 0;
            }
            int ipred = pred[i];
            // note: initial cannot be zero (checked in R function
            // validasterdata) so ipred == 0 implies predecessor nonzero
            if (ipred == 0 || result[ipred - 1] == 0) {
                int myfam = code[i];
                double mydelta[d];
                int myzeros[d];
                int myidx[d];
                for (int j = i, k = 0; j >= 0; j = revgroup[j] - 1, k++) {
                    mydelta[k] = delta[j];
                    myidx[k] = j;
                }
                astfam_is_zero(&myfam, &d, mydelta, myzeros);
                for (int k = 0; k < d; k++)
                    if (myzeros[k])
                        result[myidx[k]] = 2;
            } else {
                for (int j = i; j >= 0; j = revgroup[j] - 1)
                        result[j] = 1;
            }
        }
    for (int i = 0; i < n; i++)
        if (result[i] == 2)
            result[i] = 0;
}

void aster_revlink(int *nnode, int *group, int *revgroup)
{
    int n = nnode[0];

    for (int i = 0; i < n; i++)
        revgroup[i] = 0;

    for (int i = 0; i < n; i++)
        if (group[i] != 0) {
            int one_origin_i = i + 1;
            int one_origin_j = group[i];
            int j = one_origin_j - 1;
            revgroup[j] = one_origin_i;
        }
}

