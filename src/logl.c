
#include "aster.h"
#include "families.h"

void aster_starting_theta(int *nnode, int *group, int *code, double *theta)
{
    int n = nnode[0];
    _Bool todo[n];

    for (int i = 0; i < n; i++)
        todo[i] = 1;

    for (int i = n - 1; i >= 0; i--)
        if (todo[i]) {
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int myfam = code[i];
            double mytheta[d];
            astfam_start_theta(&myfam, &d, mytheta);
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                theta[j] = mytheta[k];
                todo[j] = 0;
            }
        }
}

