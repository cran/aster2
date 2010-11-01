
#include "aster.h"
#include "families.h"
#include <stddef.h>
#include <R.h>

void aster_theta_to_phi(int *nnode, int *deriv, int *pred, int *group,
    int *code, double *delta, double *theta, double *dtheta, double *phi,
    double *dphi)
{
    int n = nnode[0];
    int myderiv = deriv[0];
    _Bool todo[n];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++) {
        phi[i] = theta[i];
        if (myderiv == 1)
            dphi[i] = dtheta[i];
        todo[i] = 1;
    }

    for (int i = n - 1; i >= 0; i--)
        if (todo[i] && pred[i] != 0) {
            // evaluate zeroth and first derivative of cumulant function
            // for group containing i, subtract zeroth from phi[pred[i] - 1],
            // and subtract inner product of first and dtheta
            // from dphi[pred[i] - 1]
            // (conversion from 1-origin to 0-origin indexing)
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int myfam = code[i];
            double mytheta[d];
            double mydtheta[d];
            double mydelta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                mytheta[k] = theta[j];
                if (myderiv == 1)
                    mydtheta[k] = dtheta[j];
                mydelta[k] = delta[j];
                todo[j] = 0;
            }
            double zeroth;
            double first[d];
            astfam_cumulant(mytheta, &myfam, &myderiv, mydelta, &zeroth,
                first, NULL, NULL);
            phi[pred[i] - 1] -= zeroth;
            if (myderiv == 1) {
                double inner = 0.0;
                for (int k = 0; k < d; k++)
                    inner += first[k] * mydtheta[k];
                dphi[pred[i] - 1] -= inner;
            }
        }
}

void aster_phi_to_theta(int *nnode, int *deriv, int *pred, int *group,
    int *code, double *delta, double *phi, double *dphi, double *theta,
    double *dtheta)
{
    int n = nnode[0];
    int myderiv = deriv[0];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++) {
        theta[i] = phi[i];
        if (myderiv == 1)
            dtheta[i] = dphi[i];
    }

    // this one is trickier than theta-to-phi, cannot do dependence
    // group until all its successors are done, i. e., must wait until
    // last one on the list is visited, so must construct reverse links

    int revlink[n];
    for (int i = 0; i < n; i++)
        revlink[i] = 0;
    for (int i = 0; i < n; i++)
        if (group[i] != 0) {
            int one_origin_i = i + 1;
            int one_origin_j = group[i];
            int j = one_origin_j - 1;
            revlink[j] = one_origin_i;
        }

    // now revlink is just like group, except that links go in opposite
    // direction

    for (int i = n - 1; i >= 0; i--)
        if (pred[i] != 0 && group[i] == 0) {
            // evaluate zeroth and first derivative of cumulant function
            // for group containing i, add zeroth to theta[pred[i] - 1],
            // and add inner product of first and dphi to dtheta[pred[i] - 1]
            // (conversion from 1-origin to 0-origin indexing)
            int d = 0;
            for (int j = i; j >= 0; j = revlink[j] - 1)
                d++;
            int myfam = code[i];
            double mytheta[d];
            double mydtheta[d];
            double mydelta[d];
            for (int j = i, k = 0; j >= 0; j = revlink[j] - 1, k++) {
                mytheta[k] = theta[j];
                if (myderiv == 1)
                    mydtheta[k] = dtheta[j];
                mydelta[k] = delta[j];
            }
            double zeroth;
            double first[d];
            astfam_cumulant(mytheta, &myfam, &myderiv, mydelta, &zeroth,
                first, NULL, NULL);
            theta[pred[i] - 1] += zeroth;
            if (myderiv == 1) {
                double inner = 0.0;
                for (int k = 0; k < d; k++)
                    inner += first[k] * mydtheta[k];
                dtheta[pred[i] - 1] += inner;
            }
        }
}

void aster_theta_to_xi(int *nnode, int *deriv, int *group,
    int *code, double *delta, double *theta, double *dtheta, double *xi,
    double *dxi)
{
    int n = nnode[0];
    int myderiv = deriv[0];
    _Bool todo[n];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++)
        todo[i] = 1;

    for (int i = n - 1; i >= 0; i--)
        if (todo[i]) {
            // evaluate zeroth, first, and second derivatives
            // of cumulant function for group containing i.
            // first derivative is xi for group, matrix multiply
            // of second derivative times dtheta is dxi
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int myfam = code[i];
            double mytheta[d];
            double mydtheta[d];
            double mydelta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                mytheta[k] = theta[j];
                if (myderiv == 1)
                    mydtheta[k] = dtheta[j];
                mydelta[k] = delta[j];
                todo[j] = 0;
            }
            int thederiv = myderiv + 1;
            double zeroth;
            double first[d];
            double second[d * d];
            astfam_cumulant(mytheta, &myfam, &thederiv, mydelta, &zeroth,
                first, second, NULL);
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                xi[j] = first[k];
                if (myderiv == 1) {
                    double inner = 0.0;
                    for (int l = 0; l < d; l++)
                        inner += second[k * d + l] * mydtheta[l];
                    dxi[j] = inner;
                }
            }
        }
}

void aster_xi_to_theta(int *nnode, int *deriv, int *group,
    int *code, double *delta, double *xi, double *dxi, double *theta,
    double *dtheta)
{
    int n = nnode[0];
    int myderiv = deriv[0];
    _Bool todo[n];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++)
        todo[i] = 1;

    for (int i = n - 1; i >= 0; i--)
        if (todo[i]) {
            // evaluate zeroth and first derivatives
            // of link function for group containing i.
            // zeroth derivative is theta for group, matrix multiply
            // of second derivative times dxi is dtheta
            int d = 0;
            for (int j = i; j >= 0; j = group[j] - 1)
                d++;
            int myfam = code[i];
            double myxi[d];
            double mydxi[d];
            double mydelta[d];
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                myxi[k] = xi[j];
                if (myderiv == 1)
                    mydxi[k] = dxi[j];
                mydelta[k] = delta[j];
                todo[j] = 0;
            }
            double zeroth[d];
            double first[d * d];
            astfam_link(myxi, &myfam, &myderiv, mydelta, zeroth, first);
            for (int j = i, k = d - 1; j >= 0; j = group[j] - 1, k--) {
                theta[j] = zeroth[k];
                if (myderiv == 1) {
                    double inner = 0.0;
                    for (int l = 0; l < d; l++)
                        inner += first[k + d * l] * mydxi[l];
                    dtheta[j] = inner;
                }
            }
        }
}

void aster_xi_to_mu(int *nnode, int *deriv, int *pred, double *initial,
    double *xi, double *dxi, double *mu, double *dmu)
{
    int n = nnode[0];
    int myderiv = deriv[0];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++) {
        double mupred = initial[i];
        if (pred[i] != 0) {
            // conversion from 1-origin to 0-origin indexing
            mupred = mu[pred[i] - 1];
        }
        mu[i] = xi[i] * mupred;
    }

    if (myderiv == 0)
        return;

    for (int i = 0; i < n; i++)
        dmu[i] = 0.0;

    for (int i = 0; i < n; i++)
        for (int j = i; j >= 0; j = pred[j] - 1) {
            // compute term for dmu[i] involving partial of mu[i] with
            // respect to xi[j]
            double prod = 1.0;
            for (int k = i; k >= 0; k = pred[k] - 1) {
                 if (k != j)
                     prod *= xi[k];
                 else
                     prod *= dxi[k];
                 if (pred[k] == 0)
                     prod *= initial[k];
            }
            dmu[i] += prod;
        }
}

void aster_mu_to_xi(int *nnode, int *deriv, int *pred, double *initial,
    double *mu, double *dmu, double *xi, double *dxi)
{
    int n = nnode[0];
    int myderiv = deriv[0];

    if (! (myderiv == 0 || myderiv == 1))
        error("deriv must be zero or one");

    for (int i = 0; i < n; i++) {
        double mupred = initial[i];
        if (pred[i] != 0) {
            // conversion from 1-origin to 0-origin indexing
            mupred = mu[pred[i] - 1];
        }
        if (mupred <= 0.0)
            error("mu for predecessor nonpositive");
        xi[i] = mu[i] / mupred;
    }

    if (myderiv == 0)
        return;

    for (int i = 0; i < n; i++)
        dxi[i] = 0.0;

    for (int i = 0; i < n; i++) {
        // conversion from 1-origin to 0-origin indexing
        int j = pred[i] - 1;

        double mupred = initial[i];
        if (j >= 0)
            mupred = mu[j];

        dxi[i] += 1 / mupred * dmu[i];
        if (j >= 0)
            dxi[i] -= mu[i] / (mupred * mupred) * dmu[j];
    }
}

