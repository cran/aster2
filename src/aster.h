
#ifndef ASTER2_ASTER_H_
#define ASTER2_ASTER_H_

void aster_validate(int *nnode, double *resp, int *pred, int *group,
    int *code, double *initial, double *delta);

void aster_validate_theta(int *nnode, int *group, int *code,
    double *delta, double *theta);

void aster_validate_xi(int *nnode, int *group, int *code,
    double *delta, double *xi);

void aster_theta_to_phi(int *nnode, int *deriv, int *pred, int *group,
    int *code, double *delta, double *theta, double *dtheta, double *phi,
    double *dphi);

void aster_phi_to_theta(int *nnode, int *deriv, int *pred, int *group,
    int *code, double *delta, double *phi, double *dphi, double *theta,
    double *dtheta);

void aster_theta_to_xi(int *nnode, int *deriv, int *group,
    int *code, double *delta, double *theta, double *dtheta, double *xi,
    double *dxi);

void aster_xi_to_theta(int *nnode, int *deriv, int *group,
    int *code, double *delta, double *xi, double *dxi, double *theta,
    double *dtheta);

void aster_xi_to_mu(int *nnode, int *deriv, int *pred, double *initial,
    double *xi, double *dxi, double *mu, double *dmu);

void aster_mu_to_xi(int *nnode, int *deriv, int *pred, double *initial,
    double *mu, double *dmu, double *xi, double *dxi);

void aster_starting_theta(int *nnode, int *group, int *code, double *theta);

#endif /* ASTER2_ASTER_H_ */

