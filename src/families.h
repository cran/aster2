
#ifndef ASTER2_FAMILIES_H_
#define ASTER2_FAMILIES_H_

void astfam_cumulant(double *theta, int *fam, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third);

void astfam_link(double *xi, int *fam, int *deriv, double *delta,
    double *zeroth, double *first);

void astfam_dimension(int *fam, int *dimen);

void astfam_nfam(int *nfam);

void astfam_clear(void);

void astfam_set(char **name, double *hyper1, double *hyper2);

void astfam_set_tolerance(double *tolerance);

void astfam_reset_tolerance(void);

void astfam_validate_pred(int *fam, double *ypred);

void astfam_validate_delta(int *fam, int *dimen, double *delta);

void astfam_validate_resp(int *fam, int *dimen, double *delta, double *ypred,
    double *resp);

void astfam_validate_theta(int *fam, int *dimen, double *delta, double *theta);

void astfam_validate_xi(int *fam, int *dimen, double *delta, double *xi);

void astfam_constancy(int *fam, int *dimen, double *delta, int *nvec,
    double *vectors, double *rhs);

void astfam_start_theta(int *fam, int *dimen, double *theta);

void astfam_is_zero(int *fam, int *dimen, double *delta, int *zeros);

#endif /* ASTER2_FAMILIES_H_ */

