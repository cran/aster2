
#include <R.h>
#include <Rmath.h>
#include "families.h"
#include <string.h>
#include <math.h>
#include <float.h>

static const double default_tolerance = 8 * DBL_EPSILON;
static double tolerance = 8 * DBL_EPSILON;

typedef struct Family Family_t;

struct Family {
    char *name;
    double hyper1;
    double hyper2;
    int dimension;
    void (*cumulant)(double *theta, int *deriv, double *delta, double *zeroth,
        double *first, double *second, double *third, Family_t *that);
    void (*link)(double *xi, int *deriv, double *delta, double *zeroth,
        double *first, Family_t *that);
    void (*validate_pred)(double *ypred);
    void (*validate_delta)(int d, double *delta);
    void (*validate_resp)(double ypred, int d, double *delta, double *resp);
    void (*validate_theta)(int d, double *delta, double *theta);
    void (*validate_xi)(int d, double *delta, double *xi);
    void (*constancy)(int d, double *delta, int *nvec, double *vectors,
        double *rhs);
    void (*start_theta)(int d, double *theta);
};

// Bernoulli

#ifndef __GNUC__
static void bernoulli_cumulant(double *theta, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third,
    Family_t *that)
#else
static void bernoulli_cumulant(double *theta, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third,
    Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! R_finite(*theta))
        error("theta must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 3)
        error("deriv must be 0, 1, 2, or 3");
    if (*delta < 0.0) {
        *zeroth = 0.0;
        if (*deriv >= 1) *first = 0.0;
        if (*deriv >= 2) *second = 0.0;
        if (*deriv >= 3) *third = 0.0;
    }
    if (*delta > 0.0) {
        *zeroth = *theta;
        if (*deriv >= 1) *first = 1.0;
        if (*deriv >= 2) *second = 0.0;
        if (*deriv >= 3) *third = 0.0;
    }
    if (*delta == 0.0) {
        if (*theta > 0.0)
            *zeroth = *theta + log1p(exp(- *theta));
        else
            *zeroth = log1p(exp(*theta));
        if (*deriv >= 1) {
            double p = 0.0;
            double q = 0.0;
            if (*theta > 0.0) {
                double foo = exp(- *theta);
                p = 1.0 / (1.0 + foo);
                q = foo / (1.0 + foo);
            } else {
                double foo = exp(*theta);
                p = foo / (1.0 + foo);
                q = 1.0 / (1.0 + foo);
            }
            *first = p;
            if (*deriv >= 2) *second = p * q;
            if (*deriv >= 3) *third = - p * q * tanh(*theta / 2.0);
        }
    }
}

#ifndef __GNUC__
static void bernoulli_link(double *xi, int *deriv, double *delta,
    double *zeroth, double *first, Family_t *that)
#else
static void bernoulli_link(double *xi, int *deriv, double *delta,
    double *zeroth, double *first, Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! R_finite(*xi))
        error("xi must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 1)
        error("deriv must be 0 or 1");
    if (*delta != 0.0) {
        *zeroth = 0.0;
        if (*deriv >= 1)
            *first = 0.0;
    }
    if (*delta == 0.0) {
        *zeroth = log(*xi) - log1p(- *xi);
        if (*deriv >= 1)
            *first = 1.0 / (*xi * (1 - *xi));
    }
}

static void bernoulli_validate_pred(double *ypred)
{
    double x = ypred[0];
    if (x != trunc(x))
        error("Bernoulli predecessor noninteger");
    if (x < 0.0)
        error("Bernoulli predecessor negative");
}

#ifndef __GNUC__
static void do_nothing_validate_delta(int d, double *delta)
#else
static void do_nothing_validate_delta(int d __attribute__ ((unused)),
    double *delta __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    return;
}

#ifndef __GNUC__
static void bernoulli_validate_resp(double ypred, int d, double *delta,
    double *resp)
#else
static void bernoulli_validate_resp(double ypred,
    int d __attribute__ ((unused)), double *delta, double *resp)
#endif /* __GNUC__ */
{
    double x = resp[0];
    if (ypred == 0.0 && x != 0.0)
        error("predecessor zero but successor nonzero");
    if (x != trunc(x))
        error("Bernoulli successor noninteger");
    if (*delta == 0.0 && (x < 0.0 || x > ypred))
        error("Bernoulli successor not between 0 and predecessor");
    if (*delta > 0.0 && x != ypred)
        error("Bernoulli successor != predecessor for upper limit");
    if (*delta < 0.0 && x != 0.0)
        error("Bernoulli successor != 0 for lower limit");
}

#ifndef __GNUC__
static void do_nothing_validate_theta(int d, double *delta, double *theta)
#else
static void do_nothing_validate_theta(int d __attribute__ ((unused)),
    double *delta __attribute__ ((unused)),
    double *theta __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    return;
}

#ifndef __GNUC__
static void bernoulli_validate_xi(int d, double *delta, double *xi)
#else
static void bernoulli_validate_xi(int d __attribute__ ((unused)),
    double *delta, double *xi)
#endif /* __GNUC__ */
{
    if (*delta == 0.0 && (*xi <= 0.0 || *xi >= 1.0))
        error("Bernoulli xi not strictly between 0 and 1");
    if (*delta > 0.0 && *xi != 1.0)
        error("Bernoulli xi != 1 for upper limit");
    if (*delta < 0.0 && *xi != 0.0)
        error("Bernoulli xi != 0 for lower limit");
}

#ifndef __GNUC__
static void bernoulli_constancy(int d, double *delta, int *nvec,
    double *vectors, double *rhs)
#else
static void bernoulli_constancy(int d __attribute__ ((unused)),
    double *delta, int *nvec, double *vectors, double *rhs)
#endif /* __GNUC__ */
{
    if (*delta == 0.0) {
        *nvec = 0;
    }
    if (*delta < 0.0) {
        *nvec = 1;
        *vectors = 1.0;
        *rhs = 0.0;
    }
    if (*delta > 0.0) {
        *nvec = 1;
        *vectors = 1.0;
        *rhs = 1.0;
    }
}

static void zero_start_theta(int d, double *theta)
{
    for (int i = 0; i < d; i++)
        theta[i] = 0.0;
}

#ifndef __GNUC__
static Family_t bernoulli_constructor(double *hyper1, double *hyper2)
#else
static Family_t bernoulli_constructor(double *hyper1 __attribute__ ((unused)),
    double *hyper2 __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    Family_t result = {
        "bernoulli",
        R_NaN,
        R_NaN,
        1,
        bernoulli_cumulant,
        bernoulli_link,
        bernoulli_validate_pred,
        do_nothing_validate_delta,
        bernoulli_validate_resp,
        do_nothing_validate_theta,
        bernoulli_validate_xi,
        bernoulli_constancy,
        zero_start_theta
    };
    return result;
}

// Poisson

#ifndef __GNUC__
static void poisson_cumulant(double *theta, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third,
    Family_t *that)
#else
static void poisson_cumulant(double *theta, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third,
    Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! R_finite(*theta))
        error("theta must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 3)
        error("deriv must be 0, 1, 2, or 3");
    if (*delta < 0.0) {
        *zeroth = 0.0;
        if (*deriv >= 1) *first = 0.0;
        if (*deriv >= 2) *second = 0.0;
        if (*deriv >= 3) *third = 0.0;
    }
    if (*delta > 0.0)
        error("delta must nonpositive, no upper limit for Poisson");
    if (*delta == 0.0) {
        double foo = exp(*theta);
        *zeroth = foo;
        if (*deriv >= 1) *first = foo;
        if (*deriv >= 2) *second = foo;
        if (*deriv >= 3) *third = foo;
    }
}

#ifndef __GNUC__
static void poisson_link(double *xi, int *deriv, double *delta,
    double *zeroth, double *first, Family_t *that)
#else
static void poisson_link(double *xi, int *deriv, double *delta,
    double *zeroth, double *first, Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! R_finite(*xi))
        error("xi must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 1)
        error("deriv must be 0 or 1");
    if (*delta < 0.0) {
        *zeroth = 0.0;
        if (*deriv >= 1)
            *first = 0.0;
    }
    if (*delta > 0.0)
        error("delta must nonpositive, no upper limit for Poisson");
    if (*delta == 0.0) {
        *zeroth = log(*xi);
        if (*deriv >= 1)
            *first = 1.0 / *xi;
    }
}

static void poisson_validate_pred(double *ypred)
{
    double x = ypred[0];
    if (x < 0.0)
        error("Poisson predecessor negative");
}

#ifndef __GNUC__
static void poisson_validate_delta(int d, double *delta)
#else
static void poisson_validate_delta(int d __attribute__ ((unused)),
    double *delta)
#endif /* __GNUC__ */
{
    if (*delta > 0.0)
       error("delta > 0 not allowed for Poisson");
}

#ifndef __GNUC__
static void poisson_validate_resp(double ypred, int d, double *delta,
    double *resp)
#else
static void poisson_validate_resp(double ypred, int d __attribute__ ((unused)),
    double *delta, double *resp)
#endif /* __GNUC__ */
{
    double x = resp[0];
    if (ypred == 0.0 && x != 0.0)
        error("successor nonzero when predecessor zero");
    if (x != trunc(x))
        error("Poisson successor noninteger");
    if (*delta == 0.0 && x < 0.0)
        error("Poisson successor negative");
    if (*delta < 0.0 && x != 0.0)
        error("Poisson successor != 0 for lower limit");
}

#ifndef __GNUC__
static void poisson_validate_xi(int d, double *delta, double *xi)
#else
static void poisson_validate_xi(int d __attribute__ ((unused)), double *delta,
    double *xi)
#endif /* __GNUC__ */
{
    if (*delta == 0.0 && *xi <= 0.0)
        error("Poisson xi not strictly positive");
    if (*delta < 0.0 && *xi != 0.0)
        error("Poisson xi != 0 for lower limit");
}

#ifndef __GNUC__
static void poisson_constancy(int d, double *delta, int *nvec,
    double *vectors, double *rhs)
#else
static void poisson_constancy(int d __attribute__ ((unused)),
    double *delta, int *nvec, double *vectors, double *rhs)
#endif /* __GNUC__ */
{
    if (*delta == 0.0) {
        *nvec = 0;
    }
    if (*delta < 0.0) {
        *nvec = 1;
        *vectors = 1.0;
        *rhs = 0.0;
    }
}

#ifndef __GNUC__
static Family_t poisson_constructor(double *hyper1, double *hyper2)
#else
static Family_t poisson_constructor(double *hyper1 __attribute__ ((unused)),
    double *hyper2 __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    Family_t result = {
        "poisson",
        R_NaN,
        R_NaN,
        1,
        poisson_cumulant,
        poisson_link,
        poisson_validate_pred,
        poisson_validate_delta,
        poisson_validate_resp,
        do_nothing_validate_theta,
        poisson_validate_xi,
        poisson_constancy,
        zero_start_theta
    };
    return result;
}

// zero-truncated Poisson

#ifndef __GNUC__
static void zero_truncated_poisson_cumulant(double *theta, int *deriv,
    double *delta, double *zeroth, double *first, double *second,
    double *third, Family_t *that)
#else
static void zero_truncated_poisson_cumulant(double *theta, int *deriv,
    double *delta, double *zeroth, double *first, double *second,
    double *third, Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! R_finite(*theta))
        error("theta must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 3)
        error("deriv must be 0, 1, 2, or 3");
    if (*delta < 0.0) {
        *zeroth = *theta;
        if (*deriv >= 1) *first = 1.0;
        if (*deriv >= 2) *second = 0.0;
        if (*deriv >= 3) *third = 0.0;
    }
    if (*delta > 0.0)
        error("delta must nonpositive, no upper limit"
            " for zero-truncated Poisson");
    if (*delta == 0.0) {
        double m = exp(*theta);
        double mu = R_NaN;
        if (*theta > (- 4.0)) {
            *zeroth = m + log1p(- exp(- m));
            if (*deriv >= 1) mu = m / (1.0 - exp(- m));
        } else /* (*theta) <= (- 4.0) */ {
             double bar = m / 2.0 * (1.0 + m / 3.0 * (1.0 + m / 4.0 *
                 (1.0 + m / 5.0 * (1.0 + m / 6.0 * (1.0 + m / 7.0 *
                 (1.0 + m / 8.0))))));
             *zeroth = *theta + log1p(bar);
             if (*deriv >= 1) mu = (m + 1.0 / (1.0 + bar));
        }
        if (*deriv >= 1) *first = mu;
        if (*deriv >= 2) *second = mu * (1.0 + m - mu);
        if (*deriv >= 3) *third = mu * ((1.0 + m - mu) *
            (1.0 + m - 2.0 * mu) + m);
    }
}

static void zero_truncated_poisson_link(double *xi, int *deriv,
    double *delta, double *zeroth, double *first, Family_t *that)
{
    if (! R_finite(*xi))
        error("xi must be finite");
    if (! R_finite(*delta))
        error("delta must be finite");
    if (*deriv < 0 || *deriv > 1)
        error("deriv must be 0 or 1");
    if (*delta < 0.0) {
        *zeroth = 0.0;
        if (*deriv >= 1)
            *first = 0.0;
    }
    if (*delta > 0.0)
        error("delta must nonpositive, no upper limit"
            " for zero-truncated Poisson");
    if (*delta == 0.0) {
        double xii = *xi;
        double moo = xii > 2.0 ? xii : 2.0 * (xii - 1.0);
        double thetai = log(moo);
        double newton = R_PosInf;
        double newton_save = R_PosInf;
        double newton_save_save = R_PosInf;
        double newton_tolerance = sqrt(DBL_EPSILON);
        double myzeroth;
        double myfirst;
        double mysecond;
        int two = 2;
        while (fabs(newton_save_save) >= newton_tolerance) {
            that->cumulant(&thetai, &two, delta, &myzeroth, &myfirst,
                &mysecond, NULL, that);
            newton_save_save = newton_save;
            newton_save = newton;
            newton = (xii - myfirst) / mysecond;
            thetai = thetai + newton;
        }
        that->cumulant(&thetai, &two, delta, &myzeroth, &myfirst,
            &mysecond, NULL, that);
        *zeroth = thetai;
        *first = 1.0 / mysecond;
    }
}

static void zero_truncated_poisson_validate_pred(double *ypred)
{
    double x = ypred[0];
    if (x != trunc(x))
        error("zero-truncated Poisson predecessor noninteger");
    if (x < 0.0)
        error("zero-truncated Poisson predecessor negative");
}

#ifndef __GNUC__
static void zero_truncated_poisson_validate_delta(int d, double *delta)
#else
static void zero_truncated_poisson_validate_delta(
    int d __attribute__ ((unused)), double *delta)
#endif /* __GNUC__ */
{
    if (*delta > 0.0)
       error("delta > 0 not allowed for zero-truncated Poisson");
}

#ifndef __GNUC__
static void zero_truncated_poisson_validate_resp(double ypred, int d,
    double *delta, double *resp)
#else
static void zero_truncated_poisson_validate_resp(double ypred,
    int d __attribute__ ((unused)), double *delta, double *resp)
#endif /* __GNUC__ */
{
    double x = resp[0];
    if (ypred == 0.0 && x != 0.0)
        error("successor nonzero when predecessor zero");
    if (x != trunc(x))
        error("zero-truncated Poisson successor noninteger");
    if (x < 0.0)
        error("zero-truncated Poisson successor negative");
    if (*delta < 0.0 && x != ypred)
        error("zero-truncated Poisson successor != predecessor"
            " for lower limit");
}

#ifndef __GNUC__
static void zero_truncated_poisson_validate_xi(int d, double *delta,
    double *xi)
#else
static void zero_truncated_poisson_validate_xi(int d __attribute__ ((unused)),
    double *delta, double *xi)
#endif /* __GNUC__ */
{
    if (*delta == 0.0 && *xi <= 1.0)
        error("zero-truncated Poisson xi not strictly greater than 1");
    if (*delta < 0.0 && *xi != 1.0)
        error("zero-truncated Poisson xi != 1 for lower limit");
}

#ifndef __GNUC__
static void zero_truncated_poisson_constancy(int d, double *delta, int *nvec,
    double *vectors, double *rhs)
#else
static void zero_truncated_poisson_constancy(int d __attribute__ ((unused)),
    double *delta, int *nvec, double *vectors, double *rhs)
#endif /* __GNUC__ */
{
    if (*delta == 0.0) {
        *nvec = 0;
    }
    if (*delta < 0.0) {
        *nvec = 1;
        *vectors = 1.0;
        *rhs = 1.0;
    }
}

#ifndef __GNUC__
static Family_t zero_truncated_poisson_constructor(double *hyper1,
    double *hyper2)
#else
static Family_t zero_truncated_poisson_constructor(double *hyper1
    __attribute__ ((unused)), double *hyper2 __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    Family_t result = {
        "zero.truncated.poisson",
        R_NaN,
        R_NaN,
        1,
        zero_truncated_poisson_cumulant,
        zero_truncated_poisson_link,
        zero_truncated_poisson_validate_pred,
        zero_truncated_poisson_validate_delta,
        zero_truncated_poisson_validate_resp,
        do_nothing_validate_theta,
        zero_truncated_poisson_validate_xi,
        zero_truncated_poisson_constancy,
        zero_start_theta
    };
    return result;
}

// Multinomial

static void multinomial_cumulant(double *theta, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third,
    Family_t *that)
{
    int d = (*that).dimension;
    for (int i = 0; i < d; i++) {
        if (! R_finite(theta[i]))
            error("all components of theta must be finite");
        if (! R_finite(delta[i]))
            error("all components of delta must be finite");
    }
    if (*deriv < 0 || *deriv > 3)
        error("deriv must be 0, 1, 2, or 3");
    double max_delta = R_NegInf;
    for (int i = 0; i < d; i++)
        if (delta[i] > max_delta)
            max_delta = delta[i];
    double max_theta = R_NegInf;
    for (int i = 0; i < d; i++)
        if (delta[i] == max_delta && theta[i] > max_theta)
            max_theta = theta[i];
    int nmax_theta = 0;
    for (int i = 0; i < d; i++)
        if (delta[i] == max_delta && theta[i] == max_theta)
            nmax_theta++;
    double mu_unnormalized[d];
    for (int i = 0; i < d; i++)
        if (delta[i] == max_delta) {
            mu_unnormalized[i] = exp(theta[i] - max_theta);
        } else {
            mu_unnormalized[i] = 0.0;
        }
    double sum_less = 0.0;
    for (int i = 0; i < d; i++)
        if (theta[i] < max_theta)
            sum_less += mu_unnormalized[i];
    *zeroth = max_theta + log(nmax_theta) + log1p(sum_less / nmax_theta);
    if (*deriv == 0) return;
    double mu[d];
    for (int i = 0; i < d; i++) {
        mu[i] = mu_unnormalized[i] / (nmax_theta + sum_less);
        first[i] = mu[i];
    }
    if (*deriv <= 1) return;
    for (int i = 0, k = 0; i < d; i++)
        for (int j = 0; j < d; j++, k++) {
            if (i == j)
                second[k] = mu[i] * (1.0 - mu[i]);
            else
                second[k] = - mu[i] * mu[j];
        }
    if (*deriv <= 2) return;
    for (int i = 0, l = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            for (int k = 0; k < d; k++, l++) {
                if (i == j && j == k )
                    third[l] = mu[i] * (1.0 - mu[i]) * (1.0 - 2 * mu[i]);
                if (i == j && j != k)
                    third[l] = - mu[i] * mu[k] * (1.0 - 2 * mu[i]);
                if (i == k && j != k)
                    third[l] = - mu[i] * mu[j] * (1.0 - 2 * mu[i]);
                if (j == k && i != k)
                    third[l] = - mu[j] * mu[i] * (1.0 - 2 * mu[j]);
                if (i != j && j != k && k != i)
                    third[l] = 2.0 * mu[i] * mu[j] * mu[k];
            }
}

static void multinomial_link(double *xi, int *deriv, double *delta,
    double *zeroth, double *first, Family_t *that)
{
    int d = (*that).dimension;
    for (int i = 0; i < d; i++) {
        if (! R_finite(xi[i]))
            error("all components of xi must be finite");
        if (! R_finite(delta[i]))
            error("all components of delta must be finite");
    }
    if (*deriv < 0 || *deriv > 1)
        error("deriv must be 0 or 1");
    double max_delta = R_NegInf;
    for (int i = 0; i < d; i++)
        if (delta[i] > max_delta)
            max_delta = delta[i];
    int izero = 0;
    for (int i = 0; i < d; i++)
        if (delta[i] == max_delta) {
            izero = i;
            break;
        }
    for (int i = 0; i < d; i++) {
        zeroth[i] = 0.0;
        if (delta[i] == max_delta && i != izero)
            zeroth[i] = log(xi[i]) - log(xi[izero]);
    }
    if (*deriv == 0)
        return;
    for (int i = 0; i < d * d; i++)
        first[i] = 0.0;
    for (int i = 0; i < d; i++)
        if (delta[i] == max_delta && i != izero) {
            first[i + i * d] = 1.0 / xi[i];
            first[i + izero * d] = (- 1.0 / xi[izero]);
        }
}

static void multinomial_validate_pred(double *ypred)
{
    double x = ypred[0];
    if (x != trunc(x))
        error("multinomial predecessor noninteger");
    if (x < 0.0)
        error("multinomial predecessor negative");
}

static void multinomial_validate_resp(double ypred, int d, double *delta,
    double *resp)
{
    double sum = 0.0;
    double max_delta = R_NegInf;
    for (int i = 0; i < d; i++) {
        double x = resp[i];
        if (ypred == 0.0 && x != 0.0)
            error("successor nonzero when predecessor zero");
        if (x != trunc(x))
            error("multinomial successor noninteger");
        sum += x;
        if (delta[i] > max_delta)
            max_delta = delta[i];
    }
    if (sum != ypred)
        error("sum of components of multinomial successor != predecessor");
    for (int i = 0; i < d; i++)
        if (delta[i] < max_delta && resp[i] != 0.0)
            error("component of multinomial successor != 0 when"
                " so constrained by limit");
}

static void multinomial_validate_xi(int d, double *delta, double *xi)
{
    for (int i = 0; i < d; i++)
        if (xi[i] < 0)
            error("component of multinomial xi < 0");
    double sum = 0.0;
    double max_delta = R_NegInf;
    for (int i = 0; i < d; i++) {
        sum += xi[i];
        if (delta[i] > max_delta)
            max_delta = delta[i];
    }
    if (fabs(sum - 1.0) > tolerance)
        error("sum of components of multinomial xi != 1");
    for (int i = 0; i < d; i++)
        if (delta[i] < max_delta && xi[i] != 0.0)
            error("component of multinomial xi != 0 when"
                " so constrained by limit");
}

static void multinomial_constancy(int d, double *delta, int *nvec,
    double *vectors, double *rhs)
{
    double max_delta = R_NegInf;
    for (int i = 0; i < d; i++)
        if (delta[i] > max_delta)
            max_delta = delta[i];

    *nvec = 1;
    for (int i = 0; i < d; i++)
        vectors[d * i] = 1.0;
    rhs[0] = 1.0;

    for (int i = 0; i < d; i++)
        if (delta[i] < max_delta) {
            for (int j = 0; j < d; j++)
                vectors[*nvec + d * j] = 0.0;
            vectors[*nvec + d * i] = 1.0;
            rhs[*nvec] = 0.0;
            (*nvec)++;
        }
}

#ifndef __GNUC__
static Family_t multinomial_constructor(double *hyper1, double *hyper2)
#else
static Family_t multinomial_constructor(double *hyper1, double *hyper2
    __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    int d = *hyper1;
    if (! (d == (*hyper1)))
        error("hyperparameter must be integer");
    if (! (d >= 1))
        error("hyperparameter must be positive");
    Family_t result = {
        "multinomial",
        *hyper1,
        R_NaN,
        d,
        multinomial_cumulant,
        multinomial_link,
        multinomial_validate_pred,
        do_nothing_validate_delta,
        multinomial_validate_resp,
        do_nothing_validate_theta,
        multinomial_validate_xi,
        multinomial_constancy,
        zero_start_theta
    };
    return result;
}

// Normal location-scale

#ifndef __GNUC__
static void normal_location_scale_cumulant(double *theta, int *deriv,
    double *delta, double *zeroth, double *first, double *second,
    double *third, Family_t *that)
#else
static void normal_location_scale_cumulant(double *theta, int *deriv,
    double *delta, double *zeroth, double *first, double *second,
    double *third, Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! (R_finite(theta[0]) && R_finite(theta[1])))
        error("all components of theta must be finite");
    if (! (R_finite(delta[0]) && R_finite(delta[1])))
        error("all components of delta must be finite");
    if (*deriv < 0 || *deriv > 3)
        error("deriv must be 0, 1, 2, or 3");
    if (! (delta[0] == 0.0 && delta[1] == 0.0))
        error("delta must be zero vector, no limits for normal");
    if (! (theta[1] < 0.0)) {
        *zeroth = R_PosInf;
        if (*deriv == 0) return;
        for (int i = 0; i < 2; i++)
            first[i] = R_NaN;
        if (*deriv <= 1) return;
        for (int i = 0; i < 4; i++)
            second[i] = R_NaN;
        if (*deriv <= 2) return;
        for (int i = 0; i < 8; i++)
            third[i] = R_NaN;
         return;
    }
    *zeroth = (- theta[0] * theta[0]) / (4.0 * theta[1]) +
        0.5 * log(- 1.0 / (2.0 * theta[1]));
    if (*deriv == 0) return;
    first[0] = (- theta[0]) / (2.0 * theta[1]);
    first[1] = ((theta[0] * theta[0]) / (4.0 * theta[1]) - 0.5) / theta[1];
    if (*deriv <= 1) return;
    second[0] = (- 0.5 / theta[1]);
    second[1] = second[2] = (0.5 * theta[0] / (theta[1] * theta[1]));
    second[3] = ((- theta[0] * theta[0]) / theta[1] + 1.0) /
        (2.0 * theta[1] * theta[1]);
    if (*deriv <= 2) return;
    third[0] = 0.0;
    third[1] = third[2] = third[4] = 0.5 / (theta[1] * theta[1]);
    third[3] = third[5] = third[6] = (- theta[0]) /
        (theta[1] * theta[1] * theta[1]);
    third[7] = (1.5 * (theta[0] * theta[0]) / theta[1] - 1.0) /
        (theta[1] * theta[1] * theta[1]);
}

#ifndef __GNUC__
static void normal_location_scale_link(double *xi, int *deriv,
    double *delta, double *zeroth, double *first, Family_t *that)
#else
static void normal_location_scale_link(double *xi, int *deriv,
    double *delta, double *zeroth, double *first,
    Family_t *that __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    if (! (R_finite(xi[0]) && R_finite(xi[1])))
        error("all components of xi must be finite");
    if (! (R_finite(delta[0]) && R_finite(delta[1])))
        error("all components of delta must be finite");
    if (*deriv < 0 || *deriv > 1)
        error("deriv must be 0 or 1");
    if (! (delta[0] == 0.0 && delta[1] == 0.0))
        error("delta must be zero vector, no limits for normal");
    double xi1 = xi[0];
    double xi2 = xi[1];
    double sigmasq = xi2 - xi1 * xi1;
    if (sigmasq <= 0.0)
        error("must have xi[2] - xi[1]^2 > 0 in normal location-scale");
    zeroth[0] = xi1 / sigmasq;
    zeroth[1] = (- 0.5 / sigmasq);
    if (*deriv == 1) {
        first[0] = (xi2 + xi1 * xi1) / sigmasq / sigmasq;
        first[1] = first[2] = (- xi1) / sigmasq / sigmasq;
        first[3] = 0.5 / sigmasq / sigmasq;
    }
}

static void normal_location_scale_validate_pred(double *ypred)
{
    double x = ypred[0];
    if (x != trunc(x))
        error("normal location-scale predecessor noninteger");
    if (x < 0.0)
        error("normal location-scale predecessor negative");
}

#ifndef __GNUC__
static void normal_location_scale_validate_delta(int d, double *delta)
#else
static void normal_location_scale_validate_delta(
    int d __attribute__ ((unused)), double *delta)
#endif /* __GNUC__ */
{
    if (delta[0] != 0.0 || delta[1] != 0.0)
       error("delta != 0 not allowed for normal location-scale");
}

#ifndef __GNUC__
static void normal_location_scale_validate_resp(double ypred, int d,
    double *delta, double *resp)
#else
static void normal_location_scale_validate_resp(double ypred,
    int d __attribute__ ((unused)), double *delta __attribute__ ((unused)),
    double *resp)
#endif /* __GNUC__ */
{
    double y1 = resp[0];
    double y2 = resp[1];
    if (ypred == 0.0 && (y1 != 0.0 || y2 != 0.0))
        error("successor nonzero when predecessor zero");
    // DON'T COMPARE DOUBLES FOR EQUALTY !!!!!!!!!!!!!!!!!!!!
    // (maybe even above is bogus)
    double y1sq = y1 * y1;
    double qux = (ypred * y2 - y1sq) / fmax(ypred * y2, y1sq);
    if (ypred == 1.0 && fabs(qux) > tolerance)
        error ("y2 != y1^2 for normal location-scale sample size one");
    if (ypred > 1.0 && qux < (- tolerance))
        error ("ypred * y2 < y1^2 for normal location-scale");
}

#ifndef __GNUC__
static void normal_location_scale_validate_theta(int d, double *delta,
    double *theta)
#else
static void normal_location_scale_validate_theta(
    int d __attribute__ ((unused)), double *delta __attribute__ ((unused)),
    double *theta)
#endif /* __GNUC__ */
{
    if (theta[1] < 0.0)
        return;
    error("theta[2] not negative for normal location-scale");
}

#ifndef __GNUC__
static void normal_location_scale_validate_xi(int d, double *delta, double *xi)
#else
static void normal_location_scale_validate_xi(int d __attribute__ ((unused)),
    double *delta __attribute__ ((unused)), double *xi)
#endif /* __GNUC__ */
{
    double xi1 = xi[0];
    double xi2 = xi[1];
    // DON'T COMPARE DOUBLES FOR EQUALTY !!!!!!!!!!!!!!!!!!!!
    // NOPE (even more exclamation points) here we have to because
    // xi makes no statistical sense and doesn't map to valid values
    // of other parameter vectors unless this is satisfied
    double xi1sq = xi1 * xi1;
    if (xi2 - xi1sq <= 0)
        error ("xi[2] <= xi[1]^2 for normal location-scale");
}

#ifndef __GNUC__
static void normal_location_scale_constancy(int d, double *delta, int *nvec,
    double *vectors, double *rhs)
#else
static void normal_location_scale_constancy(int d __attribute__ ((unused)),
    double *delta __attribute__ ((unused)), int *nvec,
    double *vectors __attribute__ ((unused)),
    double *rhs __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    *nvec = 0;
}

#ifndef __GNUC__
static void normal_location_scale_start_theta(int d, double *theta)
#else
static void normal_location_scale_start_theta(int d __attribute__ ((unused)),
    double *theta)
#endif /* __GNUC__ */
{
    theta[0] = 0.0;
    theta[1] = (- 1.0);
}

#ifndef __GNUC__
static Family_t normal_location_scale_constructor(double *hyper1,
    double *hyper2)
#else
static Family_t normal_location_scale_constructor(double *hyper1
    __attribute__ ((unused)), double *hyper2 __attribute__ ((unused)))
#endif /* __GNUC__ */
{
    Family_t result = {
        "normal.location.scale",
        R_NaN,
        R_NaN,
        2,
        normal_location_scale_cumulant,
        normal_location_scale_link,
        normal_location_scale_validate_pred,
        normal_location_scale_validate_delta,
        normal_location_scale_validate_resp,
        normal_location_scale_validate_theta,
        normal_location_scale_validate_xi,
        normal_location_scale_constancy,
        normal_location_scale_start_theta
    };
    return result;
}

// end of family-specific stuff

struct SuperFamily {
    char *name;
    Family_t (*constructor)(double *hyper1, double *hyper2);
};

typedef struct SuperFamily SuperFamily_t;

static Family_t famtab[MAX_NFAM];
static int nfam = 0;

static SuperFamily_t superfamtab[] =
{
    {"bernoulli", bernoulli_constructor},
    {"poisson", poisson_constructor},
    {"zero.truncated.poisson", zero_truncated_poisson_constructor},
    {"multinomial", multinomial_constructor},
    {"normal.location.scale", normal_location_scale_constructor},
    {NULL, NULL},
};

void astfam_cumulant(double *theta, int *fam, int *deriv, double *delta,
    double *zeroth, double *first, double *second, double *third)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    famtab[foo].cumulant(theta, deriv, delta, zeroth, first, second, third,
        &famtab[foo]);
}

void astfam_link(double *xi, int *fam, int *deriv, double *delta,
    double *zeroth, double *first)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    famtab[foo].link(xi, deriv, delta, zeroth, first, &famtab[foo]);
}

void astfam_dimension(int *fam, int *dimen)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    *dimen = famtab[foo].dimension;
}

void astfam_clear(void)
{
    nfam = 0;
}

void astfam_set(char **name, double *hyper1, double *hyper2)
{
    if (nfam >= MAX_NFAM)
        error("not enough room in family table,\nreinstall package"
            " with larger MAX_NFAM (defined in Makevars)");

    for (int i = 0; ; i++) {
        if (i >= MAX_NFAM || superfamtab[i].name == NULL)
            error("family \"%s\" not (yet) implemented", *name);
        if (strcmp(superfamtab[i].name, *name) == 0) {
            famtab[nfam++] = superfamtab[i].constructor(hyper1, hyper2);
            break;
        }
    }
}

void astfam_nfam(int *nfam_out)
{
    *nfam_out = nfam;
}

void astfam_validate_pred(int *fam, double *ypred)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    if (! R_finite(*ypred))
        error("predecessor must be finite");
    famtab[foo].validate_pred(ypred);
}

void astfam_validate_delta(int *fam, int *dimen, double *delta)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("dimension of delta does not match dimension of family");
    for (int i = 0; i < d; i++)
        if (! R_finite(delta[i]))
            error("delta must have all components finite");
    famtab[foo].validate_delta(d, delta);
}

void astfam_validate_resp(int *fam, int *dimen, double *delta, double *ypred,
    double *resp)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("dimension of response does not match dimension of family");
    for (int i = 0; i < d; i++)
        if (! R_finite(resp[i]))
            error("response must have all components finite");
    famtab[foo].validate_pred(ypred);
    famtab[foo].validate_delta(d, delta);
    famtab[foo].validate_resp(*ypred, d, delta, resp);
}

void astfam_set_tolerance(double *tol_ptr)
{
    if (*tol_ptr <= 0.0)
        error("trying to set tolerance nonpositive");
    tolerance = *tol_ptr;
}

void astfam_reset_tolerance(void)
{
    tolerance = default_tolerance;
}

void astfam_validate_theta(int *fam, int *dimen, double *delta, double *theta)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("dimension of theta does not match dimension of family");
    for (int i = 0; i < d; i++) {
        if (! R_finite(delta[i]))
            error("delta must have all components finite");
        if (! R_finite(theta[i]))
            error("theta must have all components finite");
    }
    famtab[foo].validate_delta(d, delta);
    famtab[foo].validate_theta(d, delta, theta);
}

void astfam_validate_xi(int *fam, int *dimen, double *delta, double *xi)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("dimension of theta does not match dimension of family");
    for (int i = 0; i < d; i++) {
        if (! R_finite(delta[i]))
            error("delta must have all components finite");
        if (! R_finite(xi[i]))
            error("xi must have all components finite");
    }
    famtab[foo].validate_delta(d, delta);
    famtab[foo].validate_xi(d, delta, xi);
}

void astfam_constancy(int *fam, int *dimen, double *delta, int *nvec,
    double *vectors, double *rhs)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("dimension of vectors does not match dimension of family");
    famtab[foo].validate_delta(d, delta);
    famtab[foo].constancy(d, delta, nvec, vectors, rhs);
}

void astfam_start_theta(int *fam, int *dimen, double *theta)
{
    if (*fam < 1 || *fam > nfam)
        error("fam out of range");
    // convert from one-origin indexing (for R) to zero-origin indexing (for C)
    register int foo = (*fam) - 1;
    register int d = famtab[foo].dimension;
    if (*dimen != d)
        error("astfam_start_theta: given dimension does not match"
            " dimension of family");
    famtab[foo].start_theta(d, theta);
}

