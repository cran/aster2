
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "aster.h"
#include "constancy.h"
#include "families.h"

// astfam_cumulant

static R_NativePrimitiveArgType astfam_cumulant_types[8] = {REALSXP, INTSXP,
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle astfam_cumulant_styles[8] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT};

// astfam_link

static R_NativePrimitiveArgType astfam_link_types[6] = {REALSXP, INTSXP,
    INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle astfam_link_styles[6] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_OUT, R_ARG_OUT};

// astfam_dimension

static R_NativePrimitiveArgType astfam_dimension_types[2] = {INTSXP, INTSXP};

static R_NativeArgStyle astfam_dimension_styles[2] = {R_ARG_IN, R_ARG_OUT};

// astfam_set

static R_NativePrimitiveArgType astfam_set_types[3] = {STRSXP, REALSXP,
    REALSXP};

static R_NativeArgStyle astfam_set_styles[3] = {R_ARG_IN, R_ARG_IN, R_ARG_IN};

// astfam_set_tolerance

static R_NativePrimitiveArgType astfam_set_tolerance_types[1] = {REALSXP};

static R_NativeArgStyle astfam_set_tolerance_styles[1] = {R_ARG_IN};

// astfam_nfam

static R_NativePrimitiveArgType astfam_nfam_types[1] = {INTSXP};

static R_NativeArgStyle astfam_nfam_styles[1] = {R_ARG_OUT};

// astfam_validate_pred

static R_NativePrimitiveArgType astfam_validate_pred_types[2] = {INTSXP,
    REALSXP};

static R_NativeArgStyle astfam_validate_pred_styles[2] = {R_ARG_IN, R_ARG_IN};

// astfam_validate_delta

static R_NativePrimitiveArgType astfam_validate_delta_types[3] = {INTSXP,
    INTSXP, REALSXP};

static R_NativeArgStyle astfam_validate_delta_styles[3] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN};

// astfam_validate_resp

static R_NativePrimitiveArgType astfam_validate_resp_types[5] = {INTSXP,
    INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle astfam_validate_resp_styles[5] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN};

// astfam_validate_theta

static R_NativePrimitiveArgType astfam_validate_theta_types[4] = {INTSXP,
    INTSXP, REALSXP, REALSXP};

static R_NativeArgStyle astfam_validate_theta_styles[4] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN};

// astfam_constancy

static R_NativePrimitiveArgType astfam_constancy_types[6] = {INTSXP,
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP};

static R_NativeArgStyle astfam_constancy_styles[6] = {R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT};

// aster_validate

static R_NativePrimitiveArgType aster_validate_types[7] = {INTSXP, REALSXP,
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativeArgStyle aster_validate_styles[7] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN};

// aster_theta_to_phi

static R_NativePrimitiveArgType aster_theta_to_phi_types[10] = {INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP};

static R_NativeArgStyle aster_theta_to_phi_styles[10] = {R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN,
    R_ARG_OUT, R_ARG_OUT};

// aster_theta_to_xi

static R_NativePrimitiveArgType aster_theta_to_xi_types[9] = {INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle aster_theta_to_xi_styles[9] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_OUT, R_ARG_OUT};

// aster_xi_to_mu

static R_NativePrimitiveArgType aster_xi_to_mu_types[8] = {INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle aster_xi_to_mu_styles[8] = {R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_OUT, R_ARG_OUT};

// aster_validate_theta

static R_NativePrimitiveArgType aster_validate_theta_types[5] = {INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativeArgStyle aster_validate_theta_styles[5] = {R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN};

// aster_starting_theta

static R_NativePrimitiveArgType aster_starting_theta_types[4] = {INTSXP,
    INTSXP, INTSXP, REALSXP};

static R_NativeArgStyle aster_starting_theta_styles[4] = {R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_OUT};

// astfam_start_theta

static R_NativePrimitiveArgType astfam_start_theta_types[3] = {INTSXP,
    INTSXP, REALSXP};

static R_NativeArgStyle astfam_start_theta_styles[3] = {R_ARG_IN,
    R_ARG_IN, R_ARG_OUT};

// register

static R_CMethodDef cMethods[] = {
    {"astfam_cumulant", (DL_FUNC) &astfam_cumulant,
        8, astfam_cumulant_types, astfam_cumulant_styles},
    {"astfam_link", (DL_FUNC) &astfam_link,
        6, astfam_link_types, astfam_link_styles},
    {"astfam_dimension", (DL_FUNC) &astfam_dimension,
        2, astfam_dimension_types, astfam_dimension_styles},
    {"astfam_nfam", (DL_FUNC) &astfam_nfam,
        1, astfam_nfam_types, astfam_nfam_styles},
    {"astfam_clear", (DL_FUNC) &astfam_clear,
        0, NULL, NULL},
    {"astfam_set", (DL_FUNC) &astfam_set,
        3, astfam_set_types, astfam_set_styles},
    {"astfam_set_tolerance", (DL_FUNC) &astfam_set_tolerance,
        1, astfam_set_tolerance_types, astfam_set_tolerance_styles},
    {"astfam_reset_tolerance", (DL_FUNC) &astfam_reset_tolerance,
        0, NULL, NULL},
    {"astfam_validate_pred", (DL_FUNC) &astfam_validate_pred,
        2, astfam_validate_pred_types, astfam_validate_pred_styles},
    {"astfam_validate_delta", (DL_FUNC) &astfam_validate_delta,
        3, astfam_validate_delta_types, astfam_validate_delta_styles},
    {"astfam_validate_resp", (DL_FUNC) &astfam_validate_resp,
        5, astfam_validate_resp_types, astfam_validate_resp_styles},
    {"astfam_validate_theta", (DL_FUNC) &astfam_validate_theta,
        4, astfam_validate_theta_types, astfam_validate_theta_styles},
    {"astfam_validate_xi", (DL_FUNC) &astfam_validate_xi,
        4, astfam_validate_theta_types, astfam_validate_theta_styles},
    {"astfam_constancy", (DL_FUNC) &astfam_constancy,
        6, astfam_constancy_types, astfam_constancy_styles},
    {"astfam_start_theta", (DL_FUNC) &astfam_start_theta,
        3, astfam_start_theta_types, astfam_start_theta_styles},
    {"aster_validate", (DL_FUNC) &aster_validate,
        7, aster_validate_types, aster_validate_styles},
    {"aster_validate_theta", (DL_FUNC) &aster_validate_theta,
        5, aster_validate_theta_types, aster_validate_theta_styles},
    {"aster_validate_xi", (DL_FUNC) &aster_validate_xi,
        5, aster_validate_theta_types, aster_validate_theta_styles},
    {"aster_theta_to_phi", (DL_FUNC) &aster_theta_to_phi,
        10, aster_theta_to_phi_types, aster_theta_to_phi_styles},
    {"aster_phi_to_theta", (DL_FUNC) &aster_phi_to_theta,
        10, aster_theta_to_phi_types, aster_theta_to_phi_styles},
    {"aster_theta_to_xi", (DL_FUNC) &aster_theta_to_xi,
        9, aster_theta_to_xi_types, aster_theta_to_xi_styles},
    {"aster_xi_to_theta", (DL_FUNC) &aster_xi_to_theta,
        9, aster_theta_to_xi_types, aster_theta_to_xi_styles},
    {"aster_xi_to_mu", (DL_FUNC) &aster_xi_to_mu,
        8, aster_xi_to_mu_types, aster_xi_to_mu_styles},
    {"aster_mu_to_xi", (DL_FUNC) &aster_mu_to_xi,
        8, aster_xi_to_mu_types, aster_xi_to_mu_styles},
    {"aster_starting_theta", (DL_FUNC) &aster_starting_theta,
        4, aster_starting_theta_types, aster_starting_theta_styles},
    {NULL, NULL, 0, NULL, NULL}
};
 
static R_CallMethodDef callMethods[]  = {
    {"aster_constancy", (DL_FUNC) &aster_constancy, 5},
    {NULL, NULL, 0}
};

void R_init_aster2(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

