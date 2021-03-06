/*
 *  AKAR4_gvf.c
 *
 *  GSL C file for the vector field named: AKAR4
 *
 *  This file was generated by the program VFGEN, version: 2.6.0.dev1
 *  Generated on 25-Feb-2022 at 15:40
 */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/*
 *  The vector field.
 */

int AKAR4_vf(double t, const double y_[], double f_[], void *params)
{
    double AKAR4, AKAR4_C, AKAR4p, C;
    double kf_C_AKAR4, kb_C_AKAR4, kcat_AKARp;
    double reaction_1, reaction_2;
    double *p_;

    p_ = (double *) params;

    AKAR4      = y_[0];
    AKAR4_C    = y_[1];
    AKAR4p     = y_[2];
    C          = y_[3];

    kf_C_AKAR4 = p_[0];
    kb_C_AKAR4 = p_[1];
    kcat_AKARp = p_[2];

    reaction_1 = -kb_C_AKAR4*AKAR4_C+AKAR4*C*kf_C_AKAR4;
    reaction_2 = kcat_AKARp*AKAR4_C;

    f_[0] = -reaction_1;
    f_[1] =  reaction_1-reaction_2;
    f_[2] = reaction_2;
    f_[3] = -reaction_1+reaction_2;

    return GSL_SUCCESS;
}

/*
 *  The Jacobian.
 */

int AKAR4_jac(double t, const double y_[], double *jac_, double *dfdt_, void *params)
{
    double AKAR4, AKAR4_C, AKAR4p, C;
    double kf_C_AKAR4, kb_C_AKAR4, kcat_AKARp;
    double *p_;

    p_ = (double *) params;

    AKAR4      = y_[0];
    AKAR4_C    = y_[1];
    AKAR4p     = y_[2];
    C          = y_[3];

    kf_C_AKAR4 = p_[0];
    kb_C_AKAR4 = p_[1];
    kcat_AKARp = p_[2];

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(jac_,4,4);
    gsl_matrix *m_ = &dfdy_mat.matrix;

    gsl_matrix_set(m_, 0, 0, -C*kf_C_AKAR4);
    gsl_matrix_set(m_, 0, 1, kb_C_AKAR4);
    gsl_matrix_set(m_, 0, 2, 0.0);
    gsl_matrix_set(m_, 0, 3, -AKAR4*kf_C_AKAR4);
    gsl_matrix_set(m_, 1, 0, C*kf_C_AKAR4);
    gsl_matrix_set(m_, 1, 1, -kb_C_AKAR4-kcat_AKARp);
    gsl_matrix_set(m_, 1, 2, 0.0);
    gsl_matrix_set(m_, 1, 3, AKAR4*kf_C_AKAR4);
    gsl_matrix_set(m_, 2, 0, 0.0);
    gsl_matrix_set(m_, 2, 1, kcat_AKARp);
    gsl_matrix_set(m_, 2, 2, 0.0);
    gsl_matrix_set(m_, 2, 3, 0.0);
    gsl_matrix_set(m_, 3, 0, -C*kf_C_AKAR4);
    gsl_matrix_set(m_, 3, 1,  kb_C_AKAR4+kcat_AKARp);
    gsl_matrix_set(m_, 3, 2, 0.0);
    gsl_matrix_set(m_, 3, 3, -AKAR4*kf_C_AKAR4);

    dfdt_[0] = 0.0;
    dfdt_[1] = 0.0;
    dfdt_[2] = 0.0;
    dfdt_[3] = 0.0;

    return GSL_SUCCESS;
}

/*
 *  The Jacobian with respect to the parameters.
 */

int AKAR4_jacp(double t, const double y_[], double *jacp_, void *params)
{
    double AKAR4, AKAR4_C, AKAR4p, C;
    double kf_C_AKAR4, kb_C_AKAR4, kcat_AKARp;
    double *p_;

    p_ = (double *) params;

    AKAR4      = y_[0];
    AKAR4_C    = y_[1];
    AKAR4p     = y_[2];
    C          = y_[3];

    kf_C_AKAR4 = p_[0];
    kb_C_AKAR4 = p_[1];
    kcat_AKARp = p_[2];

    gsl_matrix_view dfdp_mat = gsl_matrix_view_array(jacp_,4,3);
    gsl_matrix *m_ = &dfdp_mat.matrix;

    gsl_matrix_set(m_, 0, 0, -AKAR4*C);
    gsl_matrix_set(m_, 0, 1, AKAR4_C);
    gsl_matrix_set(m_, 0, 2, 0.0);
    gsl_matrix_set(m_, 1, 0, AKAR4*C);
    gsl_matrix_set(m_, 1, 1, -AKAR4_C);
    gsl_matrix_set(m_, 1, 2, -AKAR4_C);
    gsl_matrix_set(m_, 2, 0, 0.0);
    gsl_matrix_set(m_, 2, 1, 0.0);
    gsl_matrix_set(m_, 2, 2, AKAR4_C);
    gsl_matrix_set(m_, 3, 0, -AKAR4*C);
    gsl_matrix_set(m_, 3, 1, AKAR4_C);
    gsl_matrix_set(m_, 3, 2, AKAR4_C);

    return GSL_SUCCESS;
}
