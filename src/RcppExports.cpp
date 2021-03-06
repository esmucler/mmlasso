// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// rho_marina
arma::vec rho_marina(arma::vec x, double cw);
RcppExport SEXP mmlasso_rho_marina(SEXP xSEXP, SEXP cwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type cw(cwSEXP );
        arma::vec __result = rho_marina(x, cw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rho_bisq
arma::vec rho_bisq(arma::vec x, double cw);
RcppExport SEXP mmlasso_rho_bisq(SEXP xSEXP, SEXP cwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type cw(cwSEXP );
        arma::vec __result = rho_bisq(x, cw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// psi_marina
arma::vec psi_marina(arma::vec x, double cw);
RcppExport SEXP mmlasso_psi_marina(SEXP xSEXP, SEXP cwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type cw(cwSEXP );
        arma::vec __result = psi_marina(x, cw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// weight_marina
arma::vec weight_marina(arma::vec x, double cw);
RcppExport SEXP mmlasso_weight_marina(SEXP xSEXP, SEXP cwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type cw(cwSEXP );
        arma::vec __result = weight_marina(x, cw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Mscale_mar
double Mscale_mar(arma::vec x, double b, double cc);
RcppExport SEXP mmlasso_Mscale_mar(SEXP xSEXP, SEXP bSEXP, SEXP ccSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type b(bSEXP );
        Rcpp::traits::input_parameter< double >::type cc(ccSEXP );
        double __result = Mscale_mar(x, b, cc);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// scale_tau
double scale_tau(arma::vec x);
RcppExport SEXP mmlasso_scale_tau(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        double __result = scale_tau(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// my_svdecon
List my_svdecon(arma::mat x);
RcppExport SEXP mmlasso_my_svdecon(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP );
        List __result = my_svdecon(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// SPCC
List SPCC(arma::mat x);
RcppExport SEXP mmlasso_SPCC(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP );
        List __result = SPCC(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// MMLassoCpp_ini
List MMLassoCpp_ini(arma::mat xx, arma::vec y, arma::vec beta_ini);
RcppExport SEXP mmlasso_MMLassoCpp_ini(SEXP xxSEXP, SEXP ySEXP, SEXP beta_iniSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type xx(xxSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::vec >::type beta_ini(beta_iniSEXP );
        List __result = MMLassoCpp_ini(xx, y, beta_ini);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// MMLassoCpp1
List MMLassoCpp1(arma::mat x, arma::vec y, arma::vec beta_ini, double scale_ini, double c1);
RcppExport SEXP mmlasso_MMLassoCpp1(SEXP xSEXP, SEXP ySEXP, SEXP beta_iniSEXP, SEXP scale_iniSEXP, SEXP c1SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::vec >::type beta_ini(beta_iniSEXP );
        Rcpp::traits::input_parameter< double >::type scale_ini(scale_iniSEXP );
        Rcpp::traits::input_parameter< double >::type c1(c1SEXP );
        List __result = MMLassoCpp1(x, y, beta_ini, scale_ini, c1);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// MMLassoCpp2
List MMLassoCpp2(arma::vec xjota, arma::vec yast, arma::vec beta_lars, arma::vec beta_o, arma::vec alpha);
RcppExport SEXP mmlasso_MMLassoCpp2(SEXP xjotaSEXP, SEXP yastSEXP, SEXP beta_larsSEXP, SEXP beta_oSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type xjota(xjotaSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type yast(yastSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type beta_lars(beta_larsSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type beta_o(beta_oSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP );
        List __result = MMLassoCpp2(xjota, yast, beta_lars, beta_o, alpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rr_se
List rr_se(arma::mat X, arma::vec y, double lambda2, double deltaesc, double cc_scale, int nkeep, int niter, double epsilon);
RcppExport SEXP mmlasso_rr_se(SEXP XSEXP, SEXP ySEXP, SEXP lambda2SEXP, SEXP deltaescSEXP, SEXP cc_scaleSEXP, SEXP nkeepSEXP, SEXP niterSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP );
        Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP );
        Rcpp::traits::input_parameter< double >::type deltaesc(deltaescSEXP );
        Rcpp::traits::input_parameter< double >::type cc_scale(cc_scaleSEXP );
        Rcpp::traits::input_parameter< int >::type nkeep(nkeepSEXP );
        Rcpp::traits::input_parameter< int >::type niter(niterSEXP );
        Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP );
        List __result = rr_se(X, y, lambda2, deltaesc, cc_scale, nkeep, niter, epsilon);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
