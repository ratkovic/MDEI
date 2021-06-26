// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bayesLasso
List bayesLasso(arma::vec y, arma::mat X, double alpha, double tol);
RcppExport SEXP _NithinPackage_bayesLasso(SEXP ySEXP, SEXP XSEXP, SEXP alphaSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesLasso(y, X, alpha, tol));
    return rcpp_result_gen;
END_RCPP
}
// updateEtausqinv
arma::vec updateEtausqinv(arma::vec y, arma::mat X, double alpha, arma::vec Etausqinv, double tol);
RcppExport SEXP _NithinPackage_updateEtausqinv(SEXP ySEXP, SEXP XSEXP, SEXP alphaSEXP, SEXP EtausqinvSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Etausqinv(EtausqinvSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(updateEtausqinv(y, X, alpha, Etausqinv, tol));
    return rcpp_result_gen;
END_RCPP
}
// GCV
arma::vec GCV(arma::vec y, arma::mat X, arma::vec alphas, double tol);
RcppExport SEXP _NithinPackage_GCV(SEXP ySEXP, SEXP XSEXP, SEXP alphasSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GCV(y, X, alphas, tol));
    return rcpp_result_gen;
END_RCPP
}
// bs2
arma::mat bs2(arma::vec x, unsigned int deg);
RcppExport SEXP _NithinPackage_bs2(SEXP xSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(bs2(x, deg));
    return rcpp_result_gen;
END_RCPP
}
// dbs2
arma::mat dbs2(arma::vec x, unsigned int deg);
RcppExport SEXP _NithinPackage_dbs2(SEXP xSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(dbs2(x, deg));
    return rcpp_result_gen;
END_RCPP
}
// bsme
arma::mat bsme(arma::vec x);
RcppExport SEXP _NithinPackage_bsme(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bsme(x));
    return rcpp_result_gen;
END_RCPP
}
// dbsme
arma::mat dbsme(arma::vec x);
RcppExport SEXP _NithinPackage_dbsme(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dbsme(x));
    return rcpp_result_gen;
END_RCPP
}
// makebs
arma::mat makebs(arma::mat X);
RcppExport SEXP _NithinPackage_makebs(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(makebs(X));
    return rcpp_result_gen;
END_RCPP
}
// splineBases
arma::mat splineBases(arma::mat X, int covs);
RcppExport SEXP _NithinPackage_splineBases(SEXP XSEXP, SEXP covsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type covs(covsSEXP);
    rcpp_result_gen = Rcpp::wrap(splineBases(X, covs));
    return rcpp_result_gen;
END_RCPP
}
// subSamp
arma::vec subSamp(arma::vec v);
RcppExport SEXP _NithinPackage_subSamp(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(subSamp(v));
    return rcpp_result_gen;
END_RCPP
}
// correlations
List correlations(int obs, int covs, arma::mat X, arma::vec y, arma::vec treat, long long unsigned int a);
RcppExport SEXP _NithinPackage_correlations(SEXP obsSEXP, SEXP covsSEXP, SEXP XSEXP, SEXP ySEXP, SEXP treatSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< long long unsigned int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(correlations(obs, covs, X, y, treat, a));
    return rcpp_result_gen;
END_RCPP
}
// checkcor
arma::vec checkcor(arma::mat cors, double thresh);
RcppExport SEXP _NithinPackage_checkcor(SEXP corsSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cors(corsSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(checkcor(cors, thresh));
    return rcpp_result_gen;
END_RCPP
}
// main
int main();
RcppExport SEXP _NithinPackage_main() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(main());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NithinPackage_bayesLasso", (DL_FUNC) &_NithinPackage_bayesLasso, 4},
    {"_NithinPackage_updateEtausqinv", (DL_FUNC) &_NithinPackage_updateEtausqinv, 5},
    {"_NithinPackage_GCV", (DL_FUNC) &_NithinPackage_GCV, 4},
    {"_NithinPackage_bs2", (DL_FUNC) &_NithinPackage_bs2, 2},
    {"_NithinPackage_dbs2", (DL_FUNC) &_NithinPackage_dbs2, 2},
    {"_NithinPackage_bsme", (DL_FUNC) &_NithinPackage_bsme, 1},
    {"_NithinPackage_dbsme", (DL_FUNC) &_NithinPackage_dbsme, 1},
    {"_NithinPackage_makebs", (DL_FUNC) &_NithinPackage_makebs, 1},
    {"_NithinPackage_splineBases", (DL_FUNC) &_NithinPackage_splineBases, 2},
    {"_NithinPackage_subSamp", (DL_FUNC) &_NithinPackage_subSamp, 1},
    {"_NithinPackage_correlations", (DL_FUNC) &_NithinPackage_correlations, 6},
    {"_NithinPackage_checkcor", (DL_FUNC) &_NithinPackage_checkcor, 2},
    {"_NithinPackage_main", (DL_FUNC) &_NithinPackage_main, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_NithinPackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
