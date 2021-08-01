// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// setupGCV
List setupGCV(arma::vec y, arma::mat X, arma::vec alphas);
RcppExport SEXP _NithinPackage_setupGCV(SEXP ySEXP, SEXP XSEXP, SEXP alphasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alphas(alphasSEXP);
    rcpp_result_gen = Rcpp::wrap(setupGCV(y, X, alphas));
    return rcpp_result_gen;
END_RCPP
}
// update
List update(List L, double tol);
RcppExport SEXP _NithinPackage_update(SEXP LSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update(L, tol));
    return rcpp_result_gen;
END_RCPP
}
// GCV
List GCV(arma::vec y, arma::mat X, arma::vec alphas, double tol);
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
// myrank
arma::vec myrank(arma::vec x);
RcppExport SEXP _NithinPackage_myrank(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(myrank(x));
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
// namesAndCorrs
List namesAndCorrs(arma::mat XSubsamp, std::vector<std::string> Xnames, arma::vec ySubsamp, std::vector<int> colSizes, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct, std::vector<std::string> treatNames, long long unsigned int a);
RcppExport SEXP _NithinPackage_namesAndCorrs(SEXP XSubsampSEXP, SEXP XnamesSEXP, SEXP ySubsampSEXP, SEXP colSizesSEXP, SEXP treatSubsampSEXP, SEXP XConstructSEXP, SEXP treatConstructSEXP, SEXP treatNamesSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XSubsamp(XSubsampSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type Xnames(XnamesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ySubsamp(ySubsampSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type colSizes(colSizesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type treatSubsamp(treatSubsampSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XConstruct(XConstructSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type treatConstruct(treatConstructSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type treatNames(treatNamesSEXP);
    Rcpp::traits::input_parameter< long long unsigned int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(namesAndCorrs(XSubsamp, Xnames, ySubsamp, colSizes, treatSubsamp, XConstruct, treatConstruct, treatNames, a));
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
// splineCorrs
List splineCorrs(arma::mat XSubsamp, arma::vec ySubsamp, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct, long long unsigned int a);
RcppExport SEXP _NithinPackage_splineCorrs(SEXP XSubsampSEXP, SEXP ySubsampSEXP, SEXP treatSubsampSEXP, SEXP XConstructSEXP, SEXP treatConstructSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XSubsamp(XSubsampSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ySubsamp(ySubsampSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type treatSubsamp(treatSubsampSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XConstruct(XConstructSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type treatConstruct(treatConstructSEXP);
    Rcpp::traits::input_parameter< long long unsigned int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(splineCorrs(XSubsamp, ySubsamp, treatSubsamp, XConstruct, treatConstruct, a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NithinPackage_bayesLasso", (DL_FUNC) &_NithinPackage_bayesLasso, 4},
    {"_NithinPackage_setupGCV", (DL_FUNC) &_NithinPackage_setupGCV, 3},
    {"_NithinPackage_update", (DL_FUNC) &_NithinPackage_update, 2},
    {"_NithinPackage_GCV", (DL_FUNC) &_NithinPackage_GCV, 4},
    {"_NithinPackage_myrank", (DL_FUNC) &_NithinPackage_myrank, 1},
    {"_NithinPackage_checkcor", (DL_FUNC) &_NithinPackage_checkcor, 2},
    {"_NithinPackage_subSamp", (DL_FUNC) &_NithinPackage_subSamp, 1},
    {"_NithinPackage_namesAndCorrs", (DL_FUNC) &_NithinPackage_namesAndCorrs, 9},
    {"_NithinPackage_checkcor", (DL_FUNC) &_NithinPackage_checkcor, 2},
    {"_NithinPackage_splineCorrs", (DL_FUNC) &_NithinPackage_splineCorrs, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_NithinPackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
