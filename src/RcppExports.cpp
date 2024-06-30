// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CPPlwls1d
Eigen::VectorXd CPPlwls1d(const double& bw, const std::string kernel_type, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xin, const Eigen::Map<Eigen::VectorXd>& yin, const Eigen::Map<Eigen::VectorXd>& xout, const unsigned int& npoly, const unsigned int& nder);
RcppExport SEXP _fdapace_CPPlwls1d(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP winSEXP, SEXP xinSEXP, SEXP yinSEXP, SEXP xoutSEXP, SEXP npolySEXP, SEXP nderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xin(xinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type yin(yinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type npoly(npolySEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nder(nderSEXP);
    rcpp_result_gen = Rcpp::wrap(CPPlwls1d(bw, kernel_type, win, xin, yin, xout, npoly, nder));
    return rcpp_result_gen;
END_RCPP
}
// GetIndCEScoresCPP
Rcpp::List GetIndCEScoresCPP(const Eigen::Map<Eigen::VectorXd>& yVec, const Eigen::Map<Eigen::VectorXd>& muVec, const Eigen::Map<Eigen::VectorXd>& lamVec, const Eigen::Map<Eigen::MatrixXd>& phiMat, const Eigen::Map<Eigen::MatrixXd>& SigmaYi);
RcppExport SEXP _fdapace_GetIndCEScoresCPP(SEXP yVecSEXP, SEXP muVecSEXP, SEXP lamVecSEXP, SEXP phiMatSEXP, SEXP SigmaYiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type muVec(muVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type lamVec(lamVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type phiMat(phiMatSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type SigmaYi(SigmaYiSEXP);
    rcpp_result_gen = Rcpp::wrap(GetIndCEScoresCPP(yVec, muVec, lamVec, phiMat, SigmaYi));
    return rcpp_result_gen;
END_RCPP
}
// GetIndCEScoresCPPnewInd
Rcpp::List GetIndCEScoresCPPnewInd(const Eigen::Map<Eigen::VectorXd>& yVec, const Eigen::Map<Eigen::VectorXd>& muVec, const Eigen::Map<Eigen::VectorXd>& lamVec, const Eigen::Map<Eigen::MatrixXd>& phiMat, const Eigen::Map<Eigen::MatrixXd>& SigmaYi, const Eigen::Map<Eigen::MatrixXd>& newPhi, const Eigen::Map<Eigen::VectorXd>& newMu);
RcppExport SEXP _fdapace_GetIndCEScoresCPPnewInd(SEXP yVecSEXP, SEXP muVecSEXP, SEXP lamVecSEXP, SEXP phiMatSEXP, SEXP SigmaYiSEXP, SEXP newPhiSEXP, SEXP newMuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type muVec(muVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type lamVec(lamVecSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type phiMat(phiMatSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type SigmaYi(SigmaYiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type newPhi(newPhiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type newMu(newMuSEXP);
    rcpp_result_gen = Rcpp::wrap(GetIndCEScoresCPPnewInd(yVec, muVec, lamVec, phiMat, SigmaYi, newPhi, newMu));
    return rcpp_result_gen;
END_RCPP
}
// RCPPmean
double RCPPmean(const Rcpp::NumericVector X);
RcppExport SEXP _fdapace_RCPPmean(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPmean(X));
    return rcpp_result_gen;
END_RCPP
}
// RCPPvar
double RCPPvar(const Rcpp::NumericVector X);
RcppExport SEXP _fdapace_RCPPvar(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPvar(X));
    return rcpp_result_gen;
END_RCPP
}
// RcppPseudoApprox
Eigen::VectorXd RcppPseudoApprox(const Eigen::Map<Eigen::VectorXd>& X, const Eigen::Map<Eigen::VectorXd>& Y, const Eigen::Map<Eigen::VectorXd>& X_target);
RcppExport SEXP _fdapace_RcppPseudoApprox(SEXP XSEXP, SEXP YSEXP, SEXP X_targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type X_target(X_targetSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppPseudoApprox(X, Y, X_target));
    return rcpp_result_gen;
END_RCPP
}
// Rcppsort
NumericVector Rcppsort(NumericVector v);
RcppExport SEXP _fdapace_Rcppsort(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcppsort(v));
    return rcpp_result_gen;
END_RCPP
}
// Rmullwlsk
Eigen::MatrixXd Rmullwlsk(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xgrid, const Eigen::Map<Eigen::VectorXd>& ygrid, const bool& bwCheck, const bool& transp);
RcppExport SEXP _fdapace_Rmullwlsk(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xgridSEXP, SEXP ygridSEXP, SEXP bwCheckSEXP, SEXP transpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    Rcpp::traits::input_parameter< const bool& >::type transp(transpSEXP);
    rcpp_result_gen = Rcpp::wrap(Rmullwlsk(bw, kernel_type, tPairs, cxxn, win, xgrid, ygrid, bwCheck, transp));
    return rcpp_result_gen;
END_RCPP
}
// RmullwlskCC
Eigen::MatrixXd RmullwlskCC(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xgrid, const Eigen::Map<Eigen::VectorXd>& ygrid, const bool& bwCheck);
RcppExport SEXP _fdapace_RmullwlskCC(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xgridSEXP, SEXP ygridSEXP, SEXP bwCheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    rcpp_result_gen = Rcpp::wrap(RmullwlskCC(bw, kernel_type, tPairs, cxxn, win, xgrid, ygrid, bwCheck));
    return rcpp_result_gen;
END_RCPP
}
// RmullwlskCCsort2
Eigen::MatrixXd RmullwlskCCsort2(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xgrid, const Eigen::Map<Eigen::VectorXd>& ygrid, const bool& bwCheck);
RcppExport SEXP _fdapace_RmullwlskCCsort2(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xgridSEXP, SEXP ygridSEXP, SEXP bwCheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    rcpp_result_gen = Rcpp::wrap(RmullwlskCCsort2(bw, kernel_type, tPairs, cxxn, win, xgrid, ygrid, bwCheck));
    return rcpp_result_gen;
END_RCPP
}
// RmullwlskUniversal
Eigen::MatrixXd RmullwlskUniversal(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xgrid, const Eigen::Map<Eigen::VectorXd>& ygrid, const bool& bwCheck, const bool& autoCov);
RcppExport SEXP _fdapace_RmullwlskUniversal(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xgridSEXP, SEXP ygridSEXP, SEXP bwCheckSEXP, SEXP autoCovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    Rcpp::traits::input_parameter< const bool& >::type autoCov(autoCovSEXP);
    rcpp_result_gen = Rcpp::wrap(RmullwlskUniversal(bw, kernel_type, tPairs, cxxn, win, xgrid, ygrid, bwCheck, autoCov));
    return rcpp_result_gen;
END_RCPP
}
// RmullwlskUniversalDeriv
Eigen::MatrixXd RmullwlskUniversalDeriv(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xgrid, const Eigen::Map<Eigen::VectorXd>& ygrid, const int& npoly, const int& nder1, const int& nder2, const bool& bwCheck, const bool& autoCov);
RcppExport SEXP _fdapace_RmullwlskUniversalDeriv(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xgridSEXP, SEXP ygridSEXP, SEXP npolySEXP, SEXP nder1SEXP, SEXP nder2SEXP, SEXP bwCheckSEXP, SEXP autoCovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< const int& >::type npoly(npolySEXP);
    Rcpp::traits::input_parameter< const int& >::type nder1(nder1SEXP);
    Rcpp::traits::input_parameter< const int& >::type nder2(nder2SEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    Rcpp::traits::input_parameter< const bool& >::type autoCov(autoCovSEXP);
    rcpp_result_gen = Rcpp::wrap(RmullwlskUniversalDeriv(bw, kernel_type, tPairs, cxxn, win, xgrid, ygrid, npoly, nder1, nder2, bwCheck, autoCov));
    return rcpp_result_gen;
END_RCPP
}
// Rrotatedmullwlsk
Eigen::VectorXd Rrotatedmullwlsk(const Eigen::Map<Eigen::VectorXd>& bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd>& tPairs, const Eigen::Map<Eigen::MatrixXd>& cxxn, const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::MatrixXd>& xygrid, const unsigned int npoly, const bool& bwCheck);
RcppExport SEXP _fdapace_Rrotatedmullwlsk(SEXP bwSEXP, SEXP kernel_typeSEXP, SEXP tPairsSEXP, SEXP cxxnSEXP, SEXP winSEXP, SEXP xygridSEXP, SEXP npolySEXP, SEXP bwCheckSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const std::string >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type tPairs(tPairsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type cxxn(cxxnSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type xygrid(xygridSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type npoly(npolySEXP);
    Rcpp::traits::input_parameter< const bool& >::type bwCheck(bwCheckSEXP);
    rcpp_result_gen = Rcpp::wrap(Rrotatedmullwlsk(bw, kernel_type, tPairs, cxxn, win, xygrid, npoly, bwCheck));
    return rcpp_result_gen;
END_RCPP
}
// cumtrapzRcpp
Rcpp::NumericVector cumtrapzRcpp(const Rcpp::NumericVector X, const Rcpp::NumericVector Y);
RcppExport SEXP _fdapace_cumtrapzRcpp(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(cumtrapzRcpp(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// dropZeroElementsXYWin
Eigen::MatrixXd dropZeroElementsXYWin(const Eigen::Map<Eigen::VectorXd>& win, const Eigen::Map<Eigen::VectorXd>& xin, const Eigen::Map<Eigen::VectorXd>& yin);
RcppExport SEXP _fdapace_dropZeroElementsXYWin(SEXP winSEXP, SEXP xinSEXP, SEXP yinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type win(winSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xin(xinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type yin(yinSEXP);
    rcpp_result_gen = Rcpp::wrap(dropZeroElementsXYWin(win, xin, yin));
    return rcpp_result_gen;
END_RCPP
}
// interp2lin
Eigen::VectorXd interp2lin(const Eigen::Map<Eigen::VectorXd>& xin, const Eigen::Map<Eigen::VectorXd>& yin, const Eigen::Map<Eigen::VectorXd>& zin, const Eigen::Map<Eigen::VectorXd>& xou, const Eigen::Map<Eigen::VectorXd>& you);
RcppExport SEXP _fdapace_interp2lin(SEXP xinSEXP, SEXP yinSEXP, SEXP zinSEXP, SEXP xouSEXP, SEXP youSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xin(xinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type yin(yinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type zin(zinSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type xou(xouSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type you(youSEXP);
    rcpp_result_gen = Rcpp::wrap(interp2lin(xin, yin, zin, xou, you));
    return rcpp_result_gen;
END_RCPP
}
// trapzRcpp
double trapzRcpp(const Rcpp::NumericVector X, const Rcpp::NumericVector Y);
RcppExport SEXP _fdapace_trapzRcpp(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(trapzRcpp(X, Y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fdapace_CPPlwls1d", (DL_FUNC) &_fdapace_CPPlwls1d, 8},
    {"_fdapace_GetIndCEScoresCPP", (DL_FUNC) &_fdapace_GetIndCEScoresCPP, 5},
    {"_fdapace_GetIndCEScoresCPPnewInd", (DL_FUNC) &_fdapace_GetIndCEScoresCPPnewInd, 7},
    {"_fdapace_RCPPmean", (DL_FUNC) &_fdapace_RCPPmean, 1},
    {"_fdapace_RCPPvar", (DL_FUNC) &_fdapace_RCPPvar, 1},
    {"_fdapace_RcppPseudoApprox", (DL_FUNC) &_fdapace_RcppPseudoApprox, 3},
    {"_fdapace_Rcppsort", (DL_FUNC) &_fdapace_Rcppsort, 1},
    {"_fdapace_Rmullwlsk", (DL_FUNC) &_fdapace_Rmullwlsk, 9},
    {"_fdapace_RmullwlskCC", (DL_FUNC) &_fdapace_RmullwlskCC, 8},
    {"_fdapace_RmullwlskCCsort2", (DL_FUNC) &_fdapace_RmullwlskCCsort2, 8},
    {"_fdapace_RmullwlskUniversal", (DL_FUNC) &_fdapace_RmullwlskUniversal, 9},
    {"_fdapace_RmullwlskUniversalDeriv", (DL_FUNC) &_fdapace_RmullwlskUniversalDeriv, 12},
    {"_fdapace_Rrotatedmullwlsk", (DL_FUNC) &_fdapace_Rrotatedmullwlsk, 8},
    {"_fdapace_cumtrapzRcpp", (DL_FUNC) &_fdapace_cumtrapzRcpp, 2},
    {"_fdapace_dropZeroElementsXYWin", (DL_FUNC) &_fdapace_dropZeroElementsXYWin, 3},
    {"_fdapace_interp2lin", (DL_FUNC) &_fdapace_interp2lin, 5},
    {"_fdapace_trapzRcpp", (DL_FUNC) &_fdapace_trapzRcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_fdapace(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
