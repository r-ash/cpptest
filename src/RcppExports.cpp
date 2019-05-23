// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// add
double add(double a, double b);
RcppExport SEXP _cpptest_add(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(add(a, b));
    return rcpp_result_gen;
END_RCPP
}
// do_double
double do_double(double x);
RcppExport SEXP _cpptest_do_double(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(do_double(x));
    return rcpp_result_gen;
END_RCPP
}
// modify_std_vector
std::vector<double> modify_std_vector(std::vector<double> vector);
RcppExport SEXP _cpptest_modify_std_vector(SEXP vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type vector(vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(modify_std_vector(vector));
    return rcpp_result_gen;
END_RCPP
}
// get_array
void get_array();
RcppExport SEXP _cpptest_get_array() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    get_array();
    return R_NilValue;
END_RCPP
}
// push_array
std::vector<double> push_array(std::vector<double> arr);
RcppExport SEXP _cpptest_push_array(SEXP arrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type arr(arrSEXP);
    rcpp_result_gen = Rcpp::wrap(push_array(arr));
    return rcpp_result_gen;
END_RCPP
}
// push_list_arrays
std::list< std::vector<double> > push_list_arrays(std::list< std::vector<double> > lst);
RcppExport SEXP _cpptest_push_list_arrays(SEXP lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::list< std::vector<double> > >::type lst(lstSEXP);
    rcpp_result_gen = Rcpp::wrap(push_list_arrays(lst));
    return rcpp_result_gen;
END_RCPP
}
// runModel
std::vector<double> runModel(std::vector<double> basePop, std::vector<int> ageGroupsSpan, int timeArtStart, SEXP entrantPrev, std::vector<double> vertTransLag, std::vector<double> paedSurveyLag, bool populationAdjust, std::vector<double> entrantPop, std::vector<double> birthLag, std::vector<double> cumSurv, std::vector<double> cumNetMigr, double netMigrHivProb, std::vector<double> paedSurvCd4Distrib, SEXP entrantArtCoverage, std::vector<double> paedSurvArtCd4Distrib, std::vector<double> survRate, std::vector<double> netMigr, std::vector<double> asfRate, std::vector<double> sexRatioBirth, int hivStepsPerYear, std::vector<double> cd4Prog, std::vector<double> cd4InitDist, std::vector<double> cd4Mort, std::vector<double> incrrAges, double relinfectArt, SEXP iota, std::vector<double> incrrSex, SEXP incidMod, int eppMod, int scaleCd4Mort, std::vector<double> projSteps, SEXP tsEpidemicStart, std::vector<int> artCd4EligId, std::vector<double> specPopPercentElig, std::vector<double> pregnantWomenArtElig, double who34PercentElig, SEXP rSplineRVec, SEXP rTrendBeta, SEXP rTrendTStab, SEXP rTrendR0, int timeSteps);
RcppExport SEXP _cpptest_runModel(SEXP basePopSEXP, SEXP ageGroupsSpanSEXP, SEXP timeArtStartSEXP, SEXP entrantPrevSEXP, SEXP vertTransLagSEXP, SEXP paedSurveyLagSEXP, SEXP populationAdjustSEXP, SEXP entrantPopSEXP, SEXP birthLagSEXP, SEXP cumSurvSEXP, SEXP cumNetMigrSEXP, SEXP netMigrHivProbSEXP, SEXP paedSurvCd4DistribSEXP, SEXP entrantArtCoverageSEXP, SEXP paedSurvArtCd4DistribSEXP, SEXP survRateSEXP, SEXP netMigrSEXP, SEXP asfRateSEXP, SEXP sexRatioBirthSEXP, SEXP hivStepsPerYearSEXP, SEXP cd4ProgSEXP, SEXP cd4InitDistSEXP, SEXP cd4MortSEXP, SEXP incrrAgesSEXP, SEXP relinfectArtSEXP, SEXP iotaSEXP, SEXP incrrSexSEXP, SEXP incidModSEXP, SEXP eppModSEXP, SEXP scaleCd4MortSEXP, SEXP projStepsSEXP, SEXP tsEpidemicStartSEXP, SEXP artCd4EligIdSEXP, SEXP specPopPercentEligSEXP, SEXP pregnantWomenArtEligSEXP, SEXP who34PercentEligSEXP, SEXP rSplineRVecSEXP, SEXP rTrendBetaSEXP, SEXP rTrendTStabSEXP, SEXP rTrendR0SEXP, SEXP timeStepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type basePop(basePopSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ageGroupsSpan(ageGroupsSpanSEXP);
    Rcpp::traits::input_parameter< int >::type timeArtStart(timeArtStartSEXP);
    Rcpp::traits::input_parameter< SEXP >::type entrantPrev(entrantPrevSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type vertTransLag(vertTransLagSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type paedSurveyLag(paedSurveyLagSEXP);
    Rcpp::traits::input_parameter< bool >::type populationAdjust(populationAdjustSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type entrantPop(entrantPopSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type birthLag(birthLagSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type cumSurv(cumSurvSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type cumNetMigr(cumNetMigrSEXP);
    Rcpp::traits::input_parameter< double >::type netMigrHivProb(netMigrHivProbSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type paedSurvCd4Distrib(paedSurvCd4DistribSEXP);
    Rcpp::traits::input_parameter< SEXP >::type entrantArtCoverage(entrantArtCoverageSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type paedSurvArtCd4Distrib(paedSurvArtCd4DistribSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type survRate(survRateSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type netMigr(netMigrSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type asfRate(asfRateSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type sexRatioBirth(sexRatioBirthSEXP);
    Rcpp::traits::input_parameter< int >::type hivStepsPerYear(hivStepsPerYearSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type cd4Prog(cd4ProgSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type cd4InitDist(cd4InitDistSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type cd4Mort(cd4MortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type incrrAges(incrrAgesSEXP);
    Rcpp::traits::input_parameter< double >::type relinfectArt(relinfectArtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type iota(iotaSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type incrrSex(incrrSexSEXP);
    Rcpp::traits::input_parameter< SEXP >::type incidMod(incidModSEXP);
    Rcpp::traits::input_parameter< int >::type eppMod(eppModSEXP);
    Rcpp::traits::input_parameter< int >::type scaleCd4Mort(scaleCd4MortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type projSteps(projStepsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tsEpidemicStart(tsEpidemicStartSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type artCd4EligId(artCd4EligIdSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type specPopPercentElig(specPopPercentEligSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pregnantWomenArtElig(pregnantWomenArtEligSEXP);
    Rcpp::traits::input_parameter< double >::type who34PercentElig(who34PercentEligSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rSplineRVec(rSplineRVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rTrendBeta(rTrendBetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rTrendTStab(rTrendTStabSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rTrendR0(rTrendR0SEXP);
    Rcpp::traits::input_parameter< int >::type timeSteps(timeStepsSEXP);
    rcpp_result_gen = Rcpp::wrap(runModel(basePop, ageGroupsSpan, timeArtStart, entrantPrev, vertTransLag, paedSurveyLag, populationAdjust, entrantPop, birthLag, cumSurv, cumNetMigr, netMigrHivProb, paedSurvCd4Distrib, entrantArtCoverage, paedSurvArtCd4Distrib, survRate, netMigr, asfRate, sexRatioBirth, hivStepsPerYear, cd4Prog, cd4InitDist, cd4Mort, incrrAges, relinfectArt, iota, incrrSex, incidMod, eppMod, scaleCd4Mort, projSteps, tsEpidemicStart, artCd4EligId, specPopPercentElig, pregnantWomenArtElig, who34PercentElig, rSplineRVec, rTrendBeta, rTrendTStab, rTrendR0, timeSteps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cpptest_add", (DL_FUNC) &_cpptest_add, 2},
    {"_cpptest_do_double", (DL_FUNC) &_cpptest_do_double, 1},
    {"_cpptest_modify_std_vector", (DL_FUNC) &_cpptest_modify_std_vector, 1},
    {"_cpptest_get_array", (DL_FUNC) &_cpptest_get_array, 0},
    {"_cpptest_push_array", (DL_FUNC) &_cpptest_push_array, 1},
    {"_cpptest_push_list_arrays", (DL_FUNC) &_cpptest_push_list_arrays, 1},
    {"_cpptest_runModel", (DL_FUNC) &_cpptest_runModel, 41},
    {NULL, NULL, 0}
};

RcppExport void R_init_cpptest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
