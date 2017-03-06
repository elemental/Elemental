/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
El::LPInfo<Real> DenseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }


    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;

    timer.Start();
    auto info = El::LP( problem, solution, ctrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output
    ("primal objective:       ",info.ipmInfo.primalObjective);
    El::Output
    ("dual   objective:       ",info.ipmInfo.dualObjective);
    El::Output
    ("infeasibility:          ",info.ipmInfo.infeasibilityError);
    El::Output
    ("relative objective gap: ",info.ipmInfo.relativeObjectiveGap);
    El::Output
    ("relative comp. gap:     ",info.ipmInfo.relativeComplementarityGap);
    El::Output
    ("num iterations:         ",info.ipmInfo.numIterations);
    return info;
}

template<typename Real>
El::LPInfo<Real> SparseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    // We will wait to do any presolves until a fast mechanism exists for
    // deleting entries from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;
    ctrl.ipmCtrl.zMinPivotValueLogEps = Real(2.0);

    timer.Start();
    auto info = El::LP( problem, solution, ctrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output
    ("primal objective:       ",info.ipmInfo.primalObjective);
    El::Output
    ("dual   objective:       ",info.ipmInfo.dualObjective);
    El::Output
    ("infeasibility:          ",info.ipmInfo.infeasibilityError);
    El::Output
    ("relative objective gap: ",info.ipmInfo.relativeObjectiveGap);
    El::Output
    ("relative comp. gap:     ",info.ipmInfo.relativeComplementarityGap);
    El::Output
    ("num iterations:         ",info.ipmInfo.numIterations);
    return info;
}

template<typename Real>
void SparseNetlibLPData
( const std::string& directory,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    El::Output
    ("Will test netlib LP_data suite at ",directory,
     " using El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    std::vector<std::pair<std::string,Real>> problemObjectives; 
    // The netlib-documented value for 25fv47 is +5.5018458883e+03.
    problemObjectives.emplace_back( "25fv47",
      Real(+5.501845888286744794581232588391644178145257e+03) );
    // TODO(poulson): Extend the remaining values with the QuadDouble output.
    problemObjectives.emplace_back( "80bau3b",  Real(+9.8723216072e+05) );
    problemObjectives.emplace_back( "adlittle", Real(+2.2549496316e+05) );
    problemObjectives.emplace_back( "afiro",    Real(-4.6475314286e+02) );
    problemObjectives.emplace_back( "agg",      Real(-3.5991767287e+07) );
    problemObjectives.emplace_back( "agg2",     Real(-2.0239252356e+07) );
    problemObjectives.emplace_back( "agg3",     Real(+1.0312115935e+07) );
    problemObjectives.emplace_back( "bandm",    Real(-1.5862801845e+02) );
    problemObjectives.emplace_back( "beaconfd", Real(+3.3592485807e+04) );
    problemObjectives.emplace_back( "blend",    Real(-3.0812149846e+01) );
    problemObjectives.emplace_back( "bnl1",     Real(+1.9776292856e+03) );
    problemObjectives.emplace_back( "bnl2",     Real(+1.8112365404e+03) );
    problemObjectives.emplace_back( "boeing1",  Real(-3.3521356751e+02) );
    problemObjectives.emplace_back( "boeing2",  Real(-3.1501872802e+02) );
    problemObjectives.emplace_back( "bore3d",   Real(+1.3730803942e+03) );
    problemObjectives.emplace_back( "brandy",   Real(+1.5185098965e+03) );
    problemObjectives.emplace_back( "capri",    Real(+2.6900129138e+03) );
    problemObjectives.emplace_back( "cycle",    Real(-5.2263930249e+00) );
    problemObjectives.emplace_back( "czprob",   Real(+2.1851966989e+06) );
    problemObjectives.emplace_back( "d2q06c",   Real(+1.2278423615e+05) );
    problemObjectives.emplace_back( "d6cube",   Real(+3.1549166667e+02) );
    problemObjectives.emplace_back( "degen2",   Real(-1.4351780000e+03) );
    problemObjectives.emplace_back( "degen3",   Real(-9.8729400000e+02) );
    problemObjectives.emplace_back( "dfl001",   Real(+1.126639607e+07) );
    problemObjectives.emplace_back( "e226",     Real(-1.8751929066e+01) );
    problemObjectives.emplace_back( "etamacro", Real(-7.5571521774e+02) );
    problemObjectives.emplace_back( "fffff800", Real(+5.5567961165e+05) );
    problemObjectives.emplace_back( "finnis",   Real(+1.7279096547e+05) );
    problemObjectives.emplace_back( "fit1d",    Real(-9.1463780924e+03) );
    problemObjectives.emplace_back( "fit1p",    Real(+9.1463780924e+03) );
    problemObjectives.emplace_back( "fit2d",    Real(-6.8464293294e+04) );
    problemObjectives.emplace_back( "fit2p",    Real(+6.8464293232e+04) );
    problemObjectives.emplace_back( "forplan",  Real(-6.6421873953e+02) );
    problemObjectives.emplace_back( "ganges",   Real(-1.0958636356e+05) );
    problemObjectives.emplace_back( "gfrd-pnc", Real(+6.9022359995e+06) );
    problemObjectives.emplace_back( "greenbea", Real(-7.2462405908e+07) );
    problemObjectives.emplace_back( "greenbeb", Real(-4.3021476065e+06) );
    problemObjectives.emplace_back( "grow15",   Real(-1.0687094129e+08) );
    problemObjectives.emplace_back( "grow22",   Real(-1.6083433648e+08) );
    problemObjectives.emplace_back( "grow7",    Real(-4.7787811815e+07) );
    problemObjectives.emplace_back( "israel",   Real(-8.9664482186e+05) );
    problemObjectives.emplace_back( "kb2",      Real(-1.7499001299e+03) );
    problemObjectives.emplace_back( "lotfi",    Real(-2.5264706062e+01) );
    problemObjectives.emplace_back( "maros",    Real(-5.8063743701e+04) );
    problemObjectives.emplace_back( "maros-r7", Real(+1.4971851665e+06) );
    problemObjectives.emplace_back( "modszk1",  Real(+3.2061972906e+02) );
    problemObjectives.emplace_back( "nesm",     Real(+1.4076073035e+07) );
    problemObjectives.emplace_back( "perold",   Real(-9.3807580773e+03) );
    problemObjectives.emplace_back( "pilot",    Real(-5.5740430007e+02) );
    problemObjectives.emplace_back( "pilot.ja", Real(-6.1131344111e+03) );
    problemObjectives.emplace_back( "pilot.we", Real(-2.7201027439e+06) );
    problemObjectives.emplace_back( "pilot4",   Real(-2.5811392641e+03) );
    problemObjectives.emplace_back( "pilot87",  Real(+3.0171072827e+02) );
    problemObjectives.emplace_back( "pilotnov", Real(-4.4972761882e+03) );
    /*
    problemObjectives.emplace_back( "qap8",  Real(+2.0350000000e+02) );
    problemObjectives.emplace_back( "qap12", Real(+5.2289435056e+02) );
    problemObjectives.emplace_back( "qap15", Real(+1.0409940410e+03) );
    */
    problemObjectives.emplace_back( "recipe",   Real(-2.6661600000e+02) );
    problemObjectives.emplace_back( "sc105",    Real(-5.2202061212e+01) );
    problemObjectives.emplace_back( "sc205",    Real(-5.2202061212e+01) );
    problemObjectives.emplace_back( "sc50a",    Real(-6.4575077059e+01) );
    problemObjectives.emplace_back( "sc50b",    Real(-7.0000000000e+01) );
    problemObjectives.emplace_back( "scagr25",  Real(-1.4753433061e+07) );
    problemObjectives.emplace_back( "scagr7",   Real(-2.3313892548e+06) );
    problemObjectives.emplace_back( "scfxm1",   Real(+1.8416759028e+04) );
    problemObjectives.emplace_back( "scfxm2",   Real(+3.6660261565e+04) );
    problemObjectives.emplace_back( "scfxm3",   Real(+5.4901254550e+04) );
    problemObjectives.emplace_back( "scorpion", Real(+1.8781248227e+03) );
    problemObjectives.emplace_back( "scrs8",    Real(+9.0429998619e+02) );
    problemObjectives.emplace_back( "scsd1",    Real(+8.6666666743e+00) );
    problemObjectives.emplace_back( "scsd6",    Real(+5.0500000078e+01) );
    problemObjectives.emplace_back( "scsd8",    Real(+9.0499999993e+02) );
    problemObjectives.emplace_back( "sctap1",   Real(+1.4122500000e+03) );
    problemObjectives.emplace_back( "sctap2",   Real(+1.7248071429e+03) );
    problemObjectives.emplace_back( "sctap3",   Real(+1.4240000000e+03) );
    problemObjectives.emplace_back( "seba",     Real(+1.5711600000e+04) );
    problemObjectives.emplace_back( "share1b",  Real(-7.6589318579e+04) );
    problemObjectives.emplace_back( "share2b",  Real(-4.1573224074e+02) );
    problemObjectives.emplace_back( "shell",    Real(+1.2088253460e+09) );
    problemObjectives.emplace_back( "ship04l",  Real(+1.7933245380e+06) );
    problemObjectives.emplace_back( "ship04s",  Real(+1.7987147004e+06) );
    problemObjectives.emplace_back( "ship08l",  Real(+1.9090552114e+06) );
    problemObjectives.emplace_back( "ship08s",  Real(+1.9200982105e+06) );
    problemObjectives.emplace_back( "ship12l",  Real(+1.4701879193e+06) );
    problemObjectives.emplace_back( "ship12s",  Real(+1.4892361344e+06) );
    problemObjectives.emplace_back( "sierra",   Real(+1.5394362184e+07) );
    problemObjectives.emplace_back( "stair",    Real(-2.5126695119e+02) );
    problemObjectives.emplace_back( "standata", Real(+1.2576995000e+03) );
    /* Skipping 'standgub' */
    problemObjectives.emplace_back( "standmps", Real(+1.4060175000e+03) );
    problemObjectives.emplace_back( "stocfor1", Real(-4.1131976219e+04) );
    problemObjectives.emplace_back( "stocfor2", Real(-3.9024408538e+04) );
    /*
    problemObjectives.emplace_back( "stocfor3", Real(-3.9976661576e+04) );
    problemObjectives.emplace_back( "truss",    Real(+4.5881584719e+05) );
    */
    problemObjectives.emplace_back( "tuff",     Real(+2.9214776509e-01) );
    problemObjectives.emplace_back( "vtp.base", Real(+1.2983146246e+05) );
    problemObjectives.emplace_back( "wood1p",   Real(+1.4429024116e+00) );
    problemObjectives.emplace_back( "woodw",    Real(+1.3044763331e+00) );

    const Real epsilon = El::limits::Epsilon<Real>();
    const Real demandedRelativeError =
      El::Pow( epsilon, relativeObjectiveGapTolLogEps );
    std::vector<std::pair<Real,std::string>> problemRelErrors;
    std::vector<std::pair<El::Int,std::string>> problemIterations;
    for( const auto& problemObjective : problemObjectives )
    {
        const std::string& name = problemObjective.first;
        const Real& objective = problemObjective.second;
        El::Output("Testing ",name);
        const std::string filename = directory + "/" + name + ".mps";
        auto info =
          SparseLoadAndSolve<Real>
          ( filename, metadataSummary, print, progress, outerEquil,
            infeasibilityTolLogEps,
            relativeObjectiveGapTolLogEps,
            relativeComplementarityGapTolLogEps,
            xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
            xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
            lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
            forceSameStep, compositeNewton );
        const Real relativeError =
          El::Abs(objective - info.ipmInfo.primalObjective) /
          El::Abs(objective);
        if( relativeError > demandedRelativeError )
        {
            El::Output
            ("WARNING: Only solved ",name," to objective ",
             info.ipmInfo.primalObjective," (vs. ",objective,")");
        }
       problemRelErrors.emplace_back( relativeError, name );
       problemIterations.emplace_back( info.ipmInfo.numIterations, name );
    }
    El::Output("");

    const El::Int numWorstErrors = 10;
    std::sort( problemRelErrors.rbegin(), problemRelErrors.rend() );
    El::Output(numWorstErrors," problems with highest error:");
    for( El::Int i=0; i<El::Min(numWorstErrors,problemRelErrors.size()); ++i )
    {
        const auto& problemRelError = problemRelErrors[i];
        El::Output
        (problemRelError.second,
         " was only solved to relative accuracy of ",problemRelError.first);
    }
    El::Output("");

    const El::Int numMostIters = 10;
    std::sort( problemIterations.rbegin(), problemIterations.rend() );
    El::Output(numMostIters," problems taking the most iterations:");
    for( El::Int i=0; i<El::Min(numMostIters,problemIterations.size()); ++i )
    {
        const auto& problemIteration = problemIterations[i];
        El::Output
        (problemIteration.second," took ",problemIteration.first," iterations");
    }
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const std::string filename =
          El::Input
          ("--filename","MPS filename",
           std::string("../data/optimization/lp_data/share1b.mps"));
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize MPS metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool testArbitrary =
          El::Input("--testArbitrary","test arbitrary precision?",false);
        const bool print = El::Input("--print","print matrices?",false);
        const bool progress = El::Input("--progress","IPM progress?",true);
        const bool outerEquil =
          El::Input("--outerEquil","outer equilibration?",true);
        const double infeasibilityTolLogEps =
          El::Input
          ("--infeasibilityTolLogEps","log_eps(infeasibilityTol)",0.45);
        const double relativeObjectiveGapTolLogEps =
          El::Input
          ("--relativeObjectiveGapTolLogEps",
           "log_eps(relativeObjectiveGapTol)",0.05);
        const double relativeComplementarityGapTolLogEps =
          El::Input
          ("--relativeComplementarityGapTolLogEps",
           "log_eps(relativeComplementarityGapTol)",0.3);
        const double xRegSmallLogEps =
          El::Input("--xRegSmallLogEps","log_eps(xRegSmall)",0.8);
        const double yRegSmallLogEps =
          El::Input("--yRegSmallLogEps","log_eps(yRegSmall)",0.8);
        const double zRegSmallLogEps =
          El::Input("--zRegSmallLogEps","log_eps(zRegSmall)",0.8);
        const double xRegLargeLogEps =
          El::Input("--xRegLargeLogEps","log_eps(xRegLarge)",0.6);
        const double yRegLargeLogEps =
          El::Input("--yRegLargeLogEps","log_eps(yRegLarge)",0.6);
        const double zRegLargeLogEps =
          El::Input("--zRegLargeLogEps","log_eps(zRegLarge)",0.6);
        const double lowerTargetRatioLogCompRatio =
          El::Input
          ("--lowerTargetRatioLogCompRatio","log_compratio(lowerTargetRatio)",
           -0.25);
        const double upperTargetRatioLogCompRatio =
          El::Input
          ("--upperTargetRatioLogCompRatio","log_compratio(upperTargetRatio)",
           0.25);
        const bool forceSameStep =
          El::Input
          ("--forceSameStep","force same primal and dual step sizes?",false);
        const bool compositeNewton =
          El::Input("--compositeNewton","Mehrotra predictor-corrector?",false);
        const bool testNetlib =
          El::Input("--testNetlib","test entire netlib LP_data?",false);
        const std::string netlibDirectory =
          El::Input
          ("--netlibDirectory","path to netlib LP_data MPS files",
           std::string("../data/optimization/lp_data"));
        El::ProcessInput();
        El::PrintInputReport();

        if( testNetlib )
        {
            // TODO(poulson): Provide a command-line interface for testing the
            // netlib suite with other datatypes.
            SparseNetlibLPData<double>
            ( netlibDirectory, metadataSummary, print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
            return 0;
        }

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, metadataSummary,
                  print, progress, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
                  forceSameStep, compositeNewton );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
#endif
#ifdef EL_HAVE_QUAD
            DenseLoadAndSolve<El::Quad>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
#endif
#ifdef EL_HAVE_MPC
            if( testArbitrary )
            {
                // TODO(poulson): Make this configurable.
                El::mpfr::SetPrecision( 512 );
                DenseLoadAndSolve<El::BigFloat>
                ( filename, metadataSummary,
                  print, progress, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
                  forceSameStep, compositeNewton );
            }
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
#ifdef EL_HAVE_QD
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton );
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton );
#endif
#ifdef EL_HAVE_QUAD
        SparseLoadAndSolve<El::Quad>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton );
#endif
#ifdef EL_HAVE_MPC
        if( testArbitrary )
        {
            // TODO(poulson): Make this configurable
            El::mpfr::SetPrecision( 512 );
            SparseLoadAndSolve<El::BigFloat>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton );
        }
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
