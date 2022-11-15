// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   qaKFParticle.cxx
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, CERN
/// \brief  Task to test the performance of the KFParticle package
///

/// includes O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include <CCDB/BasicCCDBManager.h>
/// includes O2Physics
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "TableHelper.h"

/// includes KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

#include <KFParticle.h>
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct qaKFParticle {

  /// general steering settings
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double magneticField = 0.; 

  /// Histogram Configurables
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  /// singe track selections
  Configurable<float> nSigmaTPCMinPi{"nSigmaTPCMinPi", -3., "min number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMaxPi{"nSigmaTPCMaxPi", 3., "max number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMinKa{"nSigmaTPCMinKa", -3., "min number of sigma in the TPC for kaon tracks"};
  Configurable<float> nSigmaTPCMaxKa{"nSigmaTPCMaxKa", 3., "max number of sigma in the TPC for kaon tracks"};
  Configurable<float> pTMin{"pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> etaRange{"etaRange", 0.8, "eta Range for tracks"};
  /// D0 selections
  Configurable<float> pTMinD0{"pTMinD0", 0., "minimum momentum for D0 candidates"};
  Configurable<float> pTMaxD0{"pTMaxD0", 36., "maximum momentum for D0 candidates"};
  Configurable<float> massMinD0{"massMinD0", 1.65, "minimum mass for D0"};
  Configurable<float> massMaxD0{"massMaxD0", 2.08, "minimum mass for D0"};
  Configurable<float> d_cosPA{"d_cosPA", -1., "minimum cosine Pointing angle for D0"};
  Configurable<float> d_decayLength{"d_decayLength", 0., "minimum decay length for D0"};

  // Define which track selection should be used:
  // 0 -> No track selection is applied
  // 1 kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
  //        kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF |
  //                         kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits
  //        kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz
  //        kInAcceptanceTracks = kPtRange | kEtaRange
  // 2 kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
  // 3 kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
  // 4 kQualityTracks
  // 5 kInAcceptanceTracks
  Filter trackFilter = (trackSelection.node() == 0) ||
                      ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                      ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                      ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                      ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                      ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa,aod::pidTOFFullPi, aod::pidTOFFullKa>;
  Partition<TrackTableData> tracksFiltered = (trackSelection.node() == 0) ||
                      ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                      ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                      ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                      ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                      ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
 
  HistogramRegistry histos;

  void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
              o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
              bool isRun3)
  {

    if (mRunNumber != bc.runNumber()) {

      LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
      if (!isRun3) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
      }
      mRunNumber = bc.runNumber();
    }
  } /// end initMagneticFieldCCDB


  void init(InitContext const&)
  {
    if (!doprocessData && !doprocessMC) {
      LOGF(info, "No enabled QA, all histograms are disabled");
      return;
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;

    const AxisSpec axisVertexPosX{500, -1., 1., "X [cm]"};
    const AxisSpec axisVertexPosY{500, -1., 1., "Y [cm]"};
    const AxisSpec axisVertexPosZ{100, -20., 20., "Z [cm]"};
    const AxisSpec axisVertexNumContrib{200, 0, 200, "Number Of contributors to the PV"};
    const AxisSpec axisVertexCov{100, -0.005, 0.005};

    const AxisSpec axisParX{300, -0.5, 0.5, "#it{x} [cm]"};
    const AxisSpec axisParY{200, -0.5, 0.5, "#it{y} [cm]"};
    const AxisSpec axisParZ{200, -11., 11., "#it{z} [cm]"};
    const AxisSpec axisParPX{binsPt, "#it{p}_{x} [GeV/c]"};
    const AxisSpec axisParPY{binsPt, "#it{p}_{y} [GeV/c]"};
    const AxisSpec axisParPZ{binsPt, "#it{p}_{z} [GeV/c]"};

    /// collisions
    histos.add("Events/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

    histos.add("EventsKF/posX", "", kTH1D, {axisVertexPosX});
    histos.add("EventsKF/posY", "", kTH1D, {axisVertexPosY});
    histos.add("EventsKF/posZ", "", kTH1D, {axisVertexPosZ});
    histos.add("EventsKF/posXY", "", kTH2D, {axisVertexPosX, axisVertexPosY});
    histos.add("EventsKF/nContrib", "", kTH1D, {axisVertexNumContrib});
    histos.add("EventsKF/vertexChi2", ";#chi^{2}", kTH1D, {{100, 0, 100}});
    histos.add("EventsKF/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

    /// tracks
    histos.add("TracksKFPi/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKFPi/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKFPi/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKFPi/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKFPi/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKFPi/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKFPi/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, 0.8, 1.2}});
    histos.add("TracksKFPi/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFPi/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});
    // Add nsigma TPC

    histos.add("TracksKFKa/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKFKa/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKFKa/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKFKa/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKFKa/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKFKa/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKFKa/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, 0.8, 1.2}});
    histos.add("TracksKFKa/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFKa/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});
    // Add nsigma TPC

    /// D0 candidates
    histos.add("DZeroCand/atProductionVertex", "at production vertex;", kTH1D, {{2, 0, 1}});
    histos.add("DZeroCand/X", "X [cm]", kTH1D, {axisParX});
    histos.add("DZeroCand/Y", "Y [cm]", kTH1D, {axisParY});
    histos.add("DZeroCand/Z", "Z [cm]", kTH1D, {axisParZ});
    histos.add("DZeroCand/E", "E", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCand/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCand/NDF", "NDF", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCand/p", "momentum", kTH1D, {axisParPX});
    histos.add("DZeroCand/pt", "transverse momentum", kTH1D, {axisParPX});
    histos.add("DZeroCand/eta", "eta", kTH1D, {{100, -2., 2.}});
    histos.add("DZeroCand/phi", "phi", kTH1D, {{100, 0., 3.6}});
    histos.add("DZeroCand/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
    histos.add("DZeroCand/massvspt", "mass vs pt", kTH2D, {{axisParPX},{430, 1.65, 2.08}});
    histos.add("DZeroCand/decayLength", "decay length [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCand/decayLengthXY", "decay length in xy plane [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCand/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
    histos.add("DZeroCand/lifetime", "life time", kTH1D, {{100, 0., 0.2}});
    histos.add("DZeroCand/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
    histos.add("DZeroCand/decayLengthErr", "decay length error [cm]", kTH1D, {{200, 0., 0.2}});
    histos.add("DZeroCand/distToPV", "distance to PV", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCand/deviationToPV", "deviation to PV", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCand/distToPVXY", "distance to PV in xy plane", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCand/deviationToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 20.}});


  }

  /// Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    return true;
  }

  float CosPointingAngleFromKF(KFParticle kfp, KFParticle kfpmother)
  {
    float v[3];
    v[0] = kfp.GetX() - kfpmother.GetX();
    v[1] = kfp.GetY() - kfpmother.GetY();
    v[2] = kfp.GetZ() - kfpmother.GetZ();

    float p[3];
    p[0] = kfp.GetPx();
    p[1] = kfp.GetPy();
    p[2] = kfp.GetPz();

    float ptimesv2 = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    if ( ptimesv2<=0 ) return 0.;
    else {
      double cos = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2]) / sqrt(ptimesv2);
      if(cos >  1.0) cos =  1.0;
      if(cos < -1.0) cos = -1.0;
      return cos;
    }
  }

  /// Process function for data
  void processData(CollisionTableData const& collisions, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
      /// Set magnetic field for KF vertexing
      KFParticle::SetField(magneticField);
    }
    for (auto const& collision : collisions) {
      /// Apply event selection
      if (!isSelectedCollision(collision)) {
        continue;
      }
      auto tracks = tracksFiltered->sliceByCached(aod::track::collisionId, collision.globalIndex());
      /// set KF primary vertex
      KFPVertex kfpVertex;
      kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
      /// CAREFUL!!!!!! Covariance matrix elements yy and xz are switched until a central fix is provided!!!!!
      kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ());
      kfpVertex.SetChi2(collision.chi2());
      kfpVertex.SetNDF(2); //?? What number should I put here?
      kfpVertex.SetNContributors(collision.numContrib());

      KFParticle KFPV(kfpVertex);

      /// fill collision parameters
      histos.fill(HIST("Events/covXX"), collision.covXX());
      histos.fill(HIST("Events/covXY"), collision.covXY());
      histos.fill(HIST("Events/covXZ"), collision.covXZ());
      histos.fill(HIST("Events/covYY"), collision.covYY());
      histos.fill(HIST("Events/covYZ"), collision.covYZ());
      histos.fill(HIST("Events/covZZ"), collision.covZZ());

      histos.fill(HIST("EventsKF/posX"), kfpVertex.GetX());
      histos.fill(HIST("EventsKF/posY"), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/posZ"), kfpVertex.GetZ());
      histos.fill(HIST("EventsKF/posXY"), kfpVertex.GetX(), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/nContrib"), kfpVertex.GetNContributors());
      histos.fill(HIST("EventsKF/vertexChi2"), kfpVertex.GetChi2());
      histos.fill(HIST("EventsKF/covXX"), kfpVertex.GetCovariance(0));
      histos.fill(HIST("EventsKF/covXY"), kfpVertex.GetCovariance(1));
      histos.fill(HIST("EventsKF/covYY"), kfpVertex.GetCovariance(2));
      histos.fill(HIST("EventsKF/covXZ"), kfpVertex.GetCovariance(3));
      histos.fill(HIST("EventsKF/covYZ"), kfpVertex.GetCovariance(4));
      histos.fill(HIST("EventsKF/covZZ"), kfpVertex.GetCovariance(5));

      for (auto& [track1, track2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

        auto track1p = track1.p();
        auto track2p = track2.p();

        o2::track::TrackParametrizationWithError trackparCovPi;
        o2::track::TrackParametrizationWithError trackparCovKa; 

        bool CandD0 = false;
        bool CandD0bar = false;
        /// At the moment pT independent TPC selection. Add a minimum and Maximum momentum
        /// Apply TPC+TOF at higher momenta (TOF still uncalibrated for LHC22f).
        /// Select D0 and D0bar candidates
        if (nSigmaTPCMinPi <= track1.tpcNSigmaPi() && track1.tpcNSigmaPi() <= nSigmaTPCMaxPi && nSigmaTPCMinKa <= track2.tpcNSigmaKa() && track2.tpcNSigmaPi() <= nSigmaTPCMaxKa) {
          if (track1.sign() == 1 && track2.sign() == -1) {
            CandD0 = true;
            trackparCovPi = getTrackParCov(track1);
            trackparCovKa = getTrackParCov(track2);
          }
          else if (track1.sign() == -1 && track2.sign() == 1) {
            CandD0bar = true;
            trackparCovPi = getTrackParCov(track1);
            trackparCovKa = getTrackParCov(track2);
          }
          else {
            continue;
          }
        }
        else if(nSigmaTPCMinKa <= track1.tpcNSigmaKa() && track1.tpcNSigmaKa() <= nSigmaTPCMaxKa && nSigmaTPCMinPi <= track2.tpcNSigmaPi() && track2.tpcNSigmaPi() <= nSigmaTPCMaxPi) {
          if (track1.sign() == 1 && track2.sign() == -1) {
            CandD0bar = true;
            trackparCovPi = getTrackParCov(track2);
            trackparCovKa = getTrackParCov(track1);
          }
          else if (track1.sign() == -1 && track2.sign() == 1) {
            CandD0 = true;
            trackparCovPi = getTrackParCov(track2);
            trackparCovKa = getTrackParCov(track1);
          }
          else {
            continue;
          }
        }
        else {
          continue;
        }

        /// Apply single track cuts as a prefilter on the daughter tracks.
        if (track1p < pTMin || abs(track1.eta())>etaRange || track2p < pTMin || abs(track2.eta())>etaRange) {
          continue;
        }
      

        /// Keep in mind that in the track table the parameters are stored after propagation to the DCA to the PV.
        /// Check if we need to get another table. Or where the parameters might have to be propagated.
        array<float, 3> trkpos_parPi;
        array<float, 3> trkpos_parKa;
        array<float, 3> trkmom_parPi;
        array<float, 3> trkmom_parKa;
        array<float, 21> trk_covPi;
        array<float, 21> trk_covKa;


        trackparCovPi.getXYZGlo(trkpos_parPi);
        trackparCovPi.getPxPyPzGlo(trkmom_parPi);
        trackparCovPi.getCovXYZPxPyPzGlo(trk_covPi);
        trackparCovKa.getXYZGlo(trkpos_parKa);
        trackparCovKa.getPxPyPzGlo(trkmom_parKa);
        trackparCovKa.getCovXYZPxPyPzGlo(trk_covKa);
        float trkpar_KFPi[6] = {trkpos_parPi[0], trkpos_parPi[1], trkpos_parPi[2],
                                   trkmom_parPi[0], trkmom_parPi[1], trkmom_parPi[2]};
        float trkpar_KFKa[6] = {trkpos_parKa[0], trkpos_parKa[1], trkpos_parKa[2],
                                   trkmom_parKa[0], trkmom_parKa[1], trkmom_parKa[2]};
        float trkcov_KFPi[21];
        float trkcov_KFKa[21];
        for (int i = 0; i < 21; i++) {
          trkcov_KFPi[i] = trk_covPi[i];
          trkcov_KFKa[i] = trk_covKa[i];
        }

        KFPTrack kfpTrackPi;
        KFPTrack kfpTrackKa;
        kfpTrackPi.SetParameters(trkpar_KFPi);
        kfpTrackPi.SetCovarianceMatrix(trkcov_KFPi);
        kfpTrackKa.SetParameters(trkpar_KFKa);
        kfpTrackKa.SetCovarianceMatrix(trkcov_KFKa);

        if (CandD0) {
          kfpTrackPi.SetCharge(1);
          kfpTrackKa.SetCharge(-1);
        } 
        else if (CandD0bar) {
          kfpTrackPi.SetCharge(-1);
          kfpTrackKa.SetCharge(1);
        }
        /// Add these quantities!
        kfpTrackPi.SetNDF(1); // Which is the correct number?
        kfpTrackKa.SetNDF(1); // Which is the correct number?
        //kfpTrack.SetChi2(...);

        KFParticle KFPion(kfpTrackPi, 211);
        KFParticle KFKaon(kfpTrackKa, 321);



        /// fill track parameters
        histos.fill(HIST("TracksKFPi/x"), kfpTrackPi.GetX());
        histos.fill(HIST("TracksKFPi/y"), kfpTrackPi.GetY());
        histos.fill(HIST("TracksKFPi/z"), kfpTrackPi.GetZ());
        histos.fill(HIST("TracksKFPi/px"), kfpTrackPi.GetPx());
        histos.fill(HIST("TracksKFPi/py"), kfpTrackPi.GetPy());
        histos.fill(HIST("TracksKFPi/pz"), kfpTrackPi.GetPz());
        histos.fill(HIST("TracksKFPi/chi2perNDF"), kfpTrackPi.GetChi2perNDF());
        histos.fill(HIST("TracksKFPi/dcaXY"), KFPion.GetDistanceFromVertexXY(kfpVertex));
        histos.fill(HIST("TracksKFPi/length"), KFPion.GetDecayLength());

        histos.fill(HIST("TracksKFKa/x"), kfpTrackKa.GetX());
        histos.fill(HIST("TracksKFKa/y"), kfpTrackKa.GetY());
        histos.fill(HIST("TracksKFKa/z"), kfpTrackKa.GetZ());
        histos.fill(HIST("TracksKFKa/px"), kfpTrackKa.GetPx());
        histos.fill(HIST("TracksKFKa/py"), kfpTrackKa.GetPy());
        histos.fill(HIST("TracksKFKa/pz"), kfpTrackKa.GetPz());
        histos.fill(HIST("TracksKFKa/chi2perNDF"), kfpTrackKa.GetChi2perNDF());
        histos.fill(HIST("TracksKFPi/dcaXY"), KFKaon.GetDistanceFromVertexXY(kfpVertex));
        histos.fill(HIST("TracksKFPi/length"), KFKaon.GetDecayLength());

        KFParticle KFDZero;
        const KFParticle *D0Daughters[2] = {&KFPion, &KFKaon};
        int NDaughters = 2;
        KFDZero.SetConstructMethod(2);
        KFDZero.Construct(D0Daughters, NDaughters, &KFPV);

        float X=0., Y=0.,Z=0., E=0., Chi2=0., NDF=0., P=0., Pt=0., Eta=0., Phi=0., mass=0., decayLength=0., decayLengthxy=0., cosPA=-1., lifeTime=0., massErr=0., decayLengthErr=0.;
        bool atProductionVertex = false;
        float distToPV=0., deviationToPV=0., distToPVxy=0., deviationToPVxy=0.;

        X=KFDZero.GetX();
        Y=KFDZero.GetY();
        Z=KFDZero.GetZ();
        E=KFDZero.GetE();
        Chi2=KFDZero.GetChi2();
        NDF=KFDZero.GetNDF();
        P=KFDZero.GetP();
        Pt=KFDZero.GetPt();
        Eta=KFDZero.GetEta();
        Phi=KFDZero.GetPhi();
        mass=KFDZero.GetMass();
        decayLength=KFDZero.GetDecayLength();
        decayLengthxy=KFDZero.GetDecayLengthXY();
        cosPA=CosPointingAngleFromKF(KFPV, KFDZero);
        lifeTime=KFDZero.GetLifeTime();
        massErr=KFDZero.GetErrMass();
        decayLengthErr=KFDZero.GetErrDecayLength();
        distToPV=KFDZero.GetDistanceFromVertex(KFPV);
        deviationToPV=KFDZero.GetDeviationFromVertex(KFPV);
        distToPVxy=KFDZero.GetDistanceFromVertexXY(KFPV);
        deviationToPVxy=KFDZero.GetDeviationFromVertexXY(KFPV);
        atProductionVertex=KFDZero.GetAtProductionVertex();

        /// Remove daughter tracks from PV fit?
        /// Pt selection
        if (Pt<pTMinD0 || Pt>pTMaxD0) {
          continue;
        }
        /// Mass window selection
        if (mass<massMinD0 || mass>massMaxD0) {
          continue;
        }
        /// cosine pointing anle selection
        if (cosPA < d_cosPA) {
          continue;
        }
        /// decay length selection
        if (decayLength<d_decayLength) {
          continue;
        }
        
        histos.fill(HIST("DZeroCand/atProductionVertex"), atProductionVertex);
        histos.fill(HIST("DZeroCand/X"), X);
        histos.fill(HIST("DZeroCand/Y"), Y);
        histos.fill(HIST("DZeroCand/Z"), Z);
        histos.fill(HIST("DZeroCand/E"), E);
        histos.fill(HIST("DZeroCand/Chi2"), Chi2);
        histos.fill(HIST("DZeroCand/NDF"), NDF);
        histos.fill(HIST("DZeroCand/p"), P);
        histos.fill(HIST("DZeroCand/pt"), Pt);
        histos.fill(HIST("DZeroCand/eta"), Eta);
        histos.fill(HIST("DZeroCand/phi"), Phi);
        histos.fill(HIST("DZeroCand/mass"), mass);
        histos.fill(HIST("DZeroCand/massvspt"), Pt, mass);
        histos.fill(HIST("DZeroCand/decayLength"), decayLength);
        histos.fill(HIST("DZeroCand/decayLengthXY"), decayLengthxy);
        histos.fill(HIST("DZeroCand/cosPA"), cosPA);
        histos.fill(HIST("DZeroCand/lifetime"), lifeTime);
        histos.fill(HIST("DZeroCand/massErr"), massErr);
        histos.fill(HIST("DZeroCand/massErr"), massErr);
        histos.fill(HIST("DZeroCand/decayLengthErr"), decayLengthErr);
        histos.fill(HIST("DZeroCand/decayLengthErr"), decayLengthErr);
        histos.fill(HIST("DZeroCand/distToPV"), distToPV);
        histos.fill(HIST("DZeroCand/deviationToPV"), deviationToPV);
        histos.fill(HIST("DZeroCand/distToPVXY"), distToPVxy);
        histos.fill(HIST("DZeroCand/deviationToPVXY"), deviationToPVxy);
        


        /// Add the secondary vertex and quantities of daughter particles
        // KFDZero.TransportToDecayVertex();
        // float distToPVPi=0., distToPVKa=0. dcaBetweenTracks=0.;
        // distToPVPi = KFPion.GetDistanceFromVertex(KFPV);
        // distToPVKa = KFKaon.GetDistanceFromVertex(KFPV);
        // dcaBetweenTracks = KFPion.GetDistanceFromParticle(KFKaon);


      }


    }
  }
  PROCESS_SWITCH(qaKFParticle, processData, "process data", true);

  /// Process function for MC
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  void processMC(CollisionTableMC const& collisions)
  {
    for (auto const& collision : collisions) {
    }
  }
  PROCESS_SWITCH(qaKFParticle, processMC, "process mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaKFParticle>(cfgc)};
}
