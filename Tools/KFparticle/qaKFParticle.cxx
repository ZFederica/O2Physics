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
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
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
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "TableHelper.h"
#include "Tools/KFparticle/KFUtilities.h"
#include "Tools/KFparticle/qaKFParticle.h"
#include <iostream>
using namespace std;

#ifndef HomogeneousField

#define HomogeneousField

#endif
/// includes KFParticle
#include "KFParticle.h"
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
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 24., 36., 50.0}, ""};

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  /// singe track selections
  Configurable<float> nSigmaTPCMinPi{"nSigmaTPCMinPi", -3., "min number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMaxPi{"nSigmaTPCMaxPi", 3., "max number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMinKa{"nSigmaTPCMinKa", -3., "min number of sigma in the TPC for kaon tracks"};
  Configurable<float> nSigmaTPCMaxKa{"nSigmaTPCMaxKa", 3., "max number of sigma in the TPC for kaon tracks"};
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> d_etaRange{"d_etaRange", 0.8, "eta Range for tracks"};
  Configurable<float> d_dcaXYTrackPV{"d_dcaXYTrackPV", 2., "DCA XY of the daughter tracks to the PV"};
  Configurable<float> d_dcaZTrackPV{"d_dcaZTrackPV", 10., "DCA Z of the daughter tracks to the PV"};
  /// D0 selections
  Configurable<float> d_pTMinD0{"d_pTMinD0", 0., "minimum momentum for D0 candidates"};
  Configurable<float> d_pTMaxD0{"d_pTMaxD0", 36., "maximum momentum for D0 candidates"};
  Configurable<float> d_massMinD0{"d_massMinD0", 1.65, "minimum mass for D0"};
  Configurable<float> d_massMaxD0{"d_massMaxD0", 2.08, "minimum mass for D0"};
  Configurable<float> d_cosPA{"d_cosPA", -1., "minimum cosine Pointing angle for D0"};
  Configurable<float> d_decayLength{"d_decayLength", 0., "minimum decay length for D0"};
  Configurable<float> d_normdecayLength{"d_normdecayLength", 100., "minimum normalised decay length for D0"};
  Configurable<float> d_chi2topoD0{"d_chi2topoD0", 1000., "maximum chi2 topological of D0 to PV"};
  Configurable<float> d_dist3DSVDau{"d_dist3DSVDau", 1000., "maximum geometrical distance 3D daughter tracks at the SV"};
  Configurable<float> d_cosThetaStarPi{"d_cosThetaStarPi", 1000., "maximum cosine theta star of the pion from D0"};
  Configurable<float> d_cosThetaStarKa{"d_cosThetaStarKa", 1000., "maximum cosine theta star of the kaon from D0"};
  Configurable<float> d_distPiToSV{"d_distPiToSV", 1000., "maximum distance Pi to SV"};
  Configurable<float> d_distKaToSV{"d_distKaToSV", 1000., "maximum distance Ka to SV"};
  Configurable<float> d_d0pid0ka{"d_d0pid0ka", -100000., "maximum product of impact parameters of daughters to the PV"};
  /// Option to write D0 variables in a tree
  Configurable<bool> writeTree{"writeTree", false, "write D0 variables in a tree"}; 

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

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;
  // Partition<TrackTableData> tracksFiltered = (trackSelection.node() == 0) ||
  //                     ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
  //                     ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
  //                     ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
  //                     ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
  //                     ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  HistogramRegistry histos;
  /// Table to be produced
  Produces<o2::aod::TreeDZeroKF> rowDZeroTree;

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
    histos.add("TracksKFPi/length", "track length in cm;#it{Length} [cm];", kTH1D, {{100, 0, 20}});
    // Add nsigma TPC

    histos.add("TracksKFKa/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKFKa/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKFKa/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKFKa/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKFKa/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKFKa/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKFKa/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, 0.8, 1.2}});
    histos.add("TracksKFKa/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFKa/length", "track length in cm;#it{Length} [cm];", kTH1D, {{100, 0, 20}});
    // Add nsigma TPC

    /// D0 candidates
    histos.add("DZeroCandGeo/X", "X [cm]", kTH1D, {axisParX});
    histos.add("DZeroCandGeo/Y", "Y [cm]", kTH1D, {axisParY});
    histos.add("DZeroCandGeo/Z", "Z [cm]", kTH1D, {axisParZ});
    histos.add("DZeroCandGeo/E", "E", kTH1D, {{100, 0., 50.}});
    histos.add("DZeroCandGeo/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandGeo/NDF", "NDF", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandGeo/Chi2OverNDF", "NDF", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandGeo/p", "momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandGeo/pt", "transverse momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandGeo/eta", "eta", kTH1D, {{100, -2., 2.}});
    histos.add("DZeroCandGeo/phi", "phi", kTH1D, {{100, 0., 3.6}});
    histos.add("DZeroCandGeo/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
    histos.add("DZeroCandGeo/massvspt", "mass vs pt", kTH2D, {{axisParPX}, {430, 1.65, 2.08}});
    histos.add("DZeroCandGeo/decayLength", "decay length [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCandGeo/decayLengthXY", "decay length in xy plane [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCandGeo/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
    histos.add("DZeroCandGeo/lifetime", "life time", kTH1D, {{100, 0., 0.2}});
    histos.add("DZeroCandGeo/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
    histos.add("DZeroCandGeo/decayLengthErr", "decay length error [cm]", kTH1D, {{200, 0., 0.2}});
    histos.add("DZeroCandGeo/distToPV", "distance to PV", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandGeo/deviationToPV", "deviation to PV", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCandGeo/distToPVXY", "distance to PV in xy plane", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandGeo/deviationToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 20.}});


    histos.add("DZeroCandTopo/X", "X [cm]", kTH1D, {axisParX});
    histos.add("DZeroCandTopo/Y", "Y [cm]", kTH1D, {axisParY});
    histos.add("DZeroCandTopo/Z", "Z [cm]", kTH1D, {axisParZ});
    histos.add("DZeroCandTopo/E", "E", kTH1D, {{100, 0., 50.}});
    histos.add("DZeroCandTopo/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandTopo/NDF", "NDF", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandTopo/Chi2OverNDF", "NDF", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandTopo/p", "momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandTopo/pt", "transverse momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandTopo/eta", "eta", kTH1D, {{100, -2., 2.}});
    histos.add("DZeroCandTopo/phi", "phi", kTH1D, {{100, 0., 3.6}});
    histos.add("DZeroCandTopo/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
    histos.add("DZeroCandTopo/massvspt", "mass vs pt", kTH2D, {{axisParPX}, {430, 1.65, 2.08}});
    histos.add("DZeroCandTopo/decayLength", "decay length [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCandTopo/decayLengthXY", "decay length in xy plane [cm]", kTH1D, {{200, 0., 2.}});
    histos.add("DZeroCandTopo/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
    histos.add("DZeroCandTopo/lifetime", "life time", kTH1D, {{100, 0., 0.2}});
    histos.add("DZeroCandTopo/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
    histos.add("DZeroCandTopo/decayLengthErr", "decay length error [cm]", kTH1D, {{200, 0., 0.2}});
    histos.add("DZeroCandTopo/distToPV", "distance to PV", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandTopo/deviationToPV", "deviation to PV", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCandTopo/distToPVXY", "distance to PV in xy plane", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandTopo/deviationToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCandTopo/deviationDaugtherTracks", "chi2 in 3D of daughter tracks at the SV", kTH1D, {{200, 0., 0.2}});
    histos.add("DZeroCandTopo/distanceDaugtherTracks", "distance in 3D of daughter tracks at the SV", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandTopo/cosThetaStarPion", "cosine theta star of the pion from D0", kTH1D, {{100, -10., 10.}});
    histos.add("DZeroCandTopo/cosThetaStarKaon", "cosine theta star of the kaon from D0", kTH1D, {{100, -10., 10}});
    histos.add("DZeroCandTopo/deviationPiToSV", "deviation of Pi to SV", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCandTopo/distPiToSV", "distance of Pi to SV", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandTopo/deviationKaToSV", "deviation of Ka to SV", kTH1D, {{200, 0., 20.}});
    histos.add("DZeroCandTopo/distKaToSV", "distance of Ka to SV", kTH1D, {{100, 0., 1.}});
    histos.add("DZeroCandTopo/d0pid0ka", "product of impact parameters of daughters to the PV", kTH1D, {{100, -0.01, 0.01}});
  }

  /// Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    /// Trigger selection
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    /// Reject collisions with negative covariance matrix elemts on the digonal
    if (collision.covXX() < 0. || collision.covYY() < 0. || collision.covZZ() < 0.) {
      return false;
    }
    return true;
  }

  /// Process function for data
  void processData(CollisionTableData::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
      /// Set magnetic field for KF vertexing
      KFParticle::SetField(magneticField);
    }
    /// Apply event selection
    if (!isSelectedCollision(collision)) {
      return;
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = CreateKFPVertexFromCollision(collision);
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

      /// Apply single track cuts as a prefilter on the daughter tracks.
      /// Transverse momentum range
      if (track1p < d_pTMin || track2p < d_pTMin) {
        continue;
      }
      /// Eta range
      if (abs(track1.eta()) > d_etaRange || abs(track2.eta()) > d_etaRange) {
        continue;
      }
      /// DCA XY of the daughter tracks to the primaty vertex
      if (track1.dcaXY() > d_dcaXYTrackPV || track2.dcaXY() > d_dcaXYTrackPV) {
        continue;
      }
      /// DCA Z of the daughter tracks to the primaty vertex
      if (track1.dcaZ() > d_dcaZTrackPV || track2.dcaZ() > d_dcaZTrackPV) {
        continue;
      }
      /// reject if the tracks have the same sign
      if (track1.sign() == track2.sign()) {
        continue;
      }
      

      KFPTrack kfpTrackPi;
      KFPTrack kfpTrackKa;

      bool CandD0 = false;
      bool CandD0bar = false;
      float nsigmaPi = 0., nsigmaKa = 0.;
      /// At the moment pT independent TPC selection. Add a minimum and Maximum momentum
      /// Apply TPC+TOF at higher momenta (TOF still uncalibrated for LHC22f).

      /// Select D0 and D0bar candidates
      if (nSigmaTPCMinPi <= track1.tpcNSigmaPi() && track1.tpcNSigmaPi() <= nSigmaTPCMaxPi && nSigmaTPCMinKa <= track2.tpcNSigmaKa() && track2.tpcNSigmaKa() <= nSigmaTPCMaxKa) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0 = true;
          kfpTrackPi = CreateKFPTrackFromTrack(track1);
          kfpTrackKa = CreateKFPTrackFromTrack(track2);
          nsigmaPi = track1.tpcNSigmaPi();
          nsigmaKa = track2.tpcNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0bar = true;
          kfpTrackPi = CreateKFPTrackFromTrack(track1);
          kfpTrackKa = CreateKFPTrackFromTrack(track2);
          nsigmaPi = track1.tpcNSigmaPi();
          nsigmaKa = track2.tpcNSigmaKa();
        } else {
          continue;
        }
      } 
      if (nSigmaTPCMinKa <= track1.tpcNSigmaKa() && track1.tpcNSigmaKa() <= nSigmaTPCMaxKa && nSigmaTPCMinPi <= track2.tpcNSigmaPi() && track2.tpcNSigmaPi() <= nSigmaTPCMaxPi) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0bar = true;
          kfpTrackPi = CreateKFPTrackFromTrack(track2);
          kfpTrackKa = CreateKFPTrackFromTrack(track1);
          nsigmaPi = track2.tpcNSigmaPi();
          nsigmaKa = track1.tpcNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0 = true;
          kfpTrackPi = CreateKFPTrackFromTrack(track2);
          kfpTrackKa = CreateKFPTrackFromTrack(track1);
          nsigmaPi = track2.tpcNSigmaPi();
          nsigmaKa = track1.tpcNSigmaKa();
        } else {
          continue;
        }
      } 
      if (!CandD0 && !CandD0bar) {
        continue;
      }

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
      histos.fill(HIST("TracksKFKa/dcaXY"), KFKaon.GetDistanceFromVertexXY(kfpVertex));
      histos.fill(HIST("TracksKFKa/length"), KFKaon.GetDecayLength());

      KFParticle KFDZero;
      const KFParticle* D0Daughters[2] = {&KFPion, &KFKaon};
      int NDaughters = 2;
      KFDZero.SetConstructMethod(2);
      KFDZero.Construct(D0Daughters, NDaughters);

      float X = 0., Y = 0., Z = 0., E = 0., Chi2 = 0., NDF = 0., P = 0., Pt = 0., PtPi = 0., PtKa = 0., Eta = 0., Phi = 0., mass = 0., decayLength = 0., decayLengthxy = 0., cosPA = -1., lifeTime = 0., massErr = 0., decayLengthErr = 0., normdecayLength = 0.;
      float distToPV = 0., deviationToPV = 0., distToPVxy = 0., deviationToPVxy = 0.;
      float deviationDaugtherTracks = 0., distanceDaugtherTracks = 0., distPiToSV = 0., deviationPiToSV = 0., distKaToSV = 0., deviationKaToSV = 0., distPiToPV = 0., distKaToPV = 0., d0pid0ka = 0., deviationPiToPV = 0., deviationKaToPV = 0.;
      float cosThetaStarPion = 0., cosThetaStarKaon = 0., decaylengthPi = 0., decaylengthKa = 0.;
      float chi2topo = 0.;

      float XGeo = 0., YGeo = 0., ZGeo = 0., EGeo = 0., Chi2Geo = 0., NDFGeo = 0., PGeo = 0., PtGeo = 0., EtaGeo = 0., PhiGeo = 0., massGeo = 0., decayLengthGeo = 0., decayLengthxyGeo = 0., lifeTimeGeo = 0., massErrGeo = 0., decayLengthErrGeo = 0., normdecayLengthGeo = 0.;
      float distToPVGeo = 0., deviationToPVGeo = 0., distToPVxyGeo = 0., deviationToPVxyGeo = 0.;
      float chi2geo = 0.;

      KFParticle KFDZero_PV = KFDZero;
      KFDZero_PV.SetProductionVertex(KFPV);

      KFParticle KFDZero_DecayVtx = KFDZero_PV;
      KFDZero_DecayVtx.TransportToDecayVertex();


      X = KFDZero_PV.GetX();
      Y = KFDZero_PV.GetY();
      Z = KFDZero_PV.GetZ();
      E = KFDZero_PV.GetE();
      Chi2 = KFDZero_PV.GetChi2();
      NDF = KFDZero_PV.GetNDF();
      P = KFDZero_PV.GetP();
      Pt = KFDZero_PV.GetPt();
      PtPi = KFPion.GetPt();
      PtKa = KFKaon.GetPt();
      Eta = KFDZero_PV.GetEta();
      Phi = KFDZero_PV.GetPhi();
      mass = KFDZero_PV.GetMass();
      decayLength = KFDZero_PV.GetDecayLength();
      decayLengthxy = KFDZero_PV.GetDecayLengthXY();
      cosPA = CosPointingAngleFromKF(KFDZero_DecayVtx, KFPV);
      lifeTime = KFDZero_PV.GetLifeTime();
      massErr = KFDZero_PV.GetErrMass();
      decayLengthErr = KFDZero_PV.GetErrDecayLength();
      normdecayLength = decayLength / decayLengthErr;
      distToPV = KFDZero_PV.GetDistanceFromVertex(KFPV);
      deviationToPV = KFDZero_PV.GetDeviationFromVertex(KFPV);
      distToPVxy = KFDZero_PV.GetDistanceFromVertexXY(KFPV);
      deviationToPVxy = KFDZero_PV.GetDeviationFromVertexXY(KFPV);
      chi2topo = KFDZero_PV.GetChi2()/KFDZero_PV.GetNDF();

      deviationDaugtherTracks = KFPion.GetDeviationFromParticle(KFKaon);
      distanceDaugtherTracks = KFPion.GetDistanceFromParticle(KFKaon);
      distPiToSV = KFPion.GetDistanceFromVertex(KFDZero_DecayVtx);
      deviationPiToSV = KFPion.GetDeviationFromVertex(KFDZero_DecayVtx);
      deviationPiToPV = KFPion.GetDeviationFromVertex(KFPV);
      deviationKaToPV = KFKaon.GetDeviationFromVertex(KFPV);
      distKaToSV = KFKaon.GetDistanceFromVertex(KFDZero_DecayVtx);
      deviationKaToSV = KFKaon.GetDeviationFromVertex(KFDZero_DecayVtx);
      distPiToPV = KFPion.GetDistanceFromVertexXY(KFPV);
      distKaToPV = KFKaon.GetDistanceFromVertexXY(KFPV);
      d0pid0ka = distPiToPV*distKaToPV;

      cosThetaStarPion = CosThetaStarFromKF(0, 421, 211, 321, KFPion, KFKaon, KFDZero_PV);
      cosThetaStarKaon = CosThetaStarFromKF(1, 421, 211, 321, KFPion, KFKaon, KFDZero_PV);
      decaylengthPi = KFPion.GetDecayLength();
      decaylengthKa = KFKaon.GetDecayLength();

      XGeo = KFDZero.GetX();
      YGeo = KFDZero.GetY();
      ZGeo = KFDZero.GetZ();
      EGeo = KFDZero.GetE();
      Chi2Geo = KFDZero.GetChi2();
      NDFGeo = KFDZero.GetNDF();
      PGeo = KFDZero.GetP();
      PtGeo = KFDZero.GetPt();
      EtaGeo = KFDZero.GetEta();
      PhiGeo = KFDZero.GetPhi();
      massGeo = KFDZero.GetMass();
      decayLengthGeo = KFDZero.GetDecayLength();
      decayLengthxyGeo = KFDZero.GetDecayLengthXY();
      lifeTimeGeo = KFDZero.GetLifeTime();
      massErrGeo = KFDZero.GetErrMass();
      decayLengthErrGeo = KFDZero.GetErrDecayLength();
      normdecayLengthGeo = decayLengthGeo / decayLengthErrGeo;
      distToPVGeo = KFDZero.GetDistanceFromVertex(KFPV);
      deviationToPVGeo = KFDZero.GetDeviationFromVertex(KFPV);
      distToPVxyGeo = KFDZero.GetDistanceFromVertexXY(KFPV);
      deviationToPVxyGeo = KFDZero.GetDeviationFromVertexXY(KFPV);
      chi2geo = KFDZero.GetChi2()/KFDZero.GetNDF();

      /// We need to get the position of the secondary vertex in a KFPVertex.
      /// Then we can also calculate the impact parameters with the KF Utils and the product if impact parameters.

      /// Remove daughter tracks from PV fit?
      /// Pt selection
      if (Pt < d_pTMinD0 || Pt > d_pTMaxD0) {
        continue;
      }
      /// Mass window selection
      if (mass < d_massMinD0 || mass > d_massMaxD0) {
        continue;
      }
      /// cosine pointing angle selection
      if (cosPA < d_cosPA) {
        continue;
      }
      /// decay length selection
      if (decayLength < d_decayLength) {
        continue;
      }
      /// decay length error selection
      if (normdecayLength < d_normdecayLength) {
        continue;
      }
      /// chi2 topological of DZero to PV
      if (chi2topo > d_chi2topoD0) {
        continue;
      }
      /// distance 3D daughter tracks at the secondary vertex
      if (distanceDaugtherTracks > d_dist3DSVDau) {
        continue;
      }
      /// cosine theta star of the pion from D0
      if (cosThetaStarPion > d_cosThetaStarPi) {
        continue;
      }
      /// cosine theta star of the kaon from D0
      if (cosThetaStarKaon > d_cosThetaStarKa) {
        continue;
      }
      /// distance Pi to SV
      if (distPiToSV > d_distPiToSV) {
        continue;
      }
      /// distance Ka to SV
      if (distKaToSV > d_distKaToSV) {
        continue;
      }
      /// product impact parameters of daughters to the PV
      if (d0pid0ka > d_d0pid0ka) {
        continue;
      }

      histos.fill(HIST("DZeroCandTopo/X"), X);
      histos.fill(HIST("DZeroCandTopo/Y"), Y);
      histos.fill(HIST("DZeroCandTopo/Z"), Z);
      histos.fill(HIST("DZeroCandTopo/E"), E);
      histos.fill(HIST("DZeroCandTopo/Chi2"), Chi2);
      histos.fill(HIST("DZeroCandTopo/NDF"), NDF);
      histos.fill(HIST("DZeroCandTopo/Chi2OverNDF"), chi2topo);
      histos.fill(HIST("DZeroCandTopo/p"), P);
      histos.fill(HIST("DZeroCandTopo/pt"), Pt);
      histos.fill(HIST("DZeroCandTopo/eta"), Eta);
      histos.fill(HIST("DZeroCandTopo/phi"), Phi);
      histos.fill(HIST("DZeroCandTopo/mass"), mass);
      histos.fill(HIST("DZeroCandTopo/massvspt"), Pt, mass);
      histos.fill(HIST("DZeroCandTopo/decayLength"), decayLength);
      histos.fill(HIST("DZeroCandTopo/decayLengthXY"), decayLengthxy);
      histos.fill(HIST("DZeroCandTopo/cosPA"), cosPA);
      histos.fill(HIST("DZeroCandTopo/lifetime"), lifeTime);
      histos.fill(HIST("DZeroCandTopo/massErr"), massErr);
      histos.fill(HIST("DZeroCandTopo/decayLengthErr"), decayLengthErr);
      histos.fill(HIST("DZeroCandTopo/distToPV"), distToPV);
      histos.fill(HIST("DZeroCandTopo/deviationToPV"), deviationToPV);
      histos.fill(HIST("DZeroCandTopo/distToPVXY"), distToPVxy);
      histos.fill(HIST("DZeroCandTopo/deviationToPVXY"), deviationToPVxy);
      histos.fill(HIST("DZeroCandTopo/deviationDaugtherTracks"), deviationDaugtherTracks);
      histos.fill(HIST("DZeroCandTopo/distanceDaugtherTracks"), distanceDaugtherTracks);
      histos.fill(HIST("DZeroCandTopo/cosThetaStarPion"), cosThetaStarPion);
      histos.fill(HIST("DZeroCandTopo/cosThetaStarKaon"), cosThetaStarKaon);
      histos.fill(HIST("DZeroCandTopo/deviationPiToSV"), deviationPiToSV);
      histos.fill(HIST("DZeroCandTopo/distPiToSV"), distPiToSV);
      histos.fill(HIST("DZeroCandTopo/deviationKaToSV"), deviationKaToSV);
      histos.fill(HIST("DZeroCandTopo/distKaToSV"), distKaToSV);
      histos.fill(HIST("DZeroCandTopo/d0pid0ka"), d0pid0ka);


      histos.fill(HIST("DZeroCandGeo/X"), XGeo);
      histos.fill(HIST("DZeroCandGeo/Y"), YGeo);
      histos.fill(HIST("DZeroCandGeo/Z"), ZGeo);
      histos.fill(HIST("DZeroCandGeo/E"), EGeo);
      histos.fill(HIST("DZeroCandGeo/Chi2"), Chi2Geo);
      histos.fill(HIST("DZeroCandGeo/NDF"), NDFGeo);
      histos.fill(HIST("DZeroCandGeo/Chi2OverNDF"), chi2geo);
      histos.fill(HIST("DZeroCandGeo/p"), PGeo);
      histos.fill(HIST("DZeroCandGeo/pt"), PtGeo);
      histos.fill(HIST("DZeroCandGeo/eta"), EtaGeo);
      histos.fill(HIST("DZeroCandGeo/phi"), PhiGeo);
      histos.fill(HIST("DZeroCandGeo/mass"), massGeo);
      histos.fill(HIST("DZeroCandGeo/massvspt"), PtGeo, massGeo);
      histos.fill(HIST("DZeroCandGeo/decayLength"), decayLengthGeo);
      histos.fill(HIST("DZeroCandGeo/decayLengthXY"), decayLengthxyGeo);
      histos.fill(HIST("DZeroCandGeo/lifetime"), lifeTimeGeo);
      histos.fill(HIST("DZeroCandGeo/massErr"), massErrGeo);
      histos.fill(HIST("DZeroCandGeo/decayLengthErr"), decayLengthErrGeo);
      histos.fill(HIST("DZeroCandGeo/distToPV"), distToPVGeo);
      histos.fill(HIST("DZeroCandGeo/deviationToPV"), deviationToPVGeo);
      histos.fill(HIST("DZeroCandGeo/distToPVXY"), distToPVxyGeo);
      histos.fill(HIST("DZeroCandGeo/deviationToPVXY"), deviationToPVxyGeo);

      if(writeTree) {
        /// Filling the D0 tree
        rowDZeroTree(PtPi,
        PtKa,
        nsigmaPi,
        nsigmaKa,
        decaylengthPi,
        decaylengthKa,
        Pt,
        mass,
        decayLength,
        decayLengthxy,
        cosPA,
        lifeTime,
        normdecayLength,
        distToPV,
        distToPVxy,
        distanceDaugtherTracks,
        distPiToSV,
        distKaToSV,
        distPiToPV,
        distKaToPV,
        d0pid0ka,
        cosThetaStarPion,
        cosThetaStarKaon,
        chi2geo,
        chi2topo,
        deviationToPV,
        deviationToPVxy,
        deviationDaugtherTracks,
        deviationPiToSV,
        deviationKaToSV,
        deviationPiToPV,
        deviationKaToPV);
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
