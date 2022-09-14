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

// includes O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
// includes O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

// includes KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct qaKFParticle {

  // general steering settings
  Configurable<bool> isRun3{"isRun3", false, "Is Run3 dataset"}; 
  Configurable<double> magneticField{"d_bz", 5., "magnetic field"}; // ToDo: Can be modified to be taken per timestamp from CCDB. For now only configurable.

  // option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  // options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

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
 
  HistogramRegistry histos;


  void init(InitContext const&)
  {
    if (!doprocessData && !doprocessMC) {
      LOGF(info, "No enabled QA, all histograms are disabled");
      return;
    }

    const AxisSpec axisVertexPosX{500, -1., 1., "X [cm]"};
    const AxisSpec axisVertexPosY{500, -1., 1., "Y [cm]"};
    const AxisSpec axisVertexPosZ{100, -20., 20., "Z [cm]"};
    const AxisSpec axisVertexNumContrib{200, 0, 200, "Number Of contributors to the PV"};
    const AxisSpec axisVertexCov{100, -0.005, 0.005};

    const AxisSpec axisParX{300, 0, 600, "#it{x} [cm]"};
    const AxisSpec axisParY{200, -0.5, 0.5, "#it{y} [cm]"};
    const AxisSpec axisParZ{200, -11., 11., "#it{z} [cm]"};
    const AxisSpec axisParPX{VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0, "#it{p}_{x} [GeV/c]"};
    const AxisSpec axisParPY{VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0, "#it{p}_{y} [GeV/c]"};
    const AxisSpec axisParPZ{VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0, "#it{p}_{z} [GeV/c]"};

    // collisions
    histos.add("Events/posX", "", kTH1D, {axisVertexPosX});
    histos.add("Events/posY", "", kTH1D, {axisVertexPosY});
    histos.add("Events/posZ", "", kTH1D, {axisVertexPosZ});
    histos.add("Events/posXY", "", kTH2D, {axisVertexPosX, axisVertexPosY});
    histos.add("Events/posXYZ", "", kTH3D, {axisVertexPosX, axisVertexPosY, axisVertexPosZ});
    histos.add("Events/nContrib", "", kTH1D, {axisVertexNumContrib});
    histos.add("Events/vertexChi2", ";#chi^{2}", kTH1D, {{100, 0, 100}});
    histos.add("Events/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("Events/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
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
    histos.add("EventsKF/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
    histos.add("EventsKF/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

    // tracks
    histos.add("Tracks/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("Tracks/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("Tracks/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("Tracks/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("Tracks/dcaZ", "distance of closest approach in #it{z};#it{dcaZ} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("Tracks/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});

    histos.add("TracksKF/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKF/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKF/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKF/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKF/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKF/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKF/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, -10., 10.}});
    histos.add("TracksKF/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKF/dcaZ", "distance of closest approach in #it{z};#it{dcaZ} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKF/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});

    // Set magnetic field for KF vertexing
    KFParticle::SetField(magneticField);
  }

  // Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    return true;
  }

  // Process function for data
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData const& collisions, soa::Filtered<TrackTableData> const& tracks)
  {
    for (auto const& collision : collisions) {
      // Apply event selection
      if (!isSelectedCollision(collision)) {
        continue;
      }

      // set KF primary vertex
      KFPVertex kfpVertex;
      kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
      kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      kfpVertex.SetChi2(collision.chi2());
      //kfpVertex.SetNDF(...);
      kfpVertex.SetNContributors(collision.numContrib());

      KFParticle KFPV(kfpVertex);

      // fill collision parameters
      histos.fill(HIST("Events/posX"), collision.GetX());
      histos.fill(HIST("Events/posY"), collision.GetY());
      histos.fill(HIST("Events/posZ"), collision.GetZ());
      histos.fill(HIST("Events/nContrib"), collision.numContrib());
      histos.fill(HIST("Events/vertexChi2"), collision.chi2());
      histos.fill(HIST("Events/covXX"), collision.covXX());
      histos.fill(HIST("Events/covXY"), collision.covXY());
      histos.fill(HIST("Events/covXZ"), collision.covXZ());
      histos.fill(HIST("Events/covYY"), collision.covYY());
      histos.fill(HIST("Events/covYZ"), collision.covYZ());
      histos.fill(HIST("Events/covZZ"), collision.covZZ());

      histos.fill(HIST("EventsKF/posX"), kfpVertex.posX());
      histos.fill(HIST("EventsKF/posY"), kfpVertex.posY());
      histos.fill(HIST("EventsKF/posZ"), kfpVertex.posZ());
      histos.fill(HIST("EventsKF/posXY"), kfpVertex.posX(), kfpVertex.posY());
      histos.fill(HIST("EventsKF/nContrib"), kfpVertex.GetNContributors());
      histos.fill(HIST("EventsKF/vertexChi2"), kfpVertex.GetChi2());
      histos.fill(HIST("EventsKF/covXX"), kfpVertex.GetCovariance(0));
      histos.fill(HIST("EventsKF/covXY"), kfpVertex.GetCovariance(1));
      histos.fill(HIST("EventsKF/covXZ"), kfpVertex.GetCovariance(2));
      histos.fill(HIST("EventsKF/covYY"), kfpVertex.GetCovariance(3));
      histos.fill(HIST("EventsKF/covYZ"), kfpVertex.GetCovariance(4));
      histos.fill(HIST("EventsKF/covZZ"), kfpVertex.GetCovariance(5));

      // iterate over filtered tracks.
      for (const auto& track : tracks) {

        array<float, 3> trkpos_par;
        array<float, 3> trkmom_par;
        array<float, 21> trk_cov;
        track.getXYZGlo(pos_prong0);
        track.getPxPyPzGlo(trkmom_par);
        track.getCovXYZPxPyPzGlo(trk_cov);
        float trkpar_KF[6] = {trkpos_par[0], trkpos_par[1], trkpos_par[2],
                                   trkmom_par[0], trkmom_par[1], trkmom_par[2]};
        float trkcov_KF[21];
        for (int i = 0; i < 21; i++) {
          trkcov_KF[i] = trk_cov[i];
        }

        KFPTrack kfpTrack;
        kfpTrack.SetParameters(trkpar_KF);
        kfpTrack.SetCovarianceMatrix(trkcov_KF);
        kfpTrack.SetCharge(track.sign());
        //kfpTrack.SetNDF(...);
        //kfpTrack.SetChi2(...);

        // fill track parameters
        histos.fill(HIST("Tracks/x"), track.x());
        histos.fill(HIST("Tracks/y"), track.y());
        histos.fill(HIST("Tracks/z"), track.z());
        histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
        histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
        histos.fill(HIST("Tracks/length"), track.length());

        histos.fill(HIST("TracksKF/x"), kfpTrack.GetX());
        histos.fill(HIST("TracksKF/y"), kfpTrack.GetY());
        histos.fill(HIST("TracksKF/z"), kfpTrack.GetZ());
        histos.fill(HIST("TracksKF/px"), kfpTrack.GetPx());
        histos.fill(HIST("TracksKF/py"), kfpTrack.GetPy());
        histos.fill(HIST("TracksKF/pz"), kfpTrack.GetPz());
        histos.fill(HIST("TracksKF/chi2perNDF"), kfpTrack.GetChi2perNDF());


        KFParticle KFParticleFromTrack(kfpTrack, track.pidForTracking());

        histos.fill(HIST("TracksKF/dcaXY"), KFParticleFromTrack.GetDistanceFromVertexXY(kfpVertex));
        // histos.fill(HIST("TracksKF/dcaZ"), KFParticleFromTrack.dcaZ());
        histos.fill(HIST("TracksKF/length"), KFParticleFromTrack.GetDecayLength());


      }


    }
  }
  PROCESS_SWITCH(qaKFParticle, processData, "process data", true);

  // Process function for MC
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
