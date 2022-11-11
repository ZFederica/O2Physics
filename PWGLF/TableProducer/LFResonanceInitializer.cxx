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

/// \file LFResonanceInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers
struct reso2initializer {

  Produces<aod::ResoCollisions> resoCollisions;
  Produces<aod::ResoDaughters> reso2tracks;
  Produces<aod::ResoDaughtersMC> reso2mctracks;

  // Configurables
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"}; // Choose if running on converted data or pilot beam
  Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0s"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC and TOF PID cut (loose, improve performance)"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> isRun2{"isRun2", 0, "if Run2: demand TPC refit"};

  /// DCA Selections for V0
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.05, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"}; // Pre-selection
  Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"}; // Pre-selection
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  Configurable<double> cMinV0Radius{"cMinV0Radius", 5.0, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};

  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hGoodTrackIndices", "hGoodTrackIndices", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
                                             {"hGoodV0Indices", "hGoodV0Indices", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
                                           },
                               OutputObjHandlingPolicy::QAObject};

  // Pre-filters for efficient process
  // Filter tofPIDFilter = aod::track::tofExpMom < 0.f || ((aod::track::tofExpMom > 0.f) && ((nabs(aod::pidtof::tofNSigmaPi) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < pidnSigmaPreSelectionCut))); // TOF
  Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta;                                                                                                                                                    // Eta cut
  Filter trackCutFilter = requireGlobalTrackInFilter();                                                                                                                                                      // Global track cuts
  Filter collisionFilter = nabs(aod::collision::posZ) < ConfEvtZvtx;

  using ResoEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using ResoEventsMC = soa::Join<ResoEvents, aod::McCollisionLabels>;
  using ResoTracks = aod::Reso2TracksPIDExt;
  using ResoTracksMC = soa::Join<ResoTracks, aod::McTrackLabels>;
  using ResoV0s = aod::V0Datas;
  using ResoV0sMC = soa::Join<ResoV0s, aod::McV0Labels>;

  Preslice<soa::Filtered<ResoTracks>> tracksbyCollisionID = aod::track::collisionId;
  Preslice<ResoV0s> v0sbyCollisionID = aod::v0data::collisionId;

  template <bool isMC, typename CollisionType, typename TrackType>
  bool IsTrackSelected(CollisionType const& collision, TrackType const& track)
  {
    // Track selection
    qaRegistry.fill(HIST("hGoodTrackIndices"), 0.5);
    // MC case can be handled here
    // if constexpr (isMC) {
    //   // MC check
    // }
    qaRegistry.fill(HIST("hGoodTrackIndices"), 1.5);
    return true;
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool IsV0Selected(CollisionType const& collision, V0Type const& v0, TrackType const& track)
  {
    // V0 selection
    qaRegistry.fill(HIST("hGoodV0Indices"), 0.5);

    auto postrack = v0.template posTrack_as<TrackType>();
    auto negtrack = v0.template negTrack_as<TrackType>();

    if (postrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (negtrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 1.5);

    if (fabs(postrack.dcaXY()) < cMinV0PosDCArToPVcut)
      return false;
    if (fabs(negtrack.dcaXY()) < cMinV0NegDCArToPVcut)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 2.5);

    if ((v0.v0radius() > cMaxV0Radius) || (v0.v0radius() < cMinV0Radius))
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 3.5);
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinV0CosPA)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 4.5);

    // MC case can be handled here
    // if constexpr (isMC) {
    //   // MC check
    // }
    return true;
  }

  // Filter for all tracks
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // Dummy lists for V0s
    // Loop over tracks
    for (auto& track : tracks) {
      // LOGF(info, "Track ID: %d", track.globalIndex());
      if (!IsTrackSelected<isMC>(collision, track))
        continue;

      // Add PID selection criteria here
      uint8_t tpcPIDselections = 0;
      uint8_t tofPIDselections = 0;
      // TPC PID
      if (std::abs(track.tpcNSigmaPi()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kPion;
      if (std::abs(track.tpcNSigmaKa()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kKaon;
      if (std::abs(track.tpcNSigmaPr()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kProton;
      // TOF PID
      if (track.hasTOF()) {
        tofPIDselections |= aod::resodaughter::PDGtype::kHasTOF;
        if (std::abs(track.tofNSigmaPi()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kPion;
        if (std::abs(track.tofNSigmaKa()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kKaon;
        if (std::abs(track.tofNSigmaPr()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kProton;
      }

      reso2tracks(resoCollisions.lastIndex(),
                  track.pt(),
                  track.px(),
                  track.py(),
                  track.pz(),
                  track.eta(),
                  track.phi(),
                  aod::resodaughter::DaughterType::kTrack,
                  track.dcaXY(),
                  childIDs,
                  track.sign(),
                  (uint8_t)track.tpcNClsCrossedRows(),
                  track.dcaXY(),
                  track.dcaZ(),
                  track.x(),
                  track.alpha(),
                  tpcPIDselections,
                  tofPIDselections,
                  track.tpcNSigmaPi(),
                  track.tpcNSigmaKa(),
                  track.tpcNSigmaPr(),
                  track.tofNSigmaPi(),
                  track.tofNSigmaKa(),
                  track.tofNSigmaPr(),
                  0, 0, 0,
                  0, 0, 0, 0);
      if constexpr (isMC) {
        fillMCTracks(track);
      }
    }
  }

  // Filter for all V0s
  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  void fillV0s(CollisionType const& collision, V0Type const& v0s, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto& v0 : v0s) {
      // LOGF(info, "V0 ID: %d", v0.globalIndex());
      if (!IsV0Selected<isMC>(collision, v0, tracks))
        continue;
      childIDs[0] = v0.posTrackId();
      childIDs[1] = v0.negTrackId();
      reso2tracks(resoCollisions.lastIndex(),
                  v0.pt(),
                  v0.px(),
                  v0.py(),
                  v0.pz(),
                  v0.eta(),
                  v0.phi(),
                  aod::resodaughter::DaughterType::kV0,
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                  childIDs,
                  0,
                  0,
                  v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                  0,
                  v0.x(), 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  v0.dcaV0daughters(), v0.mLambda(), v0.mAntiLambda(),
                  v0.v0radius(), v0.x(), v0.y(), v0.z());
    }
  }

  template <typename TrackType>
  void fillMCTracks(TrackType const& track)
  {
    std::vector<int> mothers = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    if (track.has_mcParticle()) {
      // Tempora lambda function to prevent error in build
      auto getMothersIndeces = [&](auto const& theMcParticle) {
        std::vector<int> lMothersIndeces{};
        for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
          LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
          lMothersIndeces.push_back(lMother.globalIndex());
        }
        return lMothersIndeces;
      };
      auto getDaughtersIndeces = [&](auto const& theMcParticle) {
        std::vector<int> lDaughtersIndeces{};
        for (auto& lMother : theMcParticle.template daughters_as<aod::McParticles>()) {
          LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
          lDaughtersIndeces.push_back(lMother.globalIndex());
        }
        return lDaughtersIndeces;
      };
      //
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers())
        mothers = getMothersIndeces(particle);
      if (particle.has_daughters())
        daughters = getDaughtersIndeces(particle);
      reso2mctracks(particle.pdgCode(),
                    &mothers[0],
                    &daughters[0],
                    particle.isPhysicalPrimary(),
                    particle.producedByGenerator());
    } else
      // No MC particle associated
      reso2mctracks(0, &mothers[0], &daughters[0], 0, 0);
  }

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);
  }

  void process(soa::Filtered<ResoEvents>::iterator const& collision,
               soa::Filtered<ResoTracks> const& tracks, ResoV0s const& V0s, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
    if (ConfStoreV0) {
      fillV0s<false>(collision, V0s, tracks);
    }
  }
  PROCESS_SWITCH(reso2initializer, process, "Process for data", true);

  void processMC(soa::Filtered<soa::Join<ResoEvents, aod::McCollisionLabels>>::iterator const& collision,
                 aod::McCollisions const& mcCols, soa::Filtered<ResoTracksMC> const& tracks,
                 ResoV0sMC const& V0s, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    if (ConfStoreV0) {
      fillV0s<true>(collision, V0s, tracks);
    }
  }
  PROCESS_SWITCH(reso2initializer, processMC, "Process for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reso2initializer>(cfgc, TaskName{"lf-reso2initializer"}),
  };
}
