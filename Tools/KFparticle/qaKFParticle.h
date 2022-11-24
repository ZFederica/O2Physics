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

/// \file qaKFParticle.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace kfparticle
{
DECLARE_SOA_COLUMN(PtPi, ptpi, float);
DECLARE_SOA_COLUMN(PtKa, ptka, float);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcnsigmapi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcnsigmaka, float);
DECLARE_SOA_COLUMN(DecayLengthPi, decaylengthpi, float);
DECLARE_SOA_COLUMN(DecayLengthKa, decaylengthka, float);
DECLARE_SOA_COLUMN(PtD, ptd, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(DecayLengthD, decaylengthd, float);
DECLARE_SOA_COLUMN(DecayLengthDXY, decaylengthdxy, float);
DECLARE_SOA_COLUMN(CosPa, cospa, float);
DECLARE_SOA_COLUMN(Lifetime, lifetime, float);
DECLARE_SOA_COLUMN(NormDecayLength, normdecaylength, float);
DECLARE_SOA_COLUMN(DistDPV, distdpv, float);
DECLARE_SOA_COLUMN(DistDPVXY, distdpvxy, float);
DECLARE_SOA_COLUMN(DistDau, distdau, float);
DECLARE_SOA_COLUMN(DistPiSV, distpisv, float);
DECLARE_SOA_COLUMN(DistKaSV, distkasv, float);
DECLARE_SOA_COLUMN(DistPiPV, distpipv, float);
DECLARE_SOA_COLUMN(DistKaPV, distkapv, float);
DECLARE_SOA_COLUMN(D0PiD0Ka, d0pid0ka, float);
DECLARE_SOA_COLUMN(CosThetaStartPi, costhetastarpi, float);
DECLARE_SOA_COLUMN(CosThetaStartKa, costhetastarka, float);
DECLARE_SOA_COLUMN(Chi2Geo, chi2geo, float);
DECLARE_SOA_COLUMN(Chi2Topo, chi2topo, float);
DECLARE_SOA_COLUMN(DeviationDToPV, deviationdtopv, float);
DECLARE_SOA_COLUMN(DeviationDToPVXY, deviationdtopvxy, float);
DECLARE_SOA_COLUMN(DeviationDau, deviationdau, float);
DECLARE_SOA_COLUMN(DeviationPiSV, deviationpisv, float);
DECLARE_SOA_COLUMN(DeviationKaSV, deviationkasv, float);
DECLARE_SOA_COLUMN(DeviationPiPV, deviationpipv, float);
DECLARE_SOA_COLUMN(DeviationKaPV, deviationkapv, float);





} // namespace kfparticle
DECLARE_SOA_TABLE(TreeDZeroKF, "AOD", "TREEDZEROKF",
                  kfparticle::PtPi,
                  kfparticle::PtKa,
                  kfparticle::TPCNSigmaPi,
                  kfparticle::TPCNSigmaKa,
                  kfparticle::DecayLengthPi,
                  kfparticle::DecayLengthKa,
                  kfparticle::PtD,
                  kfparticle::Mass,
                  kfparticle::DecayLengthD,
                  kfparticle::DecayLengthDXY,
                  kfparticle::CosPa,
                  kfparticle::Lifetime,
                  kfparticle::NormDecayLength,
                  kfparticle::DistDPV,
                  kfparticle::DistDPVXY,
                  kfparticle::DistDau,
                  kfparticle::DistPiSV,
                  kfparticle::DistKaSV,
                  kfparticle::DistPiPV,
                  kfparticle::DistKaPV,
                  kfparticle::D0PiD0Ka,
                  kfparticle::CosThetaStartPi,
                  kfparticle::CosThetaStartKa,
                  kfparticle::Chi2Geo,
                  kfparticle::Chi2Topo,
                  kfparticle::DeviationDToPV,
                  kfparticle::DeviationDToPVXY,
                  kfparticle::DeviationDau,
                  kfparticle::DeviationPiSV,
                  kfparticle::DeviationKaSV,
                  kfparticle::DeviationPiPV,
                  kfparticle::DeviationKaPV);

} // namespace o2::aod