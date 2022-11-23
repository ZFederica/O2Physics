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

/// \file KFUtilities.h
/// \brief Utilities needed for the KFParticle package
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt

#ifndef O2_ANALYSIS_KFUTILITIES_H_
#define O2_ANALYSIS_KFUTILITIES_H_

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#include "TVector3.h"

template <typename T>
KFPVertex CreateKFPVertexFromCollision(const T& collision) {
  KFPVertex kfpVertex;
  kfpVertex.SetXYZ(collision.posX(), collision.posY(), collision.posZ());
  /// CAREFUL!!!!!! Covariance matrix elements yy and xz are switched until a central fix is provided!!!!!
  kfpVertex.SetCovarianceMatrix(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ());
  kfpVertex.SetChi2(collision.chi2());
  kfpVertex.SetNDF(2 * collision.multNTracksPV() - 3);
  kfpVertex.SetNContributors(collision.multNTracksPV());

  return kfpVertex;
}

template <typename T>
KFPTrack CreateKFPTrackFromTrack(const T& track) {

o2::track::TrackParametrizationWithError trackparCov;
trackparCov = getTrackParCov(track);
array<float, 3> trkpos_par;
array<float, 3> trkmom_par;
array<float, 21> trk_cov;

trackparCov.getXYZGlo(trkpos_par);
trackparCov.getPxPyPzGlo(trkmom_par);
trackparCov.getCovXYZPxPyPzGlo(trk_cov);
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
kfpTrack.SetNDF(track.tpcNClsFound()); // Which is the correct number?
kfpTrack.SetChi2(track.tpcChi2NCl()*track.tpcNClsFound()); // Missing in the AO2Ds


return kfpTrack;

}


float CosPointingAngleFromKF(KFParticle kfp, KFParticle PV)
{
  float v[3];
  v[0] = kfp.GetX() - PV.GetX();
  v[1] = kfp.GetY() - PV.GetY();
  v[2] = kfp.GetZ() - PV.GetZ();

  float p[3];
  p[0] = kfp.GetPx();
  p[1] = kfp.GetPy();
  p[2] = kfp.GetPz();

  float ptimesv2 = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

  if (ptimesv2 <= 0)
    return 0.;
  else {
    double cos = (v[0] * p[0] + v[1] * p[1] + v[2] * p[2]) / sqrt(ptimesv2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}

float CosPointingAngleXYFromKF(KFParticle kfp, KFParticle PV)
{
  float v[3];
  v[0] = kfp.GetX() - PV.GetX();
  v[1] = kfp.GetY() - PV.GetY();

  float p[3];
  p[0] = kfp.GetPx();
  p[1] = kfp.GetPy();

  float ptimesv2 = (p[0] * p[0] + p[1] * p[1]) * (v[0] * v[0] + v[1] * v[1]);

  if (ptimesv2 <= 0)
    return 0.;
  else {
    double cos = (v[0] * p[0] + v[1] * p[1]) / sqrt(ptimesv2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}

float CosThetaStarFromKF(int ip, int pdgvtx, int pdgprong0, int pdgprong1, KFParticle kfpvtx, KFParticle kfpprong0, KFParticle kfpprong1)
{
  float massvtx = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass();
  float massp[2];

  massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
  massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

  float pStar = sqrt((massvtx * massvtx - massp[0] * massp[0] - massp[1] * massp[1]) * (massvtx * massvtx - massp[0] * massp[0] - massp[1] * massp[1]) - 4. * massp[0] * massp[0] * massp[1] * massp[1]) / (2. * massvtx);

  float e = kfpvtx.GetE();
  float beta = kfpvtx.GetP() / e;
  float gamma = e / massvtx;

  TVector3 mom;
  TVector3 momTot(kfpvtx.GetPx(), kfpvtx.GetPy(), kfpvtx.GetPz());

  if (ip == 0) {
    mom.SetXYZ(kfpprong0.GetPx(), kfpprong0.GetPy(), kfpprong0.GetPz());
  }
  if (ip == 1) {
    mom.SetXYZ(kfpprong1.GetPx(), kfpprong1.GetPy(), kfpprong1.GetPz());
  }

  float cts = ((mom.Dot(momTot) / momTot.Mag()) / gamma - beta * sqrt(pStar * pStar + massp[ip] * massp[ip])) / pStar;

  return cts;
}

/// Calculates impact parameter in the bending plane of the particle w.r.t. a point
float impParXYFromKF(KFParticle kfpParticle, KFParticle Vertex)
{
  
  TVector2 flightLineXY((Vertex.GetX() - kfpParticle.GetX()), (Vertex.GetY() - kfpParticle.GetY()));
  TVector2 mom(kfpParticle.GetPx(), kfpParticle.GetPy());
  TVector3 mom3(kfpParticle.GetPx(), kfpParticle.GetPy(), kfpParticle.GetPz());
  float pt2 = mom.X() * mom.X() + mom.Y() * mom.Y();
  float k = flightLineXY.X() * (mom.X()/pt2) + flightLineXY.Y() * (mom.Y()/pt2);
  TVector2 d =  flightLineXY - k * mom;
  float absImpPar = sqrt(d.X() * d.X() + d.Y() * d.Y());
  TVector3 flightLine((Vertex.GetX() - kfpParticle.GetX()), (Vertex.GetY() - kfpParticle.GetY()), (Vertex.GetZ() - kfpParticle.GetZ()));
  TVector3 cross = mom3.Cross(flightLine);
  return (cross.Z() > 0. ? absImpPar : -1. * absImpPar);
}

float ldlFromKF(KFParticle kfpParticle, KFParticle PV)
{
  float dx_particle = PV.GetX() - kfpParticle.GetX();
  float dy_particle = PV.GetY() - kfpParticle.GetY();
  float dz_particle = PV.GetZ() - kfpParticle.GetZ();
  float l_particle = sqrt(dx_particle * dx_particle + dy_particle * dy_particle + dz_particle * dz_particle);
  float dl_particle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dx_particle * dx_particle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dy_particle * dy_particle + (PV.GetCovariance(5) + kfpParticle.GetCovariance(5)) * dz_particle * dz_particle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dx_particle * dy_particle + (PV.GetCovariance(3) + kfpParticle.GetCovariance(3)) * dx_particle * dz_particle + (PV.GetCovariance(4) + kfpParticle.GetCovariance(4)) * dy_particle * dz_particle);
  if (fabs(l_particle) < 1.e-8f)
    l_particle = 1.e-8f;
  dl_particle = dl_particle < 0. ? 1.e8f : sqrt(dl_particle) / l_particle;
  if (dl_particle == 0.)
    return 9999.;
  return l_particle / dl_particle;
}

float ldlXYFromKF(KFParticle kfpParticle, KFParticle PV)
{
  float dx_particle = PV.GetX() - kfpParticle.GetX();
  float dy_particle = PV.GetY() - kfpParticle.GetY();
  float l_particle = sqrt(dx_particle * dx_particle + dy_particle * dy_particle);
  float dl_particle = (PV.GetCovariance(0) + kfpParticle.GetCovariance(0)) * dx_particle * dx_particle + (PV.GetCovariance(2) + kfpParticle.GetCovariance(2)) * dy_particle * dy_particle + 2 * ((PV.GetCovariance(1) + kfpParticle.GetCovariance(1)) * dx_particle * dy_particle);
  if (fabs(l_particle) < 1.e-8f)
    l_particle = 1.e-8f;
  dl_particle = dl_particle < 0. ? 1.e8f : sqrt(dl_particle) / l_particle;
  if (dl_particle == 0.)
    return 9999.;
  return l_particle / dl_particle;
}

#endif // O2_ANALYSIS_KFUTILITIES_H_
