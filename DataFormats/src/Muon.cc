#include "CATTools/DataFormats/interface/Muon.h"

using namespace cat;

/// default constructor
Muon::Muon() {}

Muon::Muon(const reco::LeafCandidate & aMuon) :
  Lepton( aMuon ),
  isGlobalMuon_(false),
  isSoftMuon_(false),
  normalizedChi2_(-1),
  ipsig_(-1),
  numberOfValidHits_(0),
  numberOfValidMuonHits_(0),
  numberOfMatchedStations_(0),
  numberOfValidPixelHits_(0),
  trackerLayersWithMeasurement_(0)
{}

/// destructor
Muon::~Muon() {
}
