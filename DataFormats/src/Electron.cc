#include "CATTools/DataFormats/interface/Electron.h"

using namespace cat;

/// default constructor
Electron::Electron() {
}

Electron::Electron(const reco::LeafCandidate & aElectron) :
  Lepton( aElectron ),
  smearedScale_(1),
  relIso03_(0),
  relIso04_(0),
  ipsig_(0),
  scEta_(0),
  passConversionVeto_(false),
  isGsfCtfScPixChargeConsistent_(false),
  isEB_(false),
  snuID_(0),
  isTrigMVAValid_(false)
{}

/// destructor
Electron::~Electron() {
}

float Electron::electronID(const std::string& name) const {
  for (std::vector<pat::Electron::IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
    if (it->first.find(name) != std::string::npos)
      return it->second;
  }
  cms::Exception ex("Key not found");
  ex << "cat::Electron: the ID " << name << " can't be found in this cat::Electron.\n";
  ex << "The available IDs are: ";
  for (std::vector<pat::Electron::IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
    ex << "'" << it->first << "' ";
  }
  ex << ".\n";
  throw ex;
}
