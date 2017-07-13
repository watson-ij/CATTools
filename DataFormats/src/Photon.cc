#include "CATTools/DataFormats/interface/Photon.h"

using namespace cat;

/// default constructor
Photon::Photon() {
}

Photon::Photon(const reco::LeafCandidate & aPhoton) :
  Particle( aPhoton ),
  isPF_(false),
  isTight_(false),
  isMedium_(false),
  isLoose_(false),
  passMVA_(false),
  passElVeto_(false),
  hasPixSeed_(false),
  rhoIso_(0),
  chargedHadronIso_(0),
  puChargedHadronIso_(0),
  neutralHadronIso_(0),
  photonIso_(0),
  iEtaiEta_(0),
  r9_(0),
  HoverE_(0),
  scEta_(0), scPhi_(0), scRawE_(0), scPreShE_(0),
  smearedScale_(1),
  mcMatched_(false)
{}

/// destructor
Photon::~Photon() {
}



float Photon::photonID(const std::string& name) const {
  for (std::vector<pat::Photon::IdPair>::const_iterator it = photonIDs_.begin(), ed = photonIDs_.end(); it != ed; ++it) {
    if (it->first.find(name) != std::string::npos)
      return it->second;
  }
  cms::Exception ex("Key not found");
  ex << "cat::Photon: the ID " << name << " can't be found in this cat::Photon.\n";
  ex << "The available IDs are: ";
  for (std::vector<pat::Photon::IdPair>::const_iterator it = photonIDs_.begin(), ed = photonIDs_.end(); it != ed; ++it) {
    ex << "'" << it->first << "' ";
  }
  ex << ".\n";
  throw ex;
}
