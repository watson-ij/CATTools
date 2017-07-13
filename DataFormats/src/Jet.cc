#include "CATTools/DataFormats/interface/Jet.h"
#include <unordered_map>
#include <algorithm>

using namespace cat;

/// default constructor
Jet::Jet():
  fJER_(1), fJERUp_(1), fJERDown_(1) {
}

Jet::Jet(const reco::LeafCandidate & aJet) :
  Particle( aJet ),
  looseJetID_(false),
  tightJetID_(false),
  tightLepVetoJetID_(false),
  pileupJetId_(0),
  chargedEmEnergyFraction_(0),
  vtxMass_(0),
  vtxNtracks_(0),
  vtx3DVal_(0),
  vtx3DSig_(0),
  partonFlavour_(0),
  hadronFlavour_(0),
  partonPdgId_(0),
  shiftedEnDown_(1),
  shiftedEnUp_(1),					     
  fJER_(1), fJERUp_(1), fJERDown_(1),
  qgLikelihood_(-2)
{}

/// destructor
Jet::~Jet() {
}

/// get b discriminant from label name
float Jet::bDiscriminator(const std::string & aLabel) const {
  float discriminator = -1000.;
  const std::string & theLabel = ((aLabel == "" || aLabel == "default")) ? "trackCountingHighEffBJetTags" : aLabel;
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    if(pairDiscriVector_[i].first == theLabel){
      discriminator = pairDiscriVector_[i].second;
    }
  }
  return discriminator;
}

/// print all bjet Discriminators
void Jet::bDiscriminatorPrint() const {
  for(unsigned int i=0; i!=pairDiscriVector_.size(); i++){
    std::cout << pairDiscriVector_[i].first << " = " << pairDiscriVector_[i].second << std::endl;
  }
}

float Jet::smearedRes(int direction, int era) const {
  // The era-based JER is going to be removed
  if ( era == 0 ) return direction == 0 ? fJER_ : direction > 0 ? fJERUp_ : fJERDown_;

  const auto aGenJet = this->genJet();
  if ( !aGenJet ) return 1; // No JER

  const double absEta = std::abs(this->eta());
  if (absEta >=5.0 ) return 1; // No JER
    
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  std::vector<double> etaBins = {5.0};
  std::vector<double> cJERs = {1.0};
  std::vector<double> cJERsUp = {1.0};
  std::vector<double> cJERsDn = {1.0};
  if ( era == 20160 ) {
    etaBins = {0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
    cJERs   = {1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16};
    cJERsUp = {1.109+0.008, 1.138+0.013, 1.114+0.013, 1.123+0.024, 1.084+0.011, 1.082+0.035, 1.140+0.047, 1.067+0.053, 1.177+0.041, 1.364+0.039, 1.857+0.071, 1.328+0.022, 1.16+0.029};
    cJERsDn = {1.109-0.008, 1.138-0.013, 1.114-0.013, 1.123-0.024, 1.084-0.011, 1.082-0.035, 1.140-0.047, 1.067-0.053, 1.177-0.041, 1.364-0.039, 1.857-0.071, 1.328-0.022, 1.16-0.029};
  }
  else if ( era == 2016 ) {
    etaBins = {0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
    cJERs   = {1.122, 1.167, 1.168, 1.029, 1.115, 1.041, 1.167, 1.094, 1.168, 1.266, 1.595, 0.998, 1.226};
    cJERsUp = {1.122+.026, 1.167+.048, 1.168+.046, 1.029+.066, 1.115+.030, 1.041+.062, 1.167+.086, 1.094+.093, 1.168+.120, 1.266+.132, 1.595+.175, 0.998+.066, 1.226+.145};
    cJERsDn = {1.122-.026, 1.167-.048, 1.168-.046, 1.029-.066, 1.115-.030, 1.041-.062, 1.167-.086, 1.094-.093, 1.168-.120, 1.266-.132, 1.595-.175, 0.998-.066, 1.226-.145};
  }
  else if ( era == 2015 ) {
    etaBins = {0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0};
    cJERs   = {1.095, 1.120, 1.097, 1.103, 1.118, 1.100, 1.162, 1.160, 1.161, 1.209, 1.564, 1.384, 1.216};
    cJERsUp = {1.095+.018, 1.120+.028, 1.097+.017, 1.103+.033, 1.118+0.014, 1.100+0.033, 1.162+0.044, 1.160+0.048, 1.161+0.060, 1.209+0.059, 1.564+0.321, 1.384+0.033, 1.216+0.050};
    cJERsDn = {1.095-.018, 1.120-.028, 1.097-.017, 1.103-.033, 1.118-0.014, 1.100-0.033, 1.162-0.044, 1.160-0.048, 1.161-0.060, 1.209-0.059, 1.564-0.321, 1.384-0.033, 1.216-0.050};
  }
  else if ( era == 201574 ) {
    etaBins = {0.8, 1.3, 1.9, 2.5, 3.0, 3.2, 5.0};
    cJERs   = {1.061, 1.088, 1.106, 1.126, 1.343, 1.303, 1.320};
    cJERsUp = {1.084, 1.117, 1.136, 1.220, 1.466, 1.414, 1.606};
    cJERsDn = {1.038, 1.059, 1.076, 1.032, 1.220, 1.192, 1.034};
  }
  else if (era == 2012){
    // 2012 values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
    etaBins = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
    cJERs   = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
    cJERsUp = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
    cJERsDn = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
  }
  // call lower_bound to find bin location.
  const size_t bin = std::lower_bound(etaBins.begin(), etaBins.end(), absEta) - etaBins.begin();
  const double jetPt = this->pt();
  const double genJetPt = aGenJet->pt();
  const double dPt = jetPt-genJetPt;

  double cJER = 0;
  if      ( direction == 0 ) cJER = cJERs[bin];
  else if ( direction >  0 ) cJER = cJERsUp[bin];
  else  cJER = cJERsDn[bin];
  
  const double fJER = std::max(0., (genJetPt+dPt*cJER)/jetPt);
  return fJER;
}
