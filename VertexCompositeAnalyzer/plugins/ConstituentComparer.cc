#ifndef CONSTITUENT_COMPARER
#define CONSTITUENT_COMPARER

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>

// #define DEBUG_OTTSUB

class ConstituentComparer {
    using LV=math::XYZTLorentzVector;
public:
    ConstituentComparer() : _nonMuon(false), _needDimuonIncl(false) {}
    ConstituentComparer(bool nonMuon, bool needDimuonIncl) : _nonMuon(nonMuon), _needDimuonIncl(needDimuonIncl) {}
    ~ConstituentComparer() {}

    // Check a jet for all OTT candidates and return the modified PTs
    std::vector<std::pair<LV, LV> > checkOttForEveryJet(const std::vector<pat::Jet>& jets, const reco::VertexCompositeCandidate& ott) {
        auto* newPTs = new std::vector<std::pair<LV, LV> >{};
        bool isDimuonIncl = false;
        for (const auto& jet : jets) {
            auto res = checkJetForACandidate(jet, ott, _nonMuon, isDimuonIncl);
            if( (!_needDimuonIncl) || isDimuonIncl) newPTs->push_back(std::make_pair(LV(jet.p4()),res));
        }
        return std::move(*newPTs);
    }
    std::vector<std::pair<LV, LV> > checkJetForEveryOtt(const pat::Jet& jet, const reco::VertexCompositeCandidateCollection& ottCol) {
        auto* newPTs = new std::vector<std::pair<LV, LV> >{};
        for (const auto& ott : ottCol) {
            bool isDimuonIncl = false;
            auto res = checkJetForACandidate(jet, ott, _nonMuon, isDimuonIncl);
            if( (!_needDimuonIncl) || isDimuonIncl) newPTs->push_back(std::make_pair(LV(jet.p4()),res));
        }
        return std::move(*newPTs);
    }

    // Check a jet for a single OTT candidate and return the modified PT
    LV checkJetForACandidate(const pat::Jet& jet, const reco::VertexCompositeCandidate& ott, bool nonMuon, bool& isDimuonIncl) {
        std::vector<LV> toAdd;

        // Loop over the daughters of the OTT candidate
        for (size_t i = 0; i < ott.numberOfDaughters(); ++i) {
            const auto* dau = ott.daughter(i);
            bool isIncl = false;

#ifdef DEBUG_OTTSUB
            std::cout << "Daughter : " << dau->pdgId() << ", " << dau->pt() << ", " << dau->eta() << std::endl;
#endif
            // Check if the daughter is a constituent of the jet
            for (const auto& PFcand : jet.getJetConstituentsQuick()) {
#ifdef DEBUG_OTTSUB
                std::cout << "Constituent : " << PFcand->pdgId() << ", " << PFcand->pt() << ", " << PFcand->eta() << " -> match ? " << isConstituent(*PFcand, *dau) << std::endl;;
#endif
                if (isConstituent(*PFcand, *dau)) {
                    if((*PFcand).pdgId() == 1 ) isDimuonIncl = true;
                    isIncl = true;
                    break;
                }
            }

            // If not included, add its momentum to the list
            if (!isIncl) {
                toAdd.push_back(dau->p4());
#ifdef DEBUG_OTTSUB
                std::cout << "To be added: " << dau->pt() << std::endl;
#endif                
            }
        }
#ifdef DEBUG_OTTSUB
        std::cout << "Total new adds : " << toAdd.size() << std::endl; 
#endif        

        // Adjust the jet's 4-momentum
        LV jetP = jet.p4();
        for (const auto& vec4 : toAdd) {
            jetP += vec4;
        }
        return jetP;
    }

    // Compare a jet constituent to a candidate
    bool isConstituent(const reco::Candidate& jetp, const reco::Candidate& acand) {
        if (jetp.charge() != acand.charge()) return false;
        if (std::abs(jetp.eta() - acand.eta()) > eps) return false;
        if (std::abs(jetp.phi() - acand.phi()) > eps) return false;
        if (std::abs(jetp.pt() - acand.pt()) > eps) return false;
        return true;
    }

private:
    const bool _nonMuon;     // Whether to exclude muons
    const bool _needDimuonIncl; // To check dimuon 
    const double eps = 0.03; // Tolerance for floating-point comparison
};

#endif