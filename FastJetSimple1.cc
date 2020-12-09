// -*- C++ -*-
//
// Package:    FastjetEx/FastJetSimple1
// Class:      FastJetSimple1
//
/**\class FastJetSimple1 FastJetSimple1.cc FastjetEx/FastJetSimple1/plugins/FastJetSimple1.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Balashangar Kailasapathy
//         Created:  Sat, 07 Nov 2020 01:13:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//======Additional header files==================
#include "FWCore/Framework/interface/EventSetup.h" 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" //to work with reco::GenParticle
#include "DataFormats/PatCandidates/interface/Jet.h" //to work with pat::Jets
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //to work with particles inside jets
#include <vector> 
#include <string> 
#include <map>
#include <iostream>
#include <fstream>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h" // to use TFileService
#include "CommonTools/UtilAlgos/interface/TFileService.h" // to use TFileService


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class FastJetSimple1 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit FastJetSimple1(const edm::ParameterSet&);
      ~FastJetSimple1();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
//      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
	edm::EDGetTokenT<std::vector<pat::Muon>> patmuonToken;
	edm::EDGetTokenT<std::vector<pat::Jet> > patjetToken;
	edm::Service<TFileService> fs;

/*
// For muons---------------------------------
	TH1F *hist_genpt; 
	TH1F *hist_m1pt; 
	TH1F *hist_m2pt; 	
	TH1F *hist_m3pt;  
	TH1F *hist_m4pt;
	TH1I *hist_muonsnumber;
*/
// For Jets------------------------------------
	TH1F *hist_jetspt;
    TH1F *hist_njets; 
    TH1F *hist_nmuons; 
    TH1F *hist_matchedmuons; 
    TH1F *hist_matchedmuonslj; 
    TH1F *hist_leadingjetspt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FastJetSimple1::FastJetSimple1(const edm::ParameterSet& iConfig)
 :

//patmuonToken(consumes<std::vector<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonpat"))),

patjetToken(consumes<std::vector<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetpat")))
{
   //now do what ever initialization is needed
   
/*
//For muons-----------------------------------------------------------
	hist_genpt = fs->make<TH1F>("Muonpt", "High pt muons", 200, 0.0, 200.0);
	hist_m1pt = fs->make<TH1F>("Muon1pt", "First muon pt", 100, 0.0, 100.0);
	hist_m2pt = fs->make<TH1F>("Muon2pt", "Second muon pt", 100, 0.0, 100.0);
	hist_m3pt = fs->make<TH1F>("Muon3pt", "Third muon pt", 100, 0.0, 100.0);
	hist_m4pt = fs->make<TH1F>("Muon4pt", "Muon pt", 100, 0.0, 100.0);
	hist_muonsnumber = fs->make<TH1I>("Muonsnumber", "number of reconstructed muons", 20, -0.5, 19.5);
*/
//For jets-------------------------------------------------------------
	hist_jetspt = fs->make<TH1F>("Jetspt", "Jets", 200, 0.0, 200.0);
    hist_njets = fs->make<TH1F>("NJets", "number of Jets", 10, -0.5, 9.5);
    hist_nmuons = fs->make<TH1F>("NMuons", "number of JMuons", 5, -0.5, 4.5);
    hist_leadingjetspt = fs->make<TH1F>("LeadingJets", "Leading Jets pt", 200, 0.0, 200.0);
    hist_matchedmuons = fs->make<TH1F>("MatchedMuons", "Matched muons", 5, -0.5, 4.5);
    hist_matchedmuonslj = fs->make<TH1F>("MatchedMuonsLJ", "Matched muons in leading jet", 5, -0.5,4.5);
}


FastJetSimple1::~FastJetSimple1()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FastJetSimple1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	/*
//For muons-----------------------------------------------------------	
	edm::Handle<std::vector<pat::Muon> > patmuon;
	iEvent.getByToken(patmuonToken, patmuon);
	*/
// For Jets-----------------------------------------------------------
	edm::Handle<std::vector<pat::Jet> > patjet;
	iEvent.getByToken(patjetToken, patjet);
	
/*
//For muons-----------------------------------------------------------
	std::vector<float> m_pt;
	std::vector<float> m_eta;
	std::vector<float> m_phi;
	int m_number = 0;
		for (std::vector<pat::Muon>::const_iterator itGen=patmuon->begin(); itGen!=patmuon->end(); ++itGen) 
		{
		float pt = itGen->pt();
		m_pt.push_back(pt);
		m_phi.push_back(itGen->phi());
		m_eta.push_back(itGen->eta());
		m_number=m_number+1;
    		}

	hist_muonsnumber->Fill(m_number);
	std::sort(m_pt.begin(), m_pt.end()); 
	std::reverse(m_pt.begin(), m_pt.end());

    for (unsigned int i=0; i<m_pt.size() && i<4; i++){
	      hist_genpt->Fill(m_pt[i]);
    }
    if (m_pt.size() >3) {hist_m1pt->Fill(m_pt[0]); hist_m2pt->Fill(m_pt[1]); hist_m3pt->Fill(m_pt[2]); hist_m4pt->Fill(m_pt[3]);}
    else if (m_pt.size()==3) {hist_m1pt->Fill(m_pt[0]); hist_m2pt->Fill(m_pt[1]); hist_m3pt->Fill(m_pt[2]); hist_m4pt->Fill(0);}
    else if (m_pt.size()==2) {hist_m1pt->Fill(m_pt[0]); hist_m2pt->Fill(m_pt[1]); hist_m3pt->Fill(0); hist_m4pt->Fill(0);}
    else if (m_pt.size()==1) {hist_m1pt->Fill(m_pt[0]); hist_m2pt->Fill(0); hist_m3pt->Fill(0); hist_m4pt->Fill(0);}

//End of Muons---------------------------------------
*/

/*
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
*/
//For Jets----------------------------------------------------------------------------------------------------------------------
   //count the jets
    int jets = 0;  // this counts the number of jets
    int max_jets = 0; //this is the jet with maximum pt 
    int m_matched = 0;
    int m_matchedlj = 0;
    int m_notmatched = 0;
    std::vector<float> j_pt; 
	std::vector<float> m_pt;
    float jpt =0;  // this is the maximum pt of the jet 

    for (std::vector<pat::Jet>::const_iterator itJets=patjet->begin(); itJets!=patjet->end(); ++itJets) {
       int nmuons =0; //to count the number of muons in the jets
       jets=jets+1; //counters for jets
       float jetspt = itJets->pt();
       hist_jetspt->Fill(jetspt);
       j_pt.push_back(jetspt);
       if (jetspt > jpt) { //to order jets with respect of pt 
         jpt=jetspt;
         max_jets=jets; //I use max_jets to select the leading jet if (jets==maxjets)
       }
       //else {continue;}
       
       //loop on components   
       std::vector daus(itJets->daughterPtrVector());
       std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); 
       for (unsigned int k =0; k < daus.size(); k++){
           const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[k]);
           //printf("   Jets %3d: constituent %3d: pt %6.4f, pdgId %+3d, eta %6.4f, phi %6.4f \n", jets,k,cand.pt(),cand.pdgId(), cand.eta(), cand.phi());
           if (fabs(cand.pdgId())==13) {
               nmuons=nmuons+1;
               //std::cout << "The muon has pt " << cand.pt() <<std::endl;
               //if (jets==max_jets) { //in the leading jet I control if there is one of the first muon 
                      for (int i=0; i<4; i++){
                           if (fabs(cand.pt()-m_pt[i])<0.01){
			     //std::cout << "-----------------Matched muon!!------------------" <<std::endl;
                              //printf("Our muon has  pt %6.4f, eta %6.4f, phi %6.4f \n", m_pt[i], m_eta[i], m_phi[i]);
                              //printf("Jet muon has  pt %6.4f, eta %6.4f, phi %6.4f \n", cand.pt(), cand.eta(), cand.phi());
                              m_matched = m_matched +1;
       			      if (jets==max_jets){m_matchedlj=m_matchedlj+1;}
                           }
    			   else {
			      //std::cout<< "-------no matched muons-----" <<std::endl;
 			      m_notmatched=m_notmatched+1;
			   }
                      }
               //}
           } 
       }
       //std::cout << "The number of muons in the " << jets << " jet is " << nmuons <<std::endl;
       hist_nmuons-> Fill(nmuons);
       if (nmuons>=2){
          //std::cout <<"There are at least two muons in this jets" <<std::endl;
       }
 
         /*for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
                const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
                printf("   Jets %3d: constituent %3d: pt %6.4f, pdgId %+3d, eta %6.4f, phi %6.4f \n", jets,i2,cand.pt(),cand.pdgId(), cand.eta(), cand.phi());
            }*/
       
    }
    //std::cout<<"The number of jets is " <<jets <<std::endl;
    hist_njets->Fill(jets); //histogram with the number of jets
    //std::cout << "The maximum pt is " <<jpt << " for jet " << max_jets << std::endl;
    //for (unsigned int k =0; k<j_pt.size(); k++)
	//{std::cout << "Non ordinati" << j_pt[k] << std::endl;}
    std::sort(j_pt.begin(), j_pt.end()); 
    std::reverse(j_pt.begin(), j_pt.end());
    for (unsigned int i=0; i<j_pt.size() && i<1; i++){
      hist_leadingjetspt->Fill(j_pt[i]);
    }
    //for (unsigned int k =0; k<j_pt.size(); k++)
	//{std::cout << j_pt[k] << std::endl;}
    //std::cout<<"The number of muons in jets is " <<m_matched+m_notmatched << "; " << m_matched << " are muons matched while " <<m_matchedlj << " are muons matched in leading jet" <<std::endl;
   hist_matchedmuons->Fill(m_matched);
   hist_matchedmuonslj->Fill(m_matchedlj);
}

//End of Jets----------------------------------------------------------------------------------------------------------------------

/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}
*/

// ------------ method called once each job just before starting event loop  ------------
void
FastJetSimple1::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
FastJetSimple1::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FastJetSimple1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetSimple1);
