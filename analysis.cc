//*******************************************************************************/
/** development:
    Andrea Massironi 
    Olivier Boundu
    Alexandra Oliveira                                                       **/
/*******************************************************************************/
/* to run
 cd ..
 root -l 
    > gSystem->Load("libDelphes")
 
 make analysis.exe
 ./analysis.exe 
******************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
// ROOT headers
//#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <vector>
// Delphes headers
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
// fastjet headers
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Pruner.hh"
// fastjet contrib headers
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessDefinition.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
// ---------
#include "rootHistos.h" //declare the histos
#include "cuts.h" // contain the analysis cuts
#include "analysis.h" // Three 
// Verbosity
#define DEBUG 0
//#include <boost/program_options.hpp> // to option in prompt, ble!
using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;
typedef fastjet::JetDefinition::DefaultRecombiner DefRecomb;
struct myclassMin {bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { return (i.first < j.first);}} myObjMin;
struct myclassMax { bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { return (i.first > j.first);}} myObjMax;
///////////////////////////////////////////////////
// functions
int dobranches(TTree* outtree);
bool findleptons(TClonesArray *branchMissingET ,TClonesArray *branchElectron, TClonesArray *branchMuon,
		ExRootTreeReader* treeReader, 
		bool doHwwselection, TLorentzVector & l1, TLorentzVector & l2); 
bool findjets(TClonesArray *branchJet,ExRootTreeReader* treeReader);
bool myJetCollection(TClonesArray *branchEFlowTrack, TClonesArray *branchEFlowTower, TClonesArray *branchEFlowMuon , TClonesArray *branchJet );
bool isThisJetALepton(TLorentzVector* jet, TLorentzVector* l1, TLorentzVector* l2);
//////////////////////////////////////////////////////////
// declare and save histos
int decla(int);
int save_hist(int,int,int,bool);
//////////////////////////////////////////////////////////
// description of functions
// 
// isThisJetALepton ==> check the jet is not an isolated lepton
// myJetCollection ==> take all jets constituents recluster and analyse and fill the histos == if we need substructure variables,  == if no PU
// findjets ==> does the analysis cuts and fill the histos == with normal delphes jets
/////////////////////////////////////////////////////
//namespace po = boost::program_options; // to option in prompt, ble!
////////////////////////////////////////////////////
int main(){
 
 std::string inputfile = "/Users/Xanda/Documents/codes/VVcombo/lhe/pp_RS_WW_jj1.5tev30k.root"; // to loop, to make automatic
 decla(0);   
 
 TChain *chain = new TChain("Delphes");
 chain->Add(inputfile.c_str());      
 ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 //---- objects in Delphes format 
 TClonesArray *branchJet = treeReader->UseBranch("Jet");
 TClonesArray *branchElectron = treeReader->UseBranch("Electron");
 TClonesArray *branchMuon = treeReader->UseBranch("Muon");
 TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
 TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
 //---- to access constituents
 TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
 TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
 TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
 //---- events
 Long64_t allEntries = treeReader->GetEntries();
 std::cout << "** Chain contains " << allEntries << " events" << std::endl;
 // Loop over all events
 std::cout<<"start analysis"<<std::endl;
 for(entry = 0; entry < allEntries; entry++) {
   treeReader->ReadEntry(entry);  // Load selected branches with data from specified event
   TLorentzVector l1, l2; // save leptons to compare to jets -- double counted as jets
   //TLorentzVector pho1, pho2; // save to compare to jets
   //std::vector<TLorentzVector> Jets; // all the jets 
   bool findjetSub = myJetCollection(branchEFlowTrack, branchEFlowTower, branchEFlowMuon , branchJet); 
 } // close loop entry
 save_hist(1,1,1,1);
 return 0;
} // close main
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool myJetCollection(TClonesArray *branchEFlowTrack, TClonesArray *branchEFlowTower, TClonesArray *branchEFlowMuon , TClonesArray *branchJet ){
    TObject *object;
    vector<PseudoJet> jets_CA;
    Jet *jetAll; // P4 returns a TLorentzVector
    TLorentzVector momentum;
    int some =0;
    vector<fastjet::PseudoJet> particles;
    for(i = 0; i < branchJet->GetEntriesFast(); i++) {
       jetAll = (Jet*) branchJet->At(i);  int constitu =0;
       for(int j = 0; j < jetAll->Constituents.GetEntriesFast(); ++j){
          //jet = jetentry;
          momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
          //particles[j].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
          object = jetAll->Constituents.At(j);
          if(object == 0) {continue;} 
          //        if(object->IsA() == GenParticle::Class()) 
          //		{momentum += ((GenParticle*) object)->P4(); }
          else if(object->IsA() == Track::Class()){momentum += ((Track*) object)->P4(); } 
          else if(object->IsA() == Tower::Class()) {momentum += ((Tower*) object)->P4(); }
          else if(object->IsA() == Muon::Class()) { momentum += ((Muon*) object)->P4(); }
          else if(object->IsA() == Electron::Class()){ momentum += ((Electron*) object)->P4(); }
          else if(object->IsA() == Photon::Class()) { momentum += ((Photon*) object)->P4(); }
          //else if(object->IsA() == EFlowTrack::Class()) { momentum += ((EFlowTrack*) object)->P4(); }
          //particles.push_back(jet->Constituents.At(j)); // not work
           particles.push_back(momentum); some++; constitu++;
          //std::cout<<"the hardest core! "<<momentum.Pt()<<std::endl;
      } // close for jet constituents
      // check that I take all the constituents with no double counting
      //cout<<" test taking all constituents "<<constitu<<" "<< jetAll->Constituents.GetEntriesFast()<<endl;
    } // close for jets
    ///////////////////////////////////////////////////    
    JetDefinition  CA(cambridge_algorithm, Rsb); //(antikt_algorithm, RR);
    ClusterSequence cs_ca(particles, CA);
    Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax); // put baseline here
    jets_CA = sorted_by_pt(jet_selector(cs_ca.inclusive_jets())); // first we do akt jets from particles
    // just check the number of jets after recluster sometimes fluctuates -- can be unclustered energy?
    /*
    JetDefinition  akt(antikt_algorithm, 0.5);
    ClusterSequence cs_akt(particles, akt);
    vector<PseudoJet> jets_akt;
    jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets())); // first we do akt jets from particles
    //
    cout<<"njets CA= "<<jets_CA.size()<<" njets akt= "<<jets_akt.size()<< " njets original "<< branchJet->GetEntriesFast()<<endl; 
    */
    bool passcuts=false;
    if (jets_CA.size() >1) if(abs(jets_CA.at(0).eta()- jets_CA.at(1).eta()) < 1.3&& (jets_CA.at(0) + jets_CA.at(1)).m() > 890) passcuts = true;
    if(passcuts){
        //////////////////////////////////////
        // CMS tagger
        //////////////////////////////////////
        Pruner pruner(cambridge_algorithm, zcut, Rcut_factor);
        PseudoJet pruned_jet; double mprun;
        pruned_jet = pruner(jets_CA.at(0)); mprun = pruned_jet.m();
        PrunMass->Fill( mprun );
        //
        OnePass_KT_Axes     axisMode1;
        double beta1 = 1.0, R0=0.8; 
        NormalizedMeasure measureSpec1(beta1,R0); //UnnormalizedMeasure measureSpec1(beta1);
        NsubjettinessRatio nSubRatio(2, 1, axisMode1, measureSpec1); //Nsubjettiness  nSub1_beta1(1,  axisMode1,measureSpec1);
        double tau21 = nSubRatio.result(jets_CA.at(0));
        // take the defautl instead
        //cout << "tau21 "<<tau21<<endl;
        //jetAll = (Jet*) branchJet->At(0);
        //Float_t tau21true = branchJet->Tau2(); ///jetAll->Tau1();
        //cout<<" tau21 "<<tau21<<" tau21true "<<tau21true<<endl;
        Nsub->Fill(tau21);
        ///////////////////////////////////////////
        //
        Njets_passing_kLooseID->Fill(jets_CA.size(),1);
        JJmass->Fill((jets_CA.at(0)+jets_CA.at(1)).m());
        DetaJJ->Fill(abs(jets_CA.at(0).eta()- jets_CA.at(1).eta()));
        J1eta->Fill(jets_CA.at(0).eta());
        J2eta->Fill(jets_CA.at(1).eta());
        J1pt->Fill(jets_CA.at(0).pt());
        J2pt->Fill(jets_CA.at(1).pt());
    }  
    /* if(some>0) {
        // define CA jet in the Delphes card
        // ATLAS TAG
        fastjet::PseudoJet teste = particles[0];
        //std::cout<<"the hardest core! "<<teste<<std::endl;
        //fastjet::Selector jet_selector = fastjet::SelectorPtMin(MINPTJET) && fastjet::SelectorAbsRapMax(rapmax);
        //fastjet::JetDefinition akt(antikt_algorithm, jetR);
        //fastjet::ClusterSequence cs_akt(particles, akt);
        //std::vector<PseudoJet> jets_akt;
        //jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));
        // jet definition for substructure
        fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, Rsb);
        // first recluster with some large CA (needed for mass-drop)
        fastjet::ClusterSequence cs_tmp(particles, CA10);
        // next get hardest jet -- it is already sorted by pt
        std::vector<fastjet::PseudoJet> ca_jet;
        ca_jet = fastjet::sorted_by_pt(cs_tmp.inclusive_jets()); // find the cores
        //if(ca_jet[0].pt()>0)std::cout<<"the hardest core! "<<ca_jet[0].pt()<<std::endl;
        // now run mass drop tagger / compare the hardest core with the rest of the jet
        fastjet::MassDropTagger md_tagger(mu, ycut); // define the cut on mass drop
        // mu: ratio in between mass of cores, symetric splitting
        fastjet::PseudoJet tagged_jet  = md_tagger(ca_jet)[0]; // save to check if survives mass drop .. different !
        //tagged_jet = md_tagger(ca_jet);
        //if(tagged_jet.m() > 110) {std::cout<<"tag!"<<std::endl; return true;} else return false;
        // Filter definition to improve mass resolution // after
        //Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
        //JetDefinition akt(antikt_algorithm, jetR);
        //ClusterSequence cs_akt(particles, akt);
        //
    } else return false;
    */
    //cout<<" total particles in the event"<< some<< endl;
    return true;
} // close my jet collection
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool findjets(TClonesArray *branchJet,ExRootTreeReader* treeReader){ 
    Jet *jet; // P4 returns a TLorentzVector
    bool passcuts=false; 
    //cout<< branchJet->GetEntriesFast()<<endl;
    vector<TLorentzVector> jetCol; 
    std::vector<int> JetsBtag; // use as bool for b-tag -- not needed now
    vector<PseudoJet> teste;
    for(i = 0; i < branchJet->GetEntriesFast(); i++) {
        // take jet and save 4-momentum readable
        jet = (Jet*) branchJet->At(i);
        double pts = jet->PT;
        double etas = jet->Eta;
        double phis=jet->Phi;
        double masse=jet->Mass;
        TLorentzVector jetP4;
        jetP4.SetPtEtaPhiM(pts,etas,phis,masse);
        // first how many pass basic selections
        if (jetP4.Pt()>jet_ptminfinal && jetP4.Eta() < etaj){ 
            jetCol.push_back(jetP4);
            teste.push_back(jetP4);
            JetsBtag.push_back(jet->BTag ); //else JetsFlavour.push_back(99);
            // sub variables still not working in versioned version
            // TLorentzVector PrunedP4[5]; // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta 
            //if (i==0) tau21true = (Jet*) branchJet->Tau2(); 
            //cout<<" prunned mass "<< jet->PrunedP4[0]>>endl;
        } // close if basic cuts
    } // close for each jet
    // composite cuts
    if (jetCol.size() >1) if(abs(jetCol.at(0).Eta()- jetCol.at(1).Eta()) < 1.3 && (jetCol.at(0) + jetCol.at(1)).M() > 890) passcuts = true;
    //if (Jets[0].Pt() < MINPTJET) std::cout << "We have a problem; countJets = " << countJets 
    //				<< "; ijet = " << ijet << " and jet4.Pt() = " << jet4.Pt() << std::endl; 
    //std::cout<<"number of jets "<<countJets<<" "<< branchJet->GetEntriesFast()<<std::endl;
    if(passcuts){
        Njets_passing_kLooseID->Fill(jetCol.size(),1);
        JJmass->Fill((jetCol.at(0)+jetCol.at(1)).M());
        DetaJJ->Fill(abs(jetCol.at(0).Eta()- jetCol.at(1).Eta()));
        J1eta->Fill(jetCol.at(0).Eta());
        J2eta->Fill(jetCol.at(1).Eta());
        J1pt->Fill(jetCol.at(0).Pt());
        J2pt->Fill(jetCol.at(1).Pt());
    }    
    return false;
} // end findjets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool findleptons( TClonesArray *branchMissingET, TClonesArray *branchElectron,TClonesArray *branchMuon,ExRootTreeReader* treeReader, bool doHwwselection , TLorentzVector & l1, TLorentzVector & l2){

  TLorentzVector gen_l1, gen_l2;
  ///---- take the two highest pt leptons in the event (m or e)
  //    maps are ordered in ascending order
  std::map <float, int> m_maxptleptons;
  
/*  // Loop over all electrons in event
  for(i = 0; i < branchElectron->GetEntriesFast(); i++) {
   electron = (Electron*) branchElectron->At(i);
   double pt = electron->PT;
   m_maxptleptons[-pt] = -(i+1);
  }
  
  // Loop over all muons in event
  for(i = 0; i < branchMuon->GetEntriesFast(); i++) {
   muon = (Muon*) branchMuon->At(i);
   double pt = muon->PT;
   m_maxptleptons[-pt] = (i+1);
  }

  if (doHwwselection && m_maxptleptons.size() > 1) {
   // kind = 0/1 if m/e
   
   std::map<float, int>::iterator it_type_m_lepton = m_maxptleptons.begin();
   int flav1 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   
   it_type_m_lepton++;
   int flav2 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   
   nlep = 0;
   for(it_type_m_lepton = m_maxptleptons.begin(); it_type_m_lepton != m_maxptleptons.end(); it_type_m_lepton++) 
		if ( -(it_type_m_lepton->first) > 10) nlep++;
  
  
   it_type_m_lepton = m_maxptleptons.begin(); 
   
   if (nlep >= 1) {
    if (it_type_m_lepton->second>0) { 
     l1     =                 ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();
     gen_l1 = ((GenParticle*) ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->Particle.GetObject())->P4();
    } else                            {
     l1     =                 ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();
     gen_l1 = ((GenParticle*) ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->Particle.GetObject())->P4();
    }
    
    if (nlep >= 2) {
     it_type_m_lepton++;
     if (it_type_m_lepton->second>0) { 
      l2     =                 ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();
      gen_l2 = ((GenParticle*) ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->Particle.GetObject())->P4();
     }
     else                            { 
      l2     =                 ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();
      gen_l2 = ((GenParticle*) ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->Particle.GetObject())->P4();
     }
    }
   }
  }  

  if (doHwwselection && m_maxptleptons.size() >= 2) {
   
   // kind = 0/1 if m/e
   
   std::map<float, int>::iterator it_type_m_lepton = m_maxptleptons.begin();
   int flav1 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   pt1 = - it_type_m_lepton->first;
   it_type_m_lepton++;
   int flav2 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
   pt2 = - it_type_m_lepton->first;
   
   nlep = 0;
   for(it_type_m_lepton = m_maxptleptons.begin(); it_type_m_lepton != m_maxptleptons.end(); it_type_m_lepton++) {
    if ( -(it_type_m_lepton->first) > 10) nlep++;
   }
   
   //                       ee/mm          e   m           m    e
   channel =             flav1*flav2+2*(flav1>flav2)+3*(flav1<flav2);
   
   // # mumu #    channel == 0
   // # mue #     channel == 3
   // # emu #     channel == 2
   // # ee #      channel == 1
   
   pt1 = l1.Pt();
   pt2 = l2.Pt();
   mll = (l1+l2).M();
   ptll = (l1+l2).Pt();
   pzll = (l1+l2).Pz();
   dphill = l1.DeltaPhi(l2);
   
   gen_pt1 = gen_l1.Pt();
   gen_pt2 = gen_l2.Pt();
   gen_mll = (gen_l1+gen_l2).M();
   gen_ptll = (gen_l1+gen_l2).Pt();
   gen_pzll = (gen_l1+gen_l2).Pz();
   gen_dphill = gen_l1.DeltaPhi(gen_l2);

  }
  else {
   pt1 = -99.;
   pt2 = -99.;
   nlep = m_maxptleptons.size();
   channel = -1.;
   mll = -99.;
   dphill = -99.;
   
   hww_pt = -99.;
   hww_etam = -99.;
   hww_etap = -99.;
   hww_phi = -99.;
  }

  //---- met ----  
  if(branchMissingET->GetEntriesFast() > 0) {
   met   = (MissingET*) branchMissingET->At(0);
   pfmet = met->MET;
  } else pfmet = -99;
*/
  return false;
} // end findlepton
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isThisJetALepton(TLorentzVector* jet, TLorentzVector* l1, TLorentzVector* l2){
 bool isLep = false;
 if (l1) if (jet->DeltaR(*l1) < lepiso) isLep = true;
 if (l2) if (jet->DeltaR(*l2) < lepiso) isLep = true;
 return isLep;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int isample,int reco,int sample, bool shower){
    const char* Mass;
    Mass = Form("Control_pp_RSy_WW_jj1.5tev30k%d.root",reco); 
    cout<<" saved "<< Form("Control_pp_RSy_WW_jj1.5tev30k%d.root",reco)<<endl;
    TFile f1(Mass, "recreate");
    f1.cd();
    Njets_passing_kLooseID->Write();
    Nlep_passing_kLooseID->Write();
    PrunMass->Write();
    Nsub->Write();
    JJmass->Write();
    DetaJJ->Write();
    J1eta->Write();
    J2eta->Write();
    J1pt->Write();
    J2pt->Write();
    f1.Close();    
    return 0;
}
/////////////////////////////////////////////////////////////////////////
int decla(int mass){
    
    const char* label="without btag im reco";
    
    Njets_passing_kLooseID = new TH1D("njets_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Njets_passing_kLooseID->GetYaxis()->SetTitle("");
    Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering"); 
    
    Nlep_passing_kLooseID = new TH1D("nlep_passing_kLooseID_ct4",  
                                     label, 
                                     13, -0.5, 12.5);
    Nlep_passing_kLooseID->GetYaxis()->SetTitle("");
    Nlep_passing_kLooseID->GetXaxis()->SetTitle("Nleptons after showering"); 
    
    PrunMass = new TH1D("PrunMass_ct4",  
                        label, 
                        100, 0, 200);
    PrunMass->GetYaxis()->SetTitle("");
    PrunMass->GetXaxis()->SetTitle("Prunned mass leading jet (GeV)"); 
    
    Nsub = new TH1D("Nsub_ct4",  
                    label, 
                    50, 0, 1);
    Nsub->GetYaxis()->SetTitle("");
    Nsub->GetXaxis()->SetTitle("Nsubjetiness");     
    
    J1pt = new TH1D("J1pt_ct4",  
                    label, 
                    200, 0, 1000);
    J1pt->GetYaxis()->SetTitle("");
    J1pt->GetXaxis()->SetTitle("Pt leading jet (GeV)"); 
    
    J2pt = new TH1D("J2pt_ct4",  
                    label, 
                    200, 0, 1000);
    J2pt->GetYaxis()->SetTitle("");
    J2pt->GetXaxis()->SetTitle("Pt subleading jet (GeV)"); 
    
    J1eta = new TH1D("J1eta_ct4",  
                     label, 
                     50, -6, 6);
    J1eta->GetYaxis()->SetTitle("");
    J1eta->GetXaxis()->SetTitle("eta leading jet"); 
    
    J2eta = new TH1D("J2eta_ct4",  
                     label, 
                     50, -6, 6);
    J2eta->GetYaxis()->SetTitle("");
    J2eta->GetXaxis()->SetTitle("eta subleading jet"); 
    
    DetaJJ = new TH1D("DetaJJ_ct4",  
                      label, 
                      20, -1, 2);
    DetaJJ->GetYaxis()->SetTitle("");
    DetaJJ->GetXaxis()->SetTitle("Deta JJ"); 
    
    JJmass = new TH1D("JJmass_ct4",  
                      label, 
                      200, 0, 2000);
    JJmass->GetYaxis()->SetTitle("");
    JJmass->GetXaxis()->SetTitle("JJmass (GeV)"); 
    ///////////////////////////////////////////////////////////////////////////////////
    return 0;
}
//////////////////////////////////////////////////
