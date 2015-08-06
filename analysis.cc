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
// git branch VVcombo
// git checkout VVcombo
// git push origin VVcombo
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
bool myJetCollection(TClonesArray *branchEFlowTrack, TClonesArray *brancheflowPhoton, TClonesArray *EFlowNeutralHadron , TClonesArray *branchJet, TClonesArray *branchEFlowMuon,TClonesArray *branchElectron, TClonesArray *branchMuon, TClonesArray *branchMissingET, TClonesArray * branchParticle, int file , int & counterHP , int & counterLP);
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
 string path = "/Users/Xanda/Documents/programs/Delphes-3.2.0/Delphes_AN_withSub/lhe/"; // to loop, to make automatic
 string inputfile; string end = ".root";
 for (unsigned i=0 ; i< 5; i++) {
 inputfile= path + file[i] + end;
 cout<<"reading from: "<< inputfile<< endl;
 decla(0);   
    //return 0;
 TChain *chain = new TChain("Delphes");
 chain->Add(inputfile.c_str());      
 ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 //---- objects in Delphes format 
 TClonesArray *branchJet = treeReader->UseBranch("Jet"); // have nsub
 TClonesArray *branchElectron = treeReader->UseBranch("Electron");
 TClonesArray *branchMuon = treeReader->UseBranch("Muon");
 TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
 TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
 //---- to access constituents -- not from the jet
 TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
 TClonesArray *brancheflowPhoton = treeReader->UseBranch("EFlowPhoton");
 TClonesArray *EFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
     //TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
 TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
 // to PU
//     add Branch Calorimeter/eflowTracks EFlowTrack Track
//     add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
//     add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower
     
 //---- events
 Long64_t allEntries = treeReader->GetEntries();
 std::cout << "** Chain contains " << allEntries << " events" << std::endl;
 // Loop over all events
 std::cout<<"start analysis; selection = "<<std::endl;
 int counterHP=0, counterLP=0, counterall=0;
 //int sel = selections[i];
 for(entry = 0; entry < allEntries; entry++) { //allEntries
   treeReader->ReadEntry(entry);  // Load selected branches with data from specified event
   TLorentzVector l1, l2; // save leptons to compare to jets -- double counted as jets
   //TLorentzVector pho1, pho2; // save to compare to jets
   //std::vector<TLorentzVector> Jets; // all the jets 
   bool findjetSub = myJetCollection(branchEFlowTrack, brancheflowPhoton, EFlowNeutralHadron , branchJet, branchEFlowMuon,branchElectron, branchMuon, branchMissingET, branchParticle, i, counterHP, counterLP); 
     counterall++;
     //findjetSub = findjets(branchJet, treeReader);
 } // close loop entry
 save_hist(i,1,1,1);
 cout<<counterHP<<" "<<counterLP<<" "<<counterall<<endl;
 } // close for file
 return 0;
} // close main
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool myJetCollection(TClonesArray *branchEFlowTrack, TClonesArray *brancheflowPhoton, TClonesArray *EFlowNeutralHadron , TClonesArray *branchJet, TClonesArray *branchEFlowMuon ,TClonesArray *branchElectron, TClonesArray *branchMuon , TClonesArray *branchMissingET, TClonesArray * branchParticle , int file , int & counterHP , int & counterLP){
    ////////////////////////////////////////////////////
    // pure constituents instead
    // separate the isolated leptons!
    TLorentzVector leptonT; double pfmet, pxmet, pymet; PseudoJet Lepton, Lepton2; double metcut=0; 
    vector<PseudoJet> Muons, Electrons; PseudoJet gen_met_vector; // to do distances
    for(i = 0; i < branchElectron->GetEntriesFast(); i++) {
        electron = (Electron*) branchElectron->At(i);
        leptonT = electron->P4();
        Electrons.push_back(fastjet::PseudoJet(leptonT.Px(),leptonT.Py(),leptonT.Pz(),leptonT.E())); 
    }
    for(i = 0; i < branchMuon->GetEntriesFast(); i++) {
        muon = (Muon*) branchMuon->At(i);
        leptonT = muon->P4();
        Muons.push_back(fastjet::PseudoJet(leptonT.Px(),leptonT.Py(),leptonT.Pz(),leptonT.E())); 
    }
    int nmu=Muons.size(); int nel= Electrons.size(); int nlep= nmu +nel;
    Nlep_passing_kLooseID->Fill(nlep,1);
    // take MET
    if(nlep==1 && nmu>0 && selections[file] == 1 ) Lepton = Muons.at(0);
    else if(nlep==1 && nel>0 && selections[file] == 1) Lepton = Electrons.at(0);
    else if(nlep==2 && nmu>1 && selections[file] == 2) {Lepton = Muons.at(0); Lepton2 = Muons.at(0);}
    else if(nlep==2 && nel>1 && selections[file] == 2) {Lepton = Electrons.at(0); Lepton2 = Electrons.at(1);}
    if(branchMissingET->GetEntriesFast() > 0 && nlep==1) {
        met   = (MissingET*) branchMissingET->At(0);
        pfmet = met->MET;
        double etamet = met->Eta;
        double phimet = met->Phi;

        pxmet = pfmet*(TMath::Cos(phimet));
        pymet = pfmet*(TMath::Sin(phimet));
        //gen_met_vector = fastjet::PseudoJet(leptonT.Px(),leptonT.Py(),leptonT.Pz(),leptonT.E()); // colinear
    } else pfmet =-100;
    // take the neutrino
    int countNeu=1;
    /*if(Muons.size() ==0 && Electrons.at(0).pt() > ptE && Electrons.at(0).eta() < etaE) 
    {metcut= METenuJ; Lepton = Electrons.at(0);}
    else if(Electrons.size() ==0 && Muons.at(0).pt() > ptMu && Muons.at(0).eta() < etaMu ) 
    {metcut = METmunuJ;  Lepton = Muons.at(0);} 
    else metcut = 10000;// = lepton veto
     */
    /*    for(int iPart = 0; iPart < branchParticle->GetEntriesFast(); iPart++){
        particle = (GenParticle*) branchParticle->At(i);
        int pdgCode = TMath::Abs(particle->PID);
        int IsPU = particle->IsPU;
        int status = particle->Status;
        if(IsPU == 0 && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) )cout<<"all "<<branchParticle->GetEntriesFast() <<" status "<< status<< " pdgid "<< pdgCode<<endl;
        if (IsPU == 0 && status == 2 && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) ) {
            gen_met_vectorT = gen_met_vectorT + particle->P4(); countNeu++;
        }
    }
    if(countNeu>0) gen_met_vector = fastjet::PseudoJet(gen_met_vectorT.Px(),gen_met_vectorT.Py(),gen_met_vectorT.Pz(),gen_met_vectorT.E());
 */
    PseudoJet VV;
    //////////////////////////////////////////////////////
    vector<PseudoJet> particles_pure;
    TLorentzVector particles_pureT;
    //cout<<" all tracks "<< branchEFlowTrack->GetEntriesFast() << endl; // check with isol leptons
    //Track *track; // P4 returns a TLorentzVector
    for(i = 0; i < branchEFlowTrack->GetEntriesFast(); i++) {
        // implememnt check iso lepton
        track = (Track*) branchEFlowTrack->At(i);
        particles_pureT = track->P4(); // to undesrtand
        PseudoJet trackTAKE = fastjet::PseudoJet(particles_pureT.Px(),particles_pureT.Py(),particles_pureT.Pz(),particles_pureT.E());
        //if(Lepton.pt()>0)if(trackTAKE.delta_R(Lepton) > 0.8) 
            particles_pure.push_back(trackTAKE);
    }
    //cout<< branchEFlowTrack->GetEntriesFast()<<" "<<endl;
    //Tower *tower; // P4 returns a TLorentzVector
    for(i = 0; i < brancheflowPhoton->GetEntriesFast(); i++) {
        tower = (Tower*) brancheflowPhoton->At(i);
         particles_pureT = tower->P4(); // to undesrtand
         if(particles_pureT.Pt() > const_ptmin) particles_pure.push_back(fastjet::PseudoJet(particles_pureT.Px(),particles_pureT.Py(),particles_pureT.Pz(),particles_pureT.E()));
    }
    for(i = 0; i < EFlowNeutralHadron->GetEntriesFast(); i++) {
        // implememnt check iso lepton
        tower = (Tower*) EFlowNeutralHadron->At(i);
        particles_pureT = tower->P4(); // to undesrtand
        if(particles_pureT.Pt() > const_ptmin) particles_pure.push_back(fastjet::PseudoJet(particles_pureT.Px(),particles_pureT.Py(),particles_pureT.Pz(),particles_pureT.E()));
    }    
    /*for(i = 0; i < branchEFlowMuon->GetEntriesFast(); i++) {
        // implememnt check iso lepton
        tower = (Tower*) branchEFlowMuon->At(i);
        particles_pureT = tower->P4(); // to undesrtand
        if(particles_pureT.Pt() > const_ptmin) 
            particles_pure.push_back(fastjet::PseudoJet(particles_pureT.Px(),particles_pureT.Py(),particles_pureT.Pz(),particles_pureT.E()));
    }*/
    // do the jet
    vector<PseudoJet> jets_CA_track;
    JetDefinition  CA(cambridge_algorithm, R0); 
    ClusterSequence cs_ca_track(particles_pure, CA);
    Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax); // put baseline here
    jets_CA_track = sorted_by_pt(jet_selector(cs_ca_track.inclusive_jets())); // first we do akt jets from particles
    // is this jet a lepton?
    vector<PseudoJet> jets_CA_final; //cout<<" wrong! 3 "<<Muons.size()<<" "<< Electrons.size()<<endl;
    int testnlep=0;
    for (unsigned j = 0; j<jets_CA_track.size(); j++ ) {    
        bool notLep = true;
        //if(Muons.size()>0)for (unsigned i = 0; i<Muons.size(); i++ ) if (jets_CA_track.at(j).delta_R(Muons.at(i)) < lepiso) notLep = false;
        //if(Electrons.size()>0)for (unsigned l = 0; l<Electrons.size(); l++ ) if (jets_CA_track.at(j).delta_R(Electrons.at(l)) < lepiso) notLep = false;
        if (jets_CA_track.at(j).delta_R(Lepton) < lepiso) notLep = false;
        if(selections[file] == 2) if (jets_CA_track.at(j).delta_R(Lepton2) < lepiso) notLep = false;
        if(notLep) jets_CA_final.push_back(jets_CA_track.at(j)); else testnlep++;
    }
    //cout<<"n isolated lepton removed "<<testnlep<< " njets "<< jets_CA_track.size() <<" njets cleaned "<< jets_CA_final.size()<<endl;
    //cout<<" all jets "<< jets_CA_track.size()<<" non lep = "<<jets_CA_final.size()<<endl;
    Njets_passing_kLooseID->Fill(jets_CA_final.size(),1);
    bool passcuts=false;
    double mprunmin,mprunmax, ptv;
    //if(selections[file] == 0 ) cout<<" had!"<<endl;  else   if(selections[file] == 1 )  cout<<" lep!"<<endl; 
    //int tau21LP=0,tau21HP=0;
    if(selections[file] == 0 && jets_CA_final.size() >1) { 
        mprunmin = mprunjj_min; mprunmax = mprunjj_max;
        //mprunmin = 0; mprunmax = 1000;
        if(abs(jets_CA_final.at(0).eta()- jets_CA_final.at(1).eta()) < Deltay && (jets_CA_final.at(0) + jets_CA_final.at(1)).m() > Mvvjj) passcuts = true;
    } else if(selections[file] == 1 && jets_CA_final.size() >0 && nlep ==1 ) { 
         mprunmin = mprunjlnu_min; mprunmax = mprunjlnu_max;
        //mprunmin = 0; mprunmax = 1000;
        ////////////////////////////////////////////////////////
         //// reco the neutrino
        /////////////////////////////////////////////////////////////
        double wt = (Lepton.px()*pxmet) + (Lepton.py()*pymet);
        double mu = (pow(wmass,2)/2 + wt);
        double aw = (pow(Lepton.pz(),2)-pow(Lepton.e(),2));
        double bw = 2*mu*Lepton.pz();
        double cw = pow(mu,2)-pow(Lepton.e(),2)*pow(pfmet,2);
        double discriminant = pow(bw,2) - 4*aw*cw ; 
        double pnuz;//,recowm,recowpt,recoweta,recowphi;  
        PseudoJet lepW; 
        if(discriminant>=0){
            //recotruth=1; 
            double pznu1 = (-bw - sqrt(discriminant))/(2*aw);
            double pznu2 = (-bw + sqrt(discriminant))/(2*aw);
            //cout<<"pz1 "<<pznu1 << " pz2 " <<pznu2 <<" pz "<<neutrinos.at(0).pz()<<" pzl "<<leptons.at(0).pz()<<endl;
            //cout<<"dumb "<< (neutrinos.at(0)+leptons.at(0)).m() << " calculated " <<
            //sqrt(pow(leptons.at(0).e()+neutrinos.at(0).e(),2)-wt-pow(leptons.at(0).pz()+neutrinos.at(0).pz(),2))
            //  <<endl;
            pnuz=TMath::Min(pznu1,pznu2); 
            double enu = sqrt(pow(pxmet,2)+pow(pymet,2)+pow(pnuz,2));
            double pxnu=pxmet,pynu=pymet;
            lepW = Lepton+ fastjet::PseudoJet(pxnu,pynu,pnuz,enu);
            VV = lepW + jets_CA_final.at(0);
            //cout << lepW.m()<< " "<< VV.m()<<endl;
            if( jets_CA_final.at(0).pt() > ptVlnu && pfmet > metcut  && jets_CA_final.at(0).eta() < etaJ) passcuts=true; // && VV.m() > Mvvlnuj
      } // close reco the neutrino
        
      //////////////////////////////////////////////////////////////
    } else if(selections[file] == 2 && jets_CA_final.size() >0 && nlep ==2 ){
            VV = Lepton2 + Lepton + jets_CA_final.at(0);
        mprunmin = mprunjll_min; mprunmax = mprunjll_max;
        //cout<<"here"<<endl;
        if( jets_CA_final.at(0).pt() > ptVlnu   && jets_CA_final.at(0).eta() < etaJ) passcuts=true;
    }// close selections
    if(passcuts){
             //////////////////////////////////////
             // CMS tagger
             //////////////////////////////////////
             Pruner pruner(cambridge_algorithm, zcut, Rcut_factor);
             PseudoJet pruned_jet; double mprun;
             pruned_jet = pruner(jets_CA_final.at(0)); mprun = pruned_jet.m();
             PrunMass->Fill( mprun );
             if(selections[file] == 0 ) PrunMassSub->Fill(pruner(jets_CA_final.at(1)).m());
             if(selections[file] == 1 ) DRJl->Fill(jets_CA_final.at(0).delta_R(Lepton));
             //
        double tau21; 
        if( mprun >mprunmin && mprun < mprunmax){
             OnePass_KT_Axes     axisMode1;
             //NormalizedMeasure measureSpec1(beta1,R0); //
             UnnormalizedMeasure measureSpec1(beta1);
             NsubjettinessRatio nSubRatio(2, 1, axisMode1, measureSpec1); //Nsubjettiness  nSub1_beta1(1,  axisMode1,measureSpec1);
             tau21 = nSubRatio.result(jets_CA_final.at(0));
             // take the defautl instead
             //cout << "tau21 "<<tau21<<endl;
             //jetAll = (Jet*) branchJet->At(0);
             //Float_t tau21true = branchJet->Tau2(); ///jetAll->Tau1();
             //cout<<" tau21 "<<tau21<<endl;
             Nsub->Fill(tau21);
        } else tau21 =-10; // close mprun
        passcuts=false; 
        if(tau21 <= 0.5 && tau21>0) {
            passcuts=true; counterHP++;
            if(selections[file] == 0 ) JJmassHP->Fill((jets_CA_final.at(0)+jets_CA_final.at(1)).m());
            if(selections[file] > 0 ) { // semilep
                JJmassHP->Fill(VV.m());
                
                //DRJmet->Fill(); -- solve fo the neutrino
            }
        } else if(tau21 > 0.5 && tau21 < 0.75) {
            passcuts=true; counterLP++;
             if(selections[file] == 0 ) JJmassLP->Fill((jets_CA_final.at(0)+jets_CA_final.at(1)).m());
             if(selections[file] > 0 ) JJmassLP->Fill(VV.m());
        }
        if(passcuts){
             J1eta->Fill(jets_CA_final.at(0).eta());
             J1pt->Fill(jets_CA_final.at(0).pt());
            if(selections[file] == 0 ){ 
                DetaJJ->Fill(abs(jets_CA_final.at(0).eta()- jets_CA_final.at(1).eta()));
                J2eta->Fill(jets_CA_final.at(1).eta());
               J2pt->Fill(jets_CA_final.at(1).pt());
            }//close seccond jet
        }// close nsubcuts
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
    double tau21, mprun;
    for(i = 0; i < branchJet->GetEntriesFast(); i++) {
        // take jet and save 4-momentum readable
        jet = (Jet*) branchJet->At(i);
        double pts = jet->PT;
        double etas = jet->Eta;
        double phis=jet->Phi;
        double masse=jet->Mass;
        TLorentzVector jetP4;
        vector<TLorentzVector> jetP4prun;
        jetP4.SetPtEtaPhiM(pts,etas,phis,masse);
        // first how many pass basic selections
        if (jetP4.Pt()>jet_ptmin && jetP4.Eta() < etaj){ 
            jetCol.push_back(jetP4);
            teste.push_back(jetP4);
            JetsBtag.push_back(jet->BTag ); //else JetsFlavour.push_back(99);
            if(i==0){
                tau21 = (jet->Tau2)/(jet->Tau1); 
                //jetP4prun = jet->PrunedP4->X;
            }
            // sub variables still not working in versioned version
            // TLorentzVector PrunedP4[5]; // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta 
            //if (i==0) tau21true = (Jet*) branchJet->Tau2(); 
            //cout<<" prunned mass "<< jet->PrunedP4[0]>>endl;
        } // close if basic cuts
    } // close for each jet
    // composite cuts
    //cout<<tau21<<endl;
    if (jetCol.size() >1) if(abs(jetCol.at(0).Eta()- jetCol.at(1).Eta()) < 1.3 && (jetCol.at(0) + jetCol.at(1)).M() > 890) passcuts = true;
    //if (Jets[0].Pt() < MINPTJET) std::cout << "We have a problem; countJets = " << countJets 
    //				<< "; ijet = " << ijet << " and jet4.Pt() = " << jet4.Pt() << std::endl; 
    //std::cout<<"number of jets "<<countJets<<" "<< branchJet->GetEntriesFast()<<std::endl;
    if(passcuts){
       /* Njets_passing_kLooseID->Fill(jetCol.size(),1);
        JJmass->Fill((jetCol.at(0)+jetCol.at(1)).M());
        DetaJJ->Fill(abs(jetCol.at(0).Eta()- jetCol.at(1).Eta()));
        J1eta->Fill(jetCol.at(0).Eta());
        J2eta->Fill(jetCol.at(1).Eta());
        J1pt->Fill(jetCol.at(0).Pt());
        J2pt->Fill(jetCol.at(1).Pt());
        */
        Nsubint->Fill(tau21);
        //PrunMass->Fill( mprun ); 
                       // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta );
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
    //const char* Mass;
    //const char* Masse = "teste";
    //printf("%d + %d", Mass, Masse);
    //Mass = Form("Control_pp_RSy_WW_jj1.5tev30k%d.root",isample); 
    //cout<<" saved "<< Form("Control_pp_RSy_WW_jj1.5tev30k%d.root",reco)<<endl;
    //string Mass = "teste.root";
    //TFile f1(Mass, "recreate");
    TFile  *f1 = new TFile(Mass[isample], "RECREATE");
    cout<<"saved = "<<Mass[isample]<<endl;
    f1->cd();
    Njets_passing_kLooseID->Write();
    Nlep_passing_kLooseID->Write();
    PrunMass->Write();
    PrunMassSub->Write();
    Nsub->Write();
    Nsubint->Write();
    JJmassHP->Write();
    JJmassLP->Write();
    DetaJJ->Write();
    J1eta->Write();
    J2eta->Write();
    J1pt->Write();
    J2pt->Write();
    DRJl->Write();
    f1->Close();    
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
                                     7, -0.5, 6.5);
    Nlep_passing_kLooseID->GetYaxis()->SetTitle("");
    Nlep_passing_kLooseID->GetXaxis()->SetTitle("Nleptons after showering"); 
    
    PrunMass = new TH1D("PrunMass_ct4",  
                        label, 
                        40, 0, 200);
    PrunMass->GetYaxis()->SetTitle("");
    PrunMass->GetXaxis()->SetTitle("Prunned mass leading jet (GeV)"); 

    PrunMassSub = new TH1D("PrunMassSub_ct4",  
                        label, 
                        40, 0, 200);
    PrunMassSub->GetYaxis()->SetTitle("");
    PrunMassSub->GetXaxis()->SetTitle("Prunned mass sub-leading jet (GeV)"); 
    
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
    
    JJmassHP = new TH1D("JJmassHP_ct4",  
                      label, 
                      50, 0, 3000);
    JJmassHP->GetYaxis()->SetTitle("");
    JJmassHP->GetXaxis()->SetTitle("JJmass (GeV)"); 

    JJmassLP = new TH1D("JJmassLP_ct4",  
                      label, 
                      50, 0, 3000);
    JJmassLP->GetYaxis()->SetTitle("");
    JJmassLP->GetXaxis()->SetTitle("JJmass (GeV)");     

    DRJl = new TH1D("DRJl_ct4",  
                        label, 
                        100, 0, 20);
    DRJl->GetYaxis()->SetTitle("");
    DRJl->GetXaxis()->SetTitle("DR(Jl) (GeV)");  

    DRJmet = new TH1D("DRJmet_ct4",  
                    label, 
                    100, 0, 20);
    DRJmet->GetYaxis()->SetTitle("");
    DRJmet->GetXaxis()->SetTitle("DR(J MET) (GeV)");  

    
    ///////////////////////////////////////////////////////////////////////////////////
    
    PrunMassint = new TH1D("PrunMassint_ct4",  
                        label, 
                        100, 0, 200);
    PrunMassint->GetYaxis()->SetTitle("");
    PrunMassint->GetXaxis()->SetTitle("Prunned mass leading jet (GeV) - interna;"); 
    
    Nsubint = new TH1D("Nsubint_ct4",  
                    label, 
                    50, 0, 1);
    Nsubint->GetYaxis()->SetTitle("");
    Nsubint->GetXaxis()->SetTitle("Nsubjetiness - internal");   
    
    ////////////////////////////////////////////////////////////////////////////
    
    return 0;
}
//////////////////////////////////////////////////
