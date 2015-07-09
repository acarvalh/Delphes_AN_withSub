#include <iostream>
#include <fstream>
#include <string>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <vector>
// Delphes headers
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
// fastjet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/PseudoJet.hh"
 
 GenParticle *particle;
 Electron *electron;
 Photon *photon;
 Muon *muon;
 MissingET *met;
 Track *track;
 Tower *tower;
 TObject *object;

 Float_t Eem, Ehad;
 Bool_t skip;
 Long64_t entry;
 Int_t i, j, pdgCode;
