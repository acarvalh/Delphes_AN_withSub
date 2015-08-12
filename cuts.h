/////////////////////////////////
const int nfiles=8; 
// selections: 0 = CMS JJ ; 1 = CMS Jlnu ; 2 = CMS Jll 
// selections: 3 = ATLAS JJ ; 4 = ATLAS Jlnu ; 5 = ATLAS Jll 

string file[8] = {
    "bulk/pp_bulk_WW_jj1.5tev30k_PU22_versioned", 
    "bulk/pp_bulk_ZZ_jj1.5tev20k_PU22_versioned", 
    "bulk/pp_bulk_WW_semilep1.5tev30k_PU22_versioned",
    "bulk/pp_bulk_ZZ_llj1.5tev30k_PU22_versioned",
    "RS/pp_RS_WW_jj1.5tev20k_PU22_versioned",
    "RS/pp_RS_ZZ_jj1.5tev20k_PU22_versioned",
    "RS/pp_RS_WW_jlnu1.5tev20k_PU22_versioned",    
    "RS/pp_RS_ZZ_jll1.5tev20k_PU22_versioned"
};
const char* Mass[8] = {
    "Control_pp_bulk_WW_jj1.5tev30k_PU22_versioned.root", 
    "Control_pp_bulk_ZZ_jj1.5tev20k_PU22_versioned.root", 
    "Control_pp_bulk_WW_semilep1.5tev30k_PU22_versioned.root",
    "Control_pp_bulk_ZZ_llj1.5tev30k_PU22_versioned.root",
    "Control_pp_RS_WW_jj1.5tev20k_PU22_versioned.root",
    "Control_pp_RS_ZZ_jj1.5tev20k_PU22_versioned.root",
    "Control_pp_RS_WW_jlnu1.5tev20k_PU22_versioned.root",    
    "Control_pp_RS_ZZ_jll1.5tev20k_PU22_versioned.root"
};
// double const selections[6] = {0 ,1,2,0,1,2}; // see cuts.h
 double const selections[8] = {3,3 ,4,5,0,0,1,2}; // see cuts.h
/////////////////////////////////
// cuts
//double weight =1.;///10000;//0.001;//
bool shower=true;
// To be applied only to hadron level events
//
// objets deffinition
double const jet_ptmin=30.0; // parton for jet reconstruction
double const const_ptmin=0; // parton for jet reconstruction
double const rapmax=2.5; // for jet reconstruction
double const etaj = 2.5;
double const etal = 2.5;
double const ptlepton = 20.;
double const lepiso = 0.8;
////////////////////////////////////
// analysis cuts
// CMS JJ - selection 0
double const Deltay = 1.3;//3;
double const Mvvjj =0;//890;//400;
// CMS Jlnu - selection 1
double const etaJ = 2.4; 
double const Mvvlnuj =700;//400; 
double const METenuJ =40;//400; 
double const METmunuJ =80;//400;
double const ptVlnu = 200; 
double const ptE = 90; 
double const ptMu = 50; 
double const etaE = 2.5; 
double const etaMu = 2.1; 
// CMS Jll - selection 2
double const Mvvllj =500;//400; 
double const ptVll = 80; 
///////////////////////////////////
// for substructure
// CMS
// nsubjetiness - unnormalized measure
double const beta1 = 1.0;
double const R0=0.8; 
double const tau21_LP = 0.5;
double const tau21_HP = 0.75;
// prunning
double const zcut = 0.1;
double const Rcut_factor =0.8;
//
double const mprunjj_min = 70;
double const mprunjj_max = 100;
//
double const mprunjlnu_min = 65;
double const mprunjlnu_max = 105;
//
double const mprunjll_min = 70;
double const mprunjll_max = 110;
//////////////////////////////
// mass drop
// ATLAS
double const R0A=1.2; 
double const Rr = 0.3; // CA recluster
double const mu = 1;//0.67; // mu: ratio in between mass of cores, symetric splitting
double const ycut = 0.04 ;// 0.09;
double const nf =3; // n subjets after filtering
//
//double const Rfilt = 0.1;
//int const n_subjet =3;
//
/////////////////////////////////////////
// smear
bool smear= false;
double const mean =1, resolution = 0.2;
bool JES = true;
double const JESup = 0.1; 
//////////////////////////////////////////
// on the wwbb analysis
double const wmass = 80.4;
double const bmass = 4.7;
double const tmass = 173.0;
double const higgs_mass = 125.0;
double const MeeMax = 30000.0;
double const MetMin = 0.0;
double const  wbtransmassMax = 3000;
double const  wbtransmassMin = 0;