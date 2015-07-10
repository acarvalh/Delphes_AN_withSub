/////////////////////////////////
const int nfiles=5; 
string file[5] = {"pp_RS_WW_jj1.5tev30k", 
    "pp_RS_ww_semilep_1.5tev30k",
    "pp_bulk_WW_jj1.5tev30k",
    "pp_bulk_ww_semilep_1.5tev30k",
    "pp_bulk_WW_jj1.5tev30k_PU22_uncersioned"
};
const char* Mass[5] = {"Control_pp_RS_WW_jj1.5tev30k.root", 
    "Control_pp_RS_ww_semilep_1.5tev30k.root",
    "Control_pp_bulk_WW_jj1.5tev30k.root",
    "Control_pp_bulk_ww_semilep_1.5tev30k.root",
    "pp_bulk_WW_jj1.5tev30k_PU22_uncersioned.root"
};
/////////////////////////////////
// cuts
//double weight =1.;///10000;//0.001;//
double Mjj =0;//400; 
double PTjj = 0;//400; 
double DeltayVBF = 0;//3;
double DeltaRVBF = 0;//3;
bool shower=true;
// To be applied only to hadron level events
//
double const genmasshad=2000; // genmass
double const genmasshadmin=0; // genmass
//
double const genmasslep=2000; // genmass
double const genmasslepmin=0; // genmass
// the gen-level cuts
double const bjetpt = 30.0; 
double const mbblow = 0.0; 
// basline
double const jet_ptmin=20.0; // parton for jet reconstruction
double const jet_ptminfinal=1.0; // in final jet reconstruction
double const rapmax=5.0; // for jet reconstruction
double const etab = 2.5;
double const etal = 25;
double const etaj=5;
double const RR =0.5;
double const ptlepton = 20.;
double const lepiso = 0.3;
double Deltay = 1.3;//3; // not yet
// analysis cuts
double const mblcut = 190;
int const cat =2; // minimum number of btag
////////////////////////////////////////
// weights b-tag
double const subjet2b=1;
double const fatjet2b=1;
double const normalb=1;
double const normalc=1;
double const normall=1;
double const misb=1;
/////////////////////////////////////////
///////////////////////////////////
// for substructure
// mass drop
double const Rsb = 0.8; // CA recluster
double const mu = 0.67;
double const ycut = 0.09;
double const Mfat =100;
//
double const Rfilt = 0.1;
int const n_subjet =3;
//
double const zcut = 0.1;
double const Rcut_factor =0.8;
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