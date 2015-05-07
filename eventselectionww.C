#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_GenEvent.h"

#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TDataType.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "/home/dan-bia/workspace/Build_clhep/include/CLHEP/Vector/LorentzVector.h"
//#include "/home/dan-bia/workspace/x86_64-slc5-gcc47-opt/include/CLHEP/Vector/LorentzVector.h"

//#include "CLHEP/Vector/LorentzVector.h"
#include <math.h>
#include <algorithm>
#include <list>

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

using namespace CLHEP;

/*
//==============================================================
11      -11     e-         e+
12      -12    nu_e       nu_ebar
13      -13     mu-        mu+
14      -14    nu_mu      nu_mubar
15      -15    tau-       tau+
16      -16    nu_tau     nu_taubar
17      -17    tau'-      tau'+
18      -18    nu'_tau    nu'_taubar
24      -24     W+         W-
37      -37     H+         H-
23       22     Z0         gamma
//===========================================================
*/


//===================================
// Declaration of auxiliary classes
//===================================
class IsW_Boson {
public:
  /// returns true if the GenParticle is a W
  bool operator()( const HepMC::GenParticle* p ) { 
    if ( p->pdg_id() == 24 ) return 1;
    if ( p->pdg_id() == -24 ) return 1;
    return 0;
  }
};

class Is_Photon {
public:
  /// returns true if the GenParticle is a W
  bool operator()( const HepMC::GenParticle* p ) { 
    if( abs( p->pdg_id() == 22)) return 1;
    return 0;
  }
};

class IsStateFinal {
public:
  /// returns true if the GenParticle does not decay
  bool operator()( const HepMC::GenParticle* p ) { 
    if ( !p->end_vertex() && p->status()==1 ) return 1;
    return 0;
  }
};

//===================================
//Event Selection
//===================================

class IsEventGood {
public:
  /// check this event for goodness
  bool operator()( const HepMC::GenEvent* evt ) { 
    for( HepMC::GenEvent::particle_const_iterator p 
	   = evt->particles_begin(); p != evt->particles_end(); ++p ){
      if ( !(*p)->end_vertex() && (*p)->status()==1 && 
	   abs((*p)->pdg_id()) == 13) {
        std::cout << "Event " << evt->event_number()
		  << " is a good event." << std::endl;
	(*p)->print();
        return 1;}
      
      else if ( !(*p)->end_vertex() && (*p)->status()==1 && 
		abs((*p)->pdg_id()) == 11){
        std::cout << "Event " << evt->event_number()
		  << " is a good event." << std::endl;
        (*p)->print();
        return 1;}
    }
    return 0;
  }
};

//========================================
//Análise principal
//=========================================


int main() { 
  
  //{ // begin scope of ascii_in and ascii_out
    HepMC::IO_GenEvent ascii_in("aa_wwa0wp150fpmc_13tev.dat",std::ios::in);
    
    // declare another IO_GenEvent for writing out the good events
    HepMC::IO_GenEvent ascii_out("outpomreg_wwhepmc_13tev.dat",std::ios::out);
    // Build HepPDT particle table
    const char infile[] ="/home/dan-bia/workspace/teste_heppdt/data/particle.tbl";   
    std::ifstream pdfile( infile );
    if( !pdfile ) { 
      std::cerr << ">>> Cannot open " << infile << std::endl;
      exit(-1);
    }
    HepPDT::ParticleDataTable pdt( "Particle Table" );
    {
      // Construct table builder
      HepPDT::TableBuilder tb(pdt);
      if( !addParticleTable( pdfile, tb, true ) ) { 
	std::cout << ">> Error reading PDG pdt file " << std::endl; 
      }
    } // the tb destructor fills datacol
    // Loop over Particle Data Table
    std::ostringstream oss;
    oss << std::setw(15) << "Particle Id"
	<< std::setw(22) << "Particle Name"
	<< std::setw(15) << "Three-charge" << std::endl; 
    for( HepPDT::ParticleDataTable::const_iterator p = pdt.begin(); p != pdt.end(); ++p ) {
      const HepPDT::ParticleID & id = p->first;
      int pdgId = id.pid();
      int q3 = id.threeCharge();
      //double q = id.charge();
      const std::string& name = id.PDTname();
      oss << std::setw(15) << pdgId
	  << std::setw(22) << name
	  << std::setw(15) << q3 << std::endl;
    }
    std::cout << oss.str();
    
    // Declare ROOT TTree
    
    //mu-
    int const NMUMAX = 50;
    int n_mu = 0;
    double mu_pt[NMUMAX];
    double mu_px[NMUMAX];
    double mu_py[NMUMAX];
    double mu_pz[NMUMAX];
    double mu_eta[NMUMAX];
    double mu_phi[NMUMAX];
    double mu_energy[NMUMAX];
    double mu_mass[NMUMAX];
    
    int const NVMUMAX=50;
    int n_vmu = 0;
    double vmu_pt[NVMUMAX];
    double vmu_px[NVMUMAX];
    double vmu_py[NVMUMAX];
    double vmu_pz[NVMUMAX];
    double vmu_eta[NVMUMAX];
    double vmu_phi[NVMUMAX];
    double vmu_energy[NVMUMAX];
    double vmu_mass[NVMUMAX];
    
    int const NEMAX = 50;
    int n_e = 0;
    double e_pt[NEMAX];
    double e_px[NEMAX];
    double e_py[NEMAX];
    double e_pz[NEMAX];
    double e_eta[NEMAX];
    double e_phi[NEMAX];
    double e_energy[NEMAX];
    double e_mass[NEMAX];
    
    int const NVEMAX = 50;
    int n_ve = 0;
    double ve_pt[NVEMAX];
    double ve_px[NVEMAX];
    double ve_py[NVEMAX];
    double ve_pz[NVEMAX];
    double ve_eta[NVEMAX];
    double ve_phi[NVEMAX];
    double ve_energy[NVEMAX];
    double ve_mass[NVEMAX];
    
    int const NCHGMAX =300;  
    int n_chg = 0;
    int chg_id[NCHGMAX];
    double chg_pt[NCHGMAX];
    double chg_px[NCHGMAX];
    double chg_py[NCHGMAX];
    double chg_pz[NCHGMAX];
    double chg_ch[NCHGMAX];
    double chg_eta[NCHGMAX];
    double chg_phi[NCHGMAX];
    double chg_energy[NCHGMAX];
    double chg_mass[NCHGMAX];
    
    int const NW1MAX = 50;
    int n_w1 =0;
    double w1_pt[NW1MAX];
    double w1_px[NW1MAX];
    double w1_py[NW1MAX];
    double w1_pz[NW1MAX];
    double w1_mass[NW1MAX];
    double w1_energy[NW1MAX];
    double w1_transversemass[NW1MAX];
    double w1_phi[NW1MAX];
    double w1_cos[NW1MAX];
    double w1_eta[NW1MAX];
    
    int w1_lep_id[NW1MAX];
    double w1_lep_px[NW1MAX];
    double w1_lep_py[NW1MAX];
    double w1_lep_pz[NW1MAX];
    double w1_lep_energy[NW1MAX];

    int w1_nu_id[NW1MAX];
    double w1_nu_px[NW1MAX];
    double w1_nu_py[NW1MAX];
    double w1_nu_pz[NW1MAX];
    double w1_nu_energy[NW1MAX];
    
    int const NW2MAX = 50;
    int n_w2 =0;
    double w2_pt[NW2MAX];
    double w2_px[NW2MAX];
    double w2_py[NW2MAX];
    double w2_pz[NW2MAX];
    double w2_energy[NW2MAX];
    double w2_mass[NW2MAX];
    double w2_transversemass[NW2MAX];
    double w2_phi[NW2MAX];
    double w2_cos[NW2MAX];
    double w2_eta[NW2MAX];
    
    int w2_lep_id[NW2MAX];
    double w2_lep_px[NW2MAX];
    double w2_lep_py[NW2MAX];
    double w2_lep_pz[NW2MAX];
    double w2_lep_energy[NW2MAX];

    int w2_nu_id[NW2MAX];
    double w2_nu_px[NW2MAX];
    double w2_nu_py[NW2MAX];
    double w2_nu_pz[NW2MAX];
    double w2_nu_energy[NW2MAX];

    int const GAMAMAX = 50;
    int n_a =0;
    double a_pt[GAMAMAX];
    double a_px[GAMAMAX];
    double a_fracmom[GAMAMAX];
    double a_pini[GAMAMAX];
    double a_py[GAMAMAX];
    double a_pz[GAMAMAX];
    double a_mass[GAMAMAX];
    double a_energy[GAMAMAX];



    int const PROTONMAX = 100;
    int n_proton =0;
    double proton_pt[PROTONMAX];
    double proton_px[PROTONMAX];
    double proton_py[PROTONMAX];
    double proton_pz[PROTONMAX];
    double proton_energy[PROTONMAX];


    TTree* T = new TTree("T","Tree");
    T->Branch("n_mu", &n_mu,"n_mu/I");
    T->Branch("mu_pt", &mu_pt,"mu_pt[n_mu]/D");
    T->Branch("mu_px", &mu_px,"mu_px[n_mu]/D");
    T->Branch("mu_py", &mu_py,"mu_py[n_mu]/D"); 
    T->Branch("mu_pz", &mu_pz,"mu_pz[n_mu]/D");
    T->Branch("mu_eta", &mu_eta,"mu_eta[n_mu]/D");
    T->Branch("mu_phi", &mu_phi,"mu_phi[n_mu]/D");
    T->Branch("mu_energy", &mu_energy,"mu_energy[n_mu]/D");
    T->Branch("mu_mass", &mu_mass,"mu_mass[n_mu]/D");

    //Neutrino muon
    T->Branch("n_vmu", &n_vmu,"n_vmu/I");
    T->Branch("vmu_pt", &vmu_pt,"vmu_pt[n_vmu]/D");
    T->Branch("vmu_px", &vmu_px,"vmu_px[n_vmu]/D");
    T->Branch("vmu_py", &vmu_py,"vmu_py[n_vmu]/D");
    T->Branch("vmu_pz", &vmu_pz,"vmu_pz[n_vmu]/D");
    T->Branch("vmu_eta", &vmu_eta,"vmu_eta[n_vmu]/D");
    T->Branch("vmu_phi", &vmu_phi,"vmu_phi[n_vmu]/D");
    T->Branch("vmu_energy", &vmu_energy,"vmu_energy[n_vmu]/D");
    T->Branch("vmu_mass", &vmu_mass,"vmu_mass[n_vmu]/D");
    //Neutrino
    
    T->Branch("n_e", &n_e,"n_e/I");
    T->Branch("e_pt", &e_pt,"e_pt[n_e]/D");
    T->Branch("e_px", &e_px,"e_px[n_e]/D");
    T->Branch("e_py", &e_py,"e_py[n_e]/D");
    T->Branch("e_pz", &e_pz,"e_pz[n_e]/D");
    T->Branch("e_eta", &e_eta,"e_eta[n_e]/D");
    T->Branch("e_phi", &e_phi,"e_phi[n_e]/D");
    T->Branch("e_energy", &e_energy,"e_energy[n_e]/D");
    T->Branch("e_mass", &e_mass,"e_mass[n_e]/D");

    //Neutrino muon
    T->Branch("n_ve", &n_ve,"n_ve/I");
    T->Branch("ve_pt", &ve_pt,"ve_pt[n_ve]/D");
    T->Branch("ve_px", &ve_px,"ve_px[n_ve]/D");
    T->Branch("ve_py", &ve_py,"ve_py[n_ve]/D");
    T->Branch("ve_pz", &ve_pz,"ve_pz[n_ve]/D");
    T->Branch("ve_eta", &ve_eta,"ve_eta[n_ve]/D");
    T->Branch("ve_phi", &ve_phi,"ve_phi[n_ve]/D");
    T->Branch("ve_energy", &ve_energy,"ve_energy[n_ve]/D");
    T->Branch("ve_mass", &ve_mass,"ve_mass[n_ve]/D");

    T->Branch("n_chg", &n_chg,"n_chg/I");
    T->Branch("chg_id", &chg_id,"chg_id[n_chg]/I");
    T->Branch("chg_ch", &chg_ch,"chg_ch[n_chg]/D");
    T->Branch("chg_px", &chg_px,"chg_px[n_chg]/D");
    T->Branch("chg_py", &chg_py,"chg_py[n_chg]/D");
    T->Branch("chg_pz", &chg_pz,"chg_pz[n_chg]/D");
    T->Branch("chg_pt", &chg_pt,"chg_pt[n_chg]/D");
    T->Branch("chg_eta", &chg_eta,"chg_eta[n_chg]/D");
    T->Branch("chg_phi", &chg_phi,"chg_phi[n_chg]/D");
    T->Branch("chg_energy", &chg_energy,"chg_energy[n_chg]/D");
    T->Branch("chg_mass", &chg_mass,"chg_mass[n_chg]/D");

    //W1 = W+
    T->Branch("n_w1", &n_w1,"n_w1/I");
    T->Branch("w1_pt", &w1_pt,"w1_pt[n_w1]/D");
    T->Branch("w1_px", &w1_px,"w1_px[n_w1]/D");
    T->Branch("w1_py", &w1_py,"w1_py[n_w1]/D");
    T->Branch("w1_pz", &w1_pz,"w1_pz[n_w1]/D");
    T->Branch("w1_mass", &w1_mass,"w1_mass[n_w1]/D");
    T->Branch("w1_energy", &w1_energy,"w1_energy[n_w1]/D");
    T->Branch("w1_transversemass", &w1_transversemass,"w_transversemass[n_w1]/D");
    T->Branch("w1_phi", &w1_phi,"w1_phi[n_w1]/D");
    T->Branch("w1_cos", &w1_cos,"w1_cos[n_w1]/D");
    T->Branch("w1_lep_id", &w1_lep_id,"w1_lep_id[n_w1]/I");
    T->Branch("w1_lep_px", &w1_lep_px,"w1_lep_px[n_w1]/D");
    T->Branch("w1_lep_py", &w1_lep_py,"w1_lep_py[n_w1]/D");
    T->Branch("w1_lep_pz", &w1_lep_pz,"w1_lep_pz[n_w1]/D");
    T->Branch("w1_lep_energy", &w1_lep_energy,"w1_lep_energy[n_w1]/D");
    T->Branch("w1_nu_id", &w1_nu_id,"w1_nu_id[n_w1]/I");
    T->Branch("w1_nu_px", &w1_nu_px,"w1_nu_px[n_w1]/D");
    T->Branch("w1_nu_py", &w1_nu_py,"w1_nu_py[n_w1]/D");
    T->Branch("w1_nu_pz", &w1_nu_pz,"w1_nu_pz[n_w1]/D");
    T->Branch("w1_nu_energy", &w1_nu_energy,"w1_nu_energy[n_w1]/D");

    // W2=W-
    T->Branch("n_w2", &n_w2,"n_w2/I");
    T->Branch("w2_pt", &w2_pt,"w2_pt[n_w2]/D");
    T->Branch("w2_px", &w2_px,"w2_px[n_w2]/D");
    T->Branch("w2_py", &w2_py,"w2_py[n_w2]/D");
    T->Branch("w2_pz", &w2_pz,"w2_pz[n_w2]/D");
    T->Branch("w2_energy", &w2_energy,"w2_energy[n_w2]/D");
    T->Branch("w2_mass", &w2_mass,"w2_mass[n_w2]/D");
    T->Branch("w2_transversemass", &w2_transversemass,"w2_transversemass[n_w2]/D");
    T->Branch("w2_phi", &w2_phi,"w2_phi[n_w2]/D");
    T->Branch("w2_cos", &w2_cos,"w2_cos[n_w2]/D");
    T->Branch("w2_lep_id", &w2_lep_id,"w2_lep_id[n_w2]/I");
    T->Branch("w2_lep_px", &w2_lep_px,"w2_lep_px[n_w2]/D");
    T->Branch("w2_lep_py", &w2_lep_py,"w2_lep_py[n_w2]/D");
    T->Branch("w2_lep_pz", &w2_lep_pz,"w2_lep_pz[n_w2]/D");
    T->Branch("w2_lep_energy", &w2_lep_energy,"w2_lep_energy[n_w2]/D");
    T->Branch("w2_nu_id", &w2_nu_id,"w2_nu_id[n_w2]/I");
    T->Branch("w2_nu_px", &w2_nu_px,"w2_nu_px[n_w2]/D");
    T->Branch("w2_nu_py", &w2_nu_py,"w2_nu_py[n_w2]/D");
    T->Branch("w2_nu_pz", &w2_nu_pz,"w2_nu_pz[n_w2]/D");
    T->Branch("w2_nu_energy", &w2_nu_energy,"w2_nu_energy[n_w2]/D");
      
    //Gamma
    T->Branch("n_a", &n_a,"n_a/I");
    T->Branch("a_pt", &a_pt,"a_pt[n_a]/D");
    T->Branch("a_px", &a_px,"a_px[n_a]/D");
    T->Branch("a_fracmom", &a_fracmom,"a_framom[n_a]/D");
    T->Branch("a_py", &a_py,"a_py[n_a]/D");
    T->Branch("a_pini", &a_pini,"a_pini[n_a]/D");
    T->Branch("a_pz", &a_pz,"a_pz[n_a]/D");
    T->Branch("a_energy", &a_energy,"a_energy[n_a]/D");
    T->Branch("a_mass", &a_mass,"a_mass[n_a]/D");


  //Gamma
    T->Branch("n_proton", &n_proton,"n_proton/I");
    T->Branch("proton_px", &proton_px,"proton_px[n_proton]/D");
    T->Branch("proton_py", &proton_py,"proton_py[n_proton]/D");
    T->Branch("proton_pz", &proton_pz,"proton_pz[n_proton]/D");
    T->Branch("proton_pt", &proton_pt,"proton_pt[n_proton]/D");
    T->Branch("proton_energy", &proton_energy,"proton_energy[n_proton]/D");
    


    // Declarar os histogramas no ROOT 

    TH1F* h_na = new TH1F("na","na",100,0.,50.);
    TH1F* h_a_pt = new TH1F("a_pt","a_pt",100,0.,1000.);
    TH1F* h_a_px = new TH1F("a_px","a_px",100,0.,1000.);
    TH1F* h_a_fracmom = new TH1F("a_fracmom","a_fracmom",100,0.,1.);
    TH1F* h_a_py = new TH1F("a_py","a_py",100,0.,1000.);
    TH1F* h_a_pini = new TH1F("a_pini","a_pini",100,-2.,2.);
    TH1F* h_a_pz = new TH1F("a_pz","a_pz",100,0.,1000.);
    TH1F* h_a_energy = new TH1F("a_energy","a_energy",100,0.,1000.);
    TH1F* h_a_mass = new TH1F("a_mass","a_mass",100,0.,1);

    TH1F* h_nmu = new TH1F("nmu","nmu",100,0.,50.);
    TH1F* h_mu_pt = new TH1F("mu_pt","mu_pt",100,0.,1000.);
    TH1F* h_mu_px = new TH1F("mu_px","mu_px",100,0.,1000.);
    TH1F* h_mu_py = new TH1F("mu_py","mu_py",100,0.,1000.);
    TH1F* h_mu_pz = new TH1F("mu_pz","mu_pz",100,0.,1000.);
    TH1F* h_mu_eta = new TH1F("mu_eta","mu_eta",100,-5.0,5.0);
    TH1F* h_mu_phi = new TH1F("mu_phi","mu_phi",100,-M_PI,M_PI);
    TH1F* h_mu_energy = new TH1F("mu_energy","mu_energy",100,0.,1000.);
    TH1F* h_mu_mass = new TH1F("mu_mass","mu_mass",100,0.,140.);

    TH1F* h_nvmu = new TH1F("nvmu","nvmu",100,0.,50.);
    TH1F* h_vmu_pt = new TH1F("vmu_pt","vmu_pt",100,0.,1000.);
    TH1F* h_vmu_px = new TH1F("vmu_px","vmu_px",100,0.,1000.);
    TH1F* h_vmu_py = new TH1F("vmu_py","vmu_py",100,0.,1000.);
    TH1F* h_vmu_pz = new TH1F("vmu_pz","vmu_pz",100,0.,1000.);
    TH1F* h_vmu_eta = new TH1F("vmu_eta","vmu_eta",100,-5.0,5.0);
    TH1F* h_vmu_phi = new TH1F("vmu_phi","vmu_phi",100,-M_PI,M_PI);
    TH1F* h_vmu_energy = new TH1F("vmu_energy","vmu_energy",100,0.,1000.);
    TH1F* h_vmu_mass = new TH1F("vmu_mass","vmu_mass",100,0.,140.);

    //Eletrons e neutrinos dos eletrons
    TH1F* h_ne = new TH1F("ne","ne",100,0.,50.);
    TH1F* h_e_pt = new TH1F("e_pt","e_pt",100,0.,1000.);
    TH1F* h_e_px = new TH1F("e_px","e_px",100,0.,1000.);
    TH1F* h_e_py = new TH1F("e_py","e_py",100,0.,1000.);
    TH1F* h_e_pz = new TH1F("e_pz","e_pz",100,0.,1000.);
    TH1F* h_e_eta = new TH1F("e_eta","e_eta",100,-5.0,5.0);
    TH1F* h_e_phi = new TH1F("e_phi","e_phi",100,-M_PI,M_PI);
    TH1F* h_e_energy = new TH1F("e_energy","e_energy",100,0.,1000.);
    TH1F* h_e_mass = new TH1F("e_mass","e_mass",100,0.,140.);

    TH1F* h_nve = new TH1F("nve","nve",100,0.,50.);
    TH1F* h_ve_pt = new TH1F("ve_pt","ve_pt",100,0.,1000.);
    TH1F* h_ve_px = new TH1F("ve_px","ve_px",100,0.,1000.);
    TH1F* h_ve_py = new TH1F("ve_py","ve_py",100,0.,1000.);
    TH1F* h_ve_pz = new TH1F("ve_pz","ve_pz",100,0.,1000.);
    TH1F* h_ve_eta = new TH1F("ve_eta","ve_eta",100,-5.0,5.0);
    TH1F* h_ve_phi = new TH1F("ve_phi","ve_phi",100,-M_PI,M_PI);
    TH1F* h_ve_energy = new TH1F("ve_energy","ve_energy",100,0.,1000.);
    TH1F* h_ve_mass = new TH1F("ve_mass","ve_mass",100,0.,140.);

    //Particulas carregadas

    TH1F* h_nchg = new TH1F("nchg","nchg",100,0.,100.);
    TH1F* h_chg_id = new TH1F("chg_id","chg_id",4000,-2000.,2000.);
    TH1F* h_chg_ch = new TH1F("chg_ch","chg_ch",1000,-1000.,1000.);
    TH1F* h_chg_pt = new TH1F("chg_pt","chg_pt",1000,0.,1000.);
    TH1F* h_chg_px = new TH1F("chg_px","chg_px",1000,0.,1000.);
    TH1F* h_chg_py = new TH1F("chg_py","chg_py",1000,0.,1000.);
    TH1F* h_chg_pz = new TH1F("chg_pz","chg_pz",1000,0.,1000.);
    TH1F* h_chg_eta = new TH1F("chg_eta","chg_eta",100,-5.0,5.0);
    TH1F* h_chg_phi = new TH1F("chg_phi","chg_phi",100,-M_PI,M_PI);
    TH1F* h_chg_energy = new TH1F("chg_energy","chg_energy",1000,0.,1000.);
    TH1F* h_chg_mass = new TH1F("chg_mass","chg_mass",100,0.,140.);




    TH1F* h_nproton = new TH1F("nproton","nproton",100,0.,50.);
    TH1F* h_proton_pt = new TH1F("proton_pt","proton_pt",100,0.,1000.);
    TH1F* h_proton_px = new TH1F("proton_px","proton_px",100,0.,1000.);
    TH1F* h_proton_py = new TH1F("proton_py","proton_py",100,0.,1000.);
    TH1F* h_proton_pz = new TH1F("proton_pz","proton_pz",100,0.,1000.);
    TH1F* h_proton_energy = new TH1F("proton_energy","proton_energy",1000,0.,1000.);
  

    TH1F* h_nw1 = new TH1F("nw1","nw1",100,0.,50.);
    TH1F* h_w1_pt = new TH1F("w1_pt","w1_pt",100,0.,1000.);
    TH1F* h_w1_px = new TH1F("w1_px","w1_px",100,0.,1000.);
    TH1F* h_w1_py = new TH1F("w1_py","w1_py",100,0.,1000.);
    TH1F* h_w1_pz = new TH1F("w1_pz","w1_pz",100,0.,1000.);
    TH1F* h_w1_energy = new TH1F("w1_energy","w1_energy",100,0.,1000.);
    TH1F* h_w1_mass = new TH1F("w1_mass","w1_mass",100,0.,500.);
    TH1F* h_w1_transversemass = new TH1F("w1_transversemass","w1_transversemass",100,0.,500.);
    TH1F* h_w1_phi = new TH1F("w1_phi","w1_phi",100,-M_PI,M_PI);
    TH1F* h_w1_cos = new TH1F("w1_cos","w1_cos",100,-5,5);
    TH1F* h_w1_eta = new TH1F("w1_eta","w1_eta",100,-5,5);


    TH1F* h_nw2 = new TH1F("nw2","nw2",100,0.,50.);
    TH1F* h_w2_pt = new TH1F("w2_pt","w2_pt",100,0.,1000.);
    TH1F* h_w2_px = new TH1F("w2_px","w2_px",100,0.,1000.);
    TH1F* h_w2_py = new TH1F("w2_py","w2_py",100,0.,1000.);
    TH1F* h_w2_pz = new TH1F("w2_pz","w2_pz",100,0.,1000.);
    TH1F* h_w2_energy = new TH1F("w2_energy","w2_energy",100,0.,1000.);
    TH1F* h_w2_mass = new TH1F("w2_mass","w2_mass",100,0.,500.);
    TH1F* h_w2_transversemass = new TH1F("w2_transversemass","w2_transversemass",100,0.,500.);
    TH1F* h_w2_phi = new TH1F("w2_phi","w2_phi",100,-M_PI,M_PI);
    TH1F* h_w2_cos = new TH1F("w2_cos","w2_cos",100,-5,5);
    TH1F* h_w2_eta = new TH1F("w2_eta","w2_eta",100,-5,5);
 
    TH1F* h_2lep2neu = new TH1F("2lep2neu","2lep2neu",100,0,100);


    // EVENT LOOP
    // declare an instance of the event selection predicate
    IsEventGood is_good_event;
    IsStateFinal isfinal;
    int icount=0;
    int num_good_events=0;
    HepMC::GenEvent* evt = ascii_in.read_next_event();

   
    while ( evt ) {
      icount++;
      if ( icount%50==1 ) std::cout << "Processing Event Number " << icount
				    << " its # " << evt->event_number() 
				    << std::endl;
      // Reset Tree variables per event
      n_mu = 0;
      for(int imu = 0; imu < NMUMAX; ++imu) {
	mu_pt[imu] =  -999.;
	mu_px[imu] =  -999.;
	mu_py[imu] =  -999.; 
	mu_pz[imu] =  -999.;                                                    	mu_eta[imu] = -999.;
	mu_phi[imu] = -999.;
	mu_energy[imu] = -999.;
	mu_mass[imu] = -999.;}

      n_vmu=0; 
      for(int ivmu = 0; ivmu < NVMUMAX; ++ivmu) {
	vmu_pt[ivmu] = -999.;
	vmu_px[ivmu] = -999.;
	vmu_py[ivmu] = -999.;
	vmu_pz[ivmu] = -999.;
	vmu_eta[ivmu] = -999.;
	vmu_phi[ivmu] = -999.;
	vmu_energy[ivmu] = -999.;
	vmu_mass[ivmu] = -999.; }
      
      n_e = 0;
      for(int ie = 0; ie < NEMAX; ++ie) { e_pt[ie] = -999.;
	e_eta[ie] = -999.;
	e_px[ie] = -999.;
	e_py[ie] = -999.;
	e_pz[ie] = -999.;
	e_phi[ie] = -999.; 
	e_energy[ie] = -999.;
	e_mass[ie] = -999.;}

      n_ve=0; 
      for(int ive = 0; ive < NVEMAX; ++ive) { ve_pt[ive] = -999.;
	ve_px[ive] = -999.;
	ve_py[ive] = -999.;
	ve_pz[ive] = -999.;
	ve_eta[ive] = -999.;
	ve_phi[ive] = -999.;
	ve_energy[ive] = -999.;
	ve_mass[ive] = -999.; }

      n_a=0; 
      for(int ia = 0; ia < GAMAMAX; ++ia) { a_pt[ia] = -999.;
	a_px[ia] = -999.;
	a_py[ia] = -999.;
	a_pz[ia] = -999.;
	a_fracmom[ia] = -999.;
	a_pini[ia] = -999.;
	a_energy[ia] = -999.;
	a_mass[ia] = -999.; }

      n_chg = 0;
      for(int ichg = 0; ichg < NCHGMAX; ++ichg) {chg_id[ichg] = -999;
	chg_ch[ichg] = -999.;
	chg_pt[ichg] = -999.;
	chg_px[ichg] = -999.;
	chg_py[ichg] = -999.;
	chg_pz[ichg] = -999.;
	chg_eta[ichg] = -999.;
	chg_phi[ichg] = -999.;
	chg_energy[ichg] = -999.;
	chg_mass[ichg] = -999.;}

      n_w1 = 0;
      for(int iw1 = 0; iw1 < NW1MAX; ++iw1){ w1_pt[iw1] = -999;
	w1_px[iw1] = -999;
	w1_py[iw1] = -999;
	w1_pz[iw1] = -999;
	w1_mass[iw1] = -999;
	w1_transversemass[iw1] = -999;
	w1_phi[iw1] = -999;
	w1_cos[iw1] = -999;
	w1_eta[iw1] = -999;
	w1_energy[iw1] = -999;
	w1_lep_id[iw1] = -999;
	w1_lep_px[iw1] = -999;
	w1_lep_py[iw1] = -999;
	w1_lep_pz[iw1] = -999;
	w1_lep_energy[iw1] = -999;

	w1_nu_id[iw1] = -999;
	w1_nu_px[iw1] = -999;
	w1_nu_py[iw1] = -999;
	w1_nu_pz[iw1] = -999;
	w1_nu_energy[iw1] = -999;}

      n_w2 = 0;
      for(int iw2 =0; iw2 < NW2MAX ; ++iw2){ w2_pt[iw2] = -999;
	w2_px[iw2] = -999;
	w2_py[iw2] = -999;
	w2_pz[iw2] = -999;
	w2_energy[iw2] = -999;
	w2_mass[iw2] = -999;
	w2_transversemass[iw2] =-999;
	w2_phi[iw2] = -999;
	w2_cos[iw2] = -999;
	w2_eta[iw2] = -999;
	w2_lep_id[iw2] = -999;
	w2_lep_px[iw2] = -999;
	w2_lep_py[iw2] = -999;
	w2_lep_pz[iw2] = -999;
	w2_lep_energy[iw2] = -999;

	w2_nu_id[iw2] = -999;
	w2_nu_px[iw2] = -999;
	w2_nu_py[iw2] = -999;
	w2_nu_pz[iw2] = -999;
	w2_nu_energy[iw2] = -999;}

      n_proton=0; 
      for(int iproton = 0; iproton < PROTONMAX; ++iproton) {
	proton_pt[iproton] = -999.;
	proton_px[iproton] = -999.;
	proton_py[iproton] = -999.;
	proton_pz[iproton] = -999.;
	proton_energy[iproton] = -999.;
      }
      

      if ( is_good_event(evt) ) {

	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){
	  // Looking for muons
	  //mu-
	  if ( isfinal(*p) && abs((*p)->pdg_id()) == 13  && (*p)->momentum().perp() > 10. && fabs((*p)->momentum().eta()) <= 2.4) {

	    int mupx = (*p)->momentum().px();
	    int mupy = (*p)->momentum().py();
	    int mupz = (*p)->momentum().pz();
	    int mupt = (*p)->momentum().perp();
	    int E    = (*p)->momentum().e();
	    int M    = (*p)->momentum().m(); 
		       
	    mu_px[n_mu] = mupx;
	    mu_pt[n_mu] = mupt;
	    mu_py[n_mu] = mupy;
	    mu_pz[n_mu] = mupz;
	    mu_eta[n_mu] = (*p)->momentum().eta(); 
	    mu_phi[n_mu] = (*p)->momentum().phi(); 
	    mu_energy[n_mu] = E;
	    mu_mass[n_mu] = M;

	    h_mu_pt->Fill(mupt);
	    h_mu_px->Fill(mupx);
	    h_mu_py->Fill(mupy);
	    h_mu_pz->Fill(mupz);
	    h_mu_eta->Fill(mu_eta[n_mu] );
	    h_mu_phi->Fill(mu_phi[n_mu] );
	    h_mu_energy->Fill(mu_energy[n_mu] );
	    h_mu_mass->Fill(mu_mass[n_mu] );
	    ++n_mu;
	  }

	  //neutrinos muons
	  if ( isfinal(*p) && abs((*p)->pdg_id()) == 14 ) {
	    int vmupx = (*p)->momentum().px();
	    int vmupy = (*p)->momentum().py();
	    int vmupt = (*p)->momentum().perp();
	    int vmuE  =  (*p)->momentum().e();
	    int vmuM  = (*p)->momentum().m(); 
	    // vmu_pt[n_vmu] = vmupt;
		
	    vmu_px[n_vmu] = vmupx;
	    vmu_pt[n_vmu] = vmupt;
	    vmu_py[n_vmu] = vmupy;
	    vmu_pz[n_vmu] = (*p)->momentum().pz();
	    vmu_eta[n_vmu] = (*p)->momentum().eta(); 
	    vmu_phi[n_vmu] = (*p)->momentum().phi(); 
	    vmu_energy[n_vmu] = vmuE; 
	    vmu_mass[n_vmu] = vmuM; 

	    h_vmu_pt->Fill((*p)->momentum().perp() );
	    h_vmu_px->Fill( vmu_px[n_vmu] );
	    h_vmu_py->Fill( vmu_py[n_vmu] );
	    h_vmu_pz->Fill( vmu_pz[n_vmu] );
	    h_vmu_eta->Fill( vmu_eta[n_vmu] );
	    h_vmu_phi->Fill( vmu_phi[n_vmu] );
	    h_vmu_energy->Fill(vmu_energy[n_vmu] );
	    h_vmu_mass->Fill( vmu_mass[n_vmu] );
	    ++n_vmu;
	  }

	  // Looking for electrons
	  if ( isfinal(*p) && abs( (*p)->pdg_id() ) == 11  && (*p)->momentum().perp() > 10. && fabs( (*p)->momentum().eta() ) <= 2.5) {

	    int px =(*p)->momentum().px();
	    int py = (*p)->momentum().py();
	    int pz = (*p)->momentum().pz();
	    e_pt[n_e] = (*p)->momentum().perp();
	    e_px[n_e] = px;
	    e_py[n_e] = py;
	    e_pz[n_e] = pz;
	    e_eta[n_e] = (*p)->momentum().eta(); 
	    e_phi[n_e] = (*p)->momentum().phi(); 
	    e_energy[n_e] = (*p)->momentum().e();   
	    e_mass[n_e] = (*p)->momentum().m();                      
	    
	    h_e_pt->Fill( e_pt[n_e]);
	    h_e_px->Fill( e_px[n_e] );
	    h_e_py->Fill( e_py[n_e] );
	    h_e_pz->Fill( e_pz[n_e] );
	    h_e_eta->Fill( e_eta[n_e] );
	    h_e_phi->Fill( e_phi[n_e] );
	    h_e_energy->Fill( e_pt[n_e] );
	    h_e_mass->Fill( e_pt[n_e] );
	    ++n_e;
	  }

	  //neutrinos eletrons
	  if ( isfinal(*p) && abs( (*p)->pdg_id() ) == 12 ) {
	    int vepx =(*p)->momentum().px();
	    int vepy = (*p)->momentum().py();
	    ve_pt[n_ve] = (*p)->momentum().perp();
	    ve_px[n_ve] = vepx;
	    ve_py[n_ve] = vepy;
	    ve_pz[n_ve] = (*p)->momentum().pz();
	    ve_eta[n_ve] = (*p)->momentum().eta(); 
	    ve_phi[n_ve] = (*p)->momentum().phi(); 
	    ve_energy[n_ve] = (*p)->momentum().e(); 
	    ve_mass[n_ve] = (*p)->momentum().m(); 

	    h_ve_pt->Fill( ve_pt[n_ve] );
	    h_ve_px->Fill( ve_px[n_ve] );
	    h_ve_py->Fill( ve_py[n_ve] );
	    h_ve_pz->Fill( ve_pz[n_ve] );
	    h_ve_eta->Fill( ve_eta[n_ve] );
	    h_ve_phi->Fill( ve_phi[n_ve] );
	    h_ve_energy->Fill(ve_energy[n_ve] );
	    h_ve_mass->Fill( ve_mass[n_ve] );
	    ++n_ve;
	  }

	  if ( abs( (*p)->pdg_id() ) == 24 && (*p)->production_vertex() ) {

	     if( (*p)->end_vertex() ){
		HepMC::GenVertex::particle_iterator first_children = (*p)->end_vertex()->particles_begin(HepMC::children);
		HepMC::GenVertex::particle_iterator end_children = (*p)->end_vertex()->particles_begin(HepMC::children);

		HepMC::GenParticle* leptonFromW = 0;
		HepMC::GenParticle* neutrinoFromW = 0;

		int pid = (*p)->pdg_id(); 
		for(HepMC::GenVertex::particle_iterator it_dau = first_children; it_dau != end_children; ++it_dau) {
		   if( abs( (*it_dau)->pdg_id() == 12 ) ||
		       abs( (*it_dau)->pdg_id() == 14 ) || 
		       abs( (*it_dau)->pdg_id() == 16 ) ) neutrinoFromW = (*it_dau);

		   if( abs( (*it_dau)->pdg_id() == 11 ) ||
		       abs( (*it_dau)->pdg_id() == 13 ) || 
		       abs( (*it_dau)->pdg_id() == 15 ) ) leptonFromW = (*it_dau);
		}
		if( leptonFromW && neutrinoFromW) {

		   if( pid == 24 ){  
		      int w1px =(*p)->momentum().px();
		      int w1py =(*p)->momentum().py();
		      int w1pz =(*p)->momentum().pz();
		      int E    =(*p)->momentum().e();
		      w1_mass[n_w1] = (*p)->momentum().m();
		      w1_phi[n_w1] = (*p)->momentum().phi();
		      w1_cos[n_w1] = (*p)->momentum().theta();
		      w1_eta[n_w1] = (*p)->momentum().eta();
		      //int w1mt= w1_transversemass[n_w1];

		      w1_px[n_w1] = w1px;
		      w1_py[n_w1] = w1py;
		      w1_pz[n_w1] = w1pz;
		      w1_energy[n_w1] = E;

		      w1_lep_id[n_w1] = leptonFromW->pdg_id();
		      w1_lep_px[n_w1] = leptonFromW->momentum().px();
		      w1_lep_py[n_w1] = leptonFromW->momentum().py();
		      w1_lep_pz[n_w1] = leptonFromW->momentum().pz();
		      w1_lep_energy[n_w1] = leptonFromW->momentum().e();

		      w1_nu_id[n_w1] = neutrinoFromW->pdg_id();
		      w1_nu_px[n_w1] = neutrinoFromW->momentum().px();
		      w1_nu_py[n_w1] = neutrinoFromW->momentum().py();
		      w1_nu_pz[n_w1] = neutrinoFromW->momentum().pz();
		      w1_nu_energy[n_w1] = neutrinoFromW->momentum().e();

		      //W1 =W+
		      CLHEP::HepLorentzVector w1p4(w1px,w1py,w1pz,E);
		      double w1mt = (w1p4.mt());
		      double w1pt = (w1p4.perp());
		      double w1e  = w1p4.e();

                      
		      h_w1_pt->Fill(w1pt);
		      h_w1_energy->Fill(w1e );
		      h_w1_px->Fill( w1p4.px());
		      h_w1_py->Fill( w1p4.py());
		      h_w1_pz->Fill(w1p4.pz());
		      h_w1_phi->Fill( w1_phi[n_w1] );
		      h_w1_cos->Fill( w1_cos[n_w1] );
		      h_w1_eta->Fill( w1_eta[n_w1] );
		      h_w1_transversemass->Fill(w1mt);
		      h_w1_mass->Fill(w1_mass[n_w1]);

		      ++n_w1;
		   } else {
		      int w2px =(*p)->momentum().px();
		      int w2py =(*p)->momentum().py();
		      int w2pz =(*p)->momentum().pz();
		      int E    =(*p)->momentum().e();

		      w2_px[n_w2] = w2px;
		      w2_py[n_w2] = w2py;
		      w2_pz[n_w2] = w2pz;
		      w2_energy[n_w2] = E;

		      w2_lep_id[n_w2] = leptonFromW->pdg_id();
		      w2_lep_px[n_w2] = leptonFromW->momentum().px();
		      w2_lep_py[n_w2] = leptonFromW->momentum().py();
		      w2_lep_pz[n_w2] = leptonFromW->momentum().pz();
		      w2_lep_energy[n_w2] = leptonFromW->momentum().e();

		      w2_nu_id[n_w2] = neutrinoFromW->pdg_id();
		      w2_nu_px[n_w2] = neutrinoFromW->momentum().px();
		      w2_nu_py[n_w2] = neutrinoFromW->momentum().py();
		      w2_nu_pz[n_w2] = neutrinoFromW->momentum().pz();
		      w2_nu_energy[n_w2] = neutrinoFromW->momentum().e();

		      //W2 =W-
		      CLHEP::HepLorentzVector w2p4(w2px,w2py,w2pz,E);
		      double w2mt = (w2p4.mt());
		      double w2pt = (w2p4.perp());
		      double w2e  = w2p4.e();

		      w2_pt[n_w2] = (*p)->momentum().perp();
		      h_w2_energy->Fill(w2e );
		      w2_mass[n_w2] = (*p)->momentum().m();
		      w2_phi[n_w2] = (*p)->momentum().phi();
		      w2_cos[n_w2] = (*p)->momentum().theta();
		      w2_eta[n_w2] = (*p)->momentum().eta();

		      h_w2_px->Fill( w2p4.px());
		      h_w2_py->Fill( w2p4.py());
		      h_w2_pz->Fill(w2p4.pz());
		      h_w2_pt->Fill( w2pt);
		      h_w2_phi->Fill( w2_phi[n_w2] );
		      h_w2_cos->Fill( w2_cos[n_w2] );
		      h_w2_eta->Fill( w2_eta[n_w2] );
		      h_w2_transversemass->Fill(w2mt);
		      h_w2_mass->Fill(w2_mass[n_w2]);
		      ++n_w2;
		   }
		}  
	     }
	  }                     

	  if ( abs((*p)->pdg_id() == 22) && (*p)->status() == -1 ){

	    int apx =(*p)->momentum().px();
	    int apy =(*p)->momentum().py();
	    int apz =(*p)->momentum().pz();
	    int E    =(*p)->momentum().e();

	    a_px[n_a] = apx;
	    a_py[n_a] = apy;
	    a_pz[n_a] = apz;
	    a_energy[n_a] = E;

	    CLHEP::HepLorentzVector ap4(apx,apy,apz,E);
	    double apt = (ap4.perp());
	    double ae = ap4.e();
	    a_pt[n_a] = (*p)->momentum().perp();
	    a_mass[n_a] = (*p)->momentum().m();

	    h_a_px->Fill( ap4.px());
	    h_a_py->Fill( ap4.py());
	    h_a_pz->Fill(ap4.pz());
	    h_a_fracmom->Fill(ap4.pz()/4000);
	    h_a_pini->Fill(1-ap4.pz()/4000);
	    h_a_pt->Fill( apt);
	    h_a_energy->Fill(ae);
	    h_a_mass->Fill(a_mass[n_a]);
	    ++n_a;
	  }
	  
	  if ( abs((*p)->pdg_id() == 2212) && (*p)->status() == 1 ){

	    proton_px[n_proton] = (*p)->momentum().px();
	    proton_py[n_proton] = (*p)->momentum().py();
	    proton_pz[n_proton] = (*p)->momentum().pz();
	    proton_pt[n_proton] = (*p)->momentum().perp();
	    proton_energy[n_proton] =(*p)->momentum().e();
    	    h_proton_px->Fill( proton_px[n_proton]);
	    h_proton_py->Fill( proton_py[n_proton]);
	    h_proton_pz->Fill(proton_pz[n_proton]);
	    h_proton_pt->Fill( proton_pt[n_proton]);
	    h_proton_energy->Fill(proton_energy[n_proton]);
	    ++n_proton;
	}



	  // Looking for charged particles
	  if ( isfinal(*p) ){
	    int pdg_id = (*p)->pdg_id();
	    // Could potentially be slow.
	    // Instead build map with charge per PDG id before event loop.
	    HepPDT::ParticleData * pd;
	    pd = pdt.particle( HepPDT::ParticleID( pdg_id ) );
	    double charge = pd->charge(); 
	    if( charge != 0. &&
		fabs( (*p)->momentum().eta() ) <= 2.5 &&
		(*p)->momentum().perp() >= 0.01) {
	      chg_id[n_chg] = pdg_id;
	      chg_ch[n_chg] = charge;
	      chg_pt[n_chg] = (*p)->momentum().perp();
	      chg_px[n_chg] = (*p)->momentum().px();
	      chg_py[n_chg] = (*p)->momentum().py();
	      chg_pz[n_chg] = (*p)->momentum().pz();
	      chg_eta[n_chg] = (*p)->momentum().eta(); 
	      chg_phi[n_chg] = (*p)->momentum().phi(); 
	      chg_energy[n_chg] = (*p)->momentum().e(); 
	      chg_mass[n_chg] = (*p)->momentum().m(); 

	      h_chg_id->Fill( chg_id[n_chg] );
	      h_chg_ch->Fill( chg_ch[n_chg] );
	      h_chg_pt->Fill( chg_pt[n_chg]);
	      h_chg_px->Fill( chg_px[n_chg]);
	      h_chg_py->Fill( chg_py[n_chg]);
	      h_chg_pz->Fill( chg_pz[n_chg]);
	      h_chg_eta->Fill( chg_eta[n_chg] );
	      h_chg_phi->Fill( chg_phi[n_chg] );
	      h_chg_energy->Fill( chg_energy[n_chg]);
	      h_chg_mass->Fill( chg_mass[n_chg]);
	      ++n_chg;
	    }
	  }
	}

	h_nmu->Fill( n_mu );
	h_nvmu->Fill( n_vmu );
	h_ne->Fill( n_e );
	h_nw1->Fill( n_w1 );
	h_nw2->Fill( n_w2 );
	h_nve->Fill( n_ve );
	h_na->Fill( n_a );
	h_nproton->Fill(n_proton);
	h_nchg->Fill( n_chg );

	// Here it will fill TTree only for good events
	T->Fill();
	ascii_out << evt;
	++num_good_events;
      }

      // Here it will fill TTree for all events 
      //T->Fill();
      delete evt;
      ascii_in >> evt;
    }
    //........................................PRINT RESULT
    std::cout << num_good_events << " out of " << icount 
	      << " processed events passed the cuts. Finished." << std::endl;
    T->Print();
    // Output file
    TFile* output = new TFile("aa_wwa0wp150hepmc_13tev.root","RECREATE");
    output->cd();
    // Write TTree and histograms to file
    T->Write();
    h_nmu->Write();
    h_mu_pt->Write();
    h_mu_px->Write();
    h_mu_py->Write();
    h_mu_pz->Write();
    h_mu_eta->Write();
    h_mu_phi->Write();
    h_mu_energy->Write();
    h_mu_mass->Write();

    h_na->Write();
    h_a_pt->Write();
    h_a_px->Write();
    h_a_py->Write();
    h_a_pz->Write();
    h_a_pini->Write();
    h_a_fracmom->Write();
    h_a_energy->Write();
    h_a_mass->Write();


    h_nproton->Write();
    h_proton_px->Write();
    h_proton_py->Write();
    h_proton_pz->Write();
    h_proton_pt->Write();
    h_proton_energy->Write();
    

    h_nvmu->Write();
    h_vmu_pt->Write();
    h_vmu_px->Write();
    h_vmu_py->Write();
    h_vmu_pz->Write();
    h_vmu_eta->Write();
    h_vmu_phi->Write();
    h_vmu_energy->Write();
    h_vmu_mass->Write();

    h_nve->Write();
    h_ve_pt->Write();
    h_ve_px->Write();
    h_ve_py->Write();
    h_ve_pz->Write();
    h_ve_eta->Write();
    h_ve_phi->Write();
    h_ve_energy->Write();
    h_ve_mass->Write();

    h_ne->Write();
    h_e_pt->Write();
    h_e_px->Write();
    h_e_py->Write();
    h_e_pz->Write();
    h_e_eta->Write();
    h_e_phi->Write();
    h_e_energy->Write();
    h_e_mass->Write();

    h_nchg->Write();
    h_chg_id->Write();
    h_chg_ch->Write();
    h_chg_pt->Write();
    h_chg_px->Write();
    h_chg_py->Write();
    h_chg_pz->Write();
    h_chg_eta->Write();
    h_chg_phi->Write();
    h_chg_energy->Write();
    h_chg_mass->Write();

    h_nw1->Write();
    h_w1_pt->Write();
    h_w1_px->Write();
    h_w1_py->Write();
    h_w1_pz->Write();
    h_w1_mass->Write();
    h_w1_energy->Write();
    h_w1_transversemass->Write();
    h_w1_phi->Write();
    h_w1_cos->Write();
    h_w1_eta->Write();

    h_nw2->Write();
    h_w2_pt->Write();
    h_w2_px->Write();
    h_w2_py->Write();
    h_w2_pz->Write();
    h_w2_energy->Write();
    h_w2_mass->Write();
    h_w2_transversemass->Write();
    h_w2_phi->Write();
    h_w2_cos->Write();
    h_w2_eta->Write();
    h_2lep2neu->Write();

    output->Close();

  return 0;
} //fim do int main



