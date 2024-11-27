#include "fastjet/Selector.hh" //.......... Background Sutraction event by event
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"//.......... Background Sutraction event by event
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/Recluster.hh"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"  
#include <ctime>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <valarray>


using namespace Pythia8;

double Calculate_pX( double pT, double eta, double phi){
  return(pT*TMath::Cos(phi));
}

double Calculate_pY( double pT, double eta, double phi){
  return(pT*TMath::Sin(phi));
}

double Calculate_pZ( double pT, double eta, double phi){
  return( pT*TMath::SinH(eta) );
}

double Calculate_E( double pT, double eta, double phi){
  double pZ = Calculate_pZ(pT,eta, phi);
  return( TMath::Sqrt(pT*pT + pZ*pZ) );
}

//__________________________________________________________________________

bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax) {
  if(fjJet.eta() > etaMax || fjJet.eta() < etaMin){
    return false;
  } else {
    return true;
  }
}
//_________________________________________________________________________ 
Double_t RelativePhi(Double_t mphi,Double_t vphi){
  //Get relative azimuthal angle of two particles -pi to pi
  if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
  else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

  if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
  else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

  Double_t dphi = mphi - vphi;
  if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
  else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

  return dphi;//dphi in [-Pi, Pi]
}
//__________________________________________________________________________
float fShapesVar[7];
TTree *fTreeResponse = new TTree("variables", "variables");


int CheckSDProng (fastjet::PseudoJet prong, double jetR) {
  int flagsd=0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

  try {
    std::vector<fastjet::PseudoJet> prongParticles = prong.constituents();
    fastjet::ClusterSequence prongSeq(prongParticles, jet_def);

    std::vector<fastjet::PseudoJet> prongJets = prongSeq.inclusive_jets(0);
    prongJets = sorted_by_pt(prongJets);

    if (prongJets.empty()) {
      return 0;
    } else {
      fastjet::PseudoJet jj = prongJets[0];
      fastjet::PseudoJet j1;
      fastjet::PseudoJet j2;

      while (jj.has_parents(j1, j2)){
	      if (j1.perp() < j2.perp()) std::swap(j1, j2);

	      double delta_R = j1.delta_R(j2);
	      double z = j2.perp() / (j1.perp() + j2.perp());
	      double delta = delta_R / jetR;

	      if(z>0.1*pow(delta,2)) flagsd=1;

	      jj = j1;
      }
      vector < fastjet::PseudoJet > constitjj = sorted_by_pt(jj.constituents());
      if (flagsd==0 && constitjj[0].user_index()>0) {
	      return constitjj[0].user_index();
      } else {
	      return 0;
      }
    }
  }
  catch (fastjet::Error) { /*return -1;*/ }
}



void IterativeDeclustering (fastjet::PseudoJet jet, double jetR, double weight, double flavour, int runparton){
   
  double zg=0;
  double rg=0;
  double ktg=0;
  double nsd=0;
  int flagsd=0;
  int flagl=-1;
  int flagsl=-1;
  int tagpart=-1; 
  if(flavour>0) tagpart=2;

  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

  try{
    std::vector<fastjet::PseudoJet> particles = jet.constituents();
    fastjet::ClusterSequence cs_gen(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = cs_gen.inclusive_jets(1);
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;  // subjet 1 (largest pt)
    fastjet::PseudoJet j2;  // subjet 2 (smaller pt)

    // Unclustering jet
    while(jj.has_parents(j1,j2)){       
      if(j1.perp() < j2.perp()) std::swap(j1,j2);
      // Calculate deltaR and Zg between j1 and j2
      double delta_R = j1.delta_R(j2);
      double cut=1;
      double angkt=j2.perp();
      double z = j2.perp()/(j2.perp()+j1.perp());
      double kt= j2.perp()*delta_R;
      double erad=j2.e()+j1.e();

      if(z > 0.1 && flagsd==0) {
	      zg = z;
	      rg = delta_R;
	      ktg=kt;
	      flagsd=1;
        flagl=CheckSDProng(j1,jetR/2);
	      flagsl=CheckSDProng(j2,jetR/2);
      }

      jj=j1;
    }

    fShapesVar[1]=zg;
    fShapesVar[2]=ktg;
    fShapesVar[3]=rg;
    fShapesVar[4]=weight;
    fShapesVar[5]=flagl;
    fShapesVar[6]=flagsl;
  }
  catch (fastjet::Error) { /*return -1;*/ }
}


//__________________________________________________________________________


int main(int argc, char* argv[]) {
  Int_t cislo = -1;                 //unique number for each file
  Int_t tune  = -1;                 //pythia tune

  Float_t pthatmin=5;
  Float_t pthatmax=5020; 
  Int_t flavour=2;
  Float_t jetR=0.4;
 
  if(argc!=7){  
    cout<<"Usage:"<<endl<<"./pygen <PythiaTune> <Number> <nEvts> <flavour> <jetR>"<<endl;
    return 0;
  }
  tune  = atoi(argv[1]);
  cislo = atoi(argv[2]); // number of random seed

  cout<<cislo<<endl;

  Int_t nEvents= atoi(argv[3]);   //(Int_t) 1e3 + 1.0;
 
  flavour=atoi(argv[4]);
  jetR=atof(argv[5]);

  Int_t outputnumber = atoi(argv[6]);
  cout<<nEvents<<" events"<<endl;

  const Int_t nVar = 7;
  fTreeResponse = new TTree("variables", "variables");
  TString *fShapesVarNames = new TString [nVar];

  fShapesVarNames[0] = "ptJet";
  fShapesVarNames[1] = "zg";
  fShapesVarNames[2] = "kg";
  fShapesVarNames[3] = "rg";
  fShapesVarNames[4] = "weight";
  fShapesVarNames[5] = "flagl";
  fShapesVarNames[6] = "flagsl";

  for (Int_t ivar=0; ivar < nVar; ivar++){  
    fTreeResponse->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar],Form("%s/F", fShapesVarNames[ivar].Data()));
  }

  //__________________________________________________________________________
  //ANALYSIS SETTINGS

  // double jetParameterR   = 0.4; //jet R
  double trackEtaCut     = 2.4;
  double trackLowPtCut   = 0.; //GeV
  
  TRandom3 *r3=new TRandom3();
  //__________________________________________________________________________
  //PYTHIA SETTINGS
  Int_t mypart=-1;
  Int_t mypartq=-1;
  TString name;
  
  //Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:idA = 2212"); //beam 1 proton
  pythia.readString("Beams:idB = 2212"); //beam 2 proton
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("Tune:pp = 5");  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",cislo));
  pythia.readString("HardQCD:all = on");
  pythia.readString("TimeShower:QEDshowerByOther= off");
  
  if(pthatmin<0 || pthatmax <0){
    pythia.readString("PhaseSpace:pTHatMin = 0."); // <<<<<<<<<<<<<<<<<<<<<<<
  } else {
    pythia.readString("PhaseSpace:bias2Selection = on");
    pythia.readString("PhaseSpace:bias2SelectionPow = 5.");
    pythia.readString("PhaseSpace:bias2SelectionRef = 8.");
    pythia.readString(name.Data());
  }

  if(flavour==2) pythia.readString("511:mayDecay = off");
  if(flavour==1) pythia.readString("421:mayDecay = off");
  if(flavour==2) mypart=511;
  if(flavour==1) mypart=22;
  if(flavour==2) mypartq=5;
  if(flavour==1) mypartq=4;

  pythia.init();
  
  //_________________________________________________________________________________________________
  //FASTJET  SETTINGS

  double etamin_Sig = - trackEtaCut + jetR; //signal jet eta range
  double etamax_Sig = - etamin_Sig;
  
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  
  fastjet::JetDefinition *jetDefAKT_Sig = NULL;
  
  jetDefAKT_Sig = new fastjet::JetDefinition(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
  
  fastjet::GhostedAreaSpec ghostareaspec(trackEtaCut, 1, 0.05); //ghost 
  //max rap, repeat, ghostarea default 0.01
  fastjet::AreaType areaType = fastjet::active_area_explicit_ghosts;
  fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, ghostareaspec);
  
  // Fastjet input
  std::vector<fastjet::PseudoJet> fjInputs1;
  std::vector<fastjet::PseudoJet> fjInputs1parton;
  std::vector<fastjet::PseudoJet> fjInputs2;
 
  //Begin event loop
  Int_t count=0;
  for(int i = 0; i < nEvents; i++){
    double fourvec[4];
    double fourvecp[4];
    if(!pythia.next()) continue;
    count=count+1;
    Double_t weight=pythia.info.weight();

    fjInputs1.resize(0);
    fjInputs1parton.resize(0);

    if(pythia.event.size()==0) continue; 

    for(int j = 0; j < pythia.event.size(); j++){
      if(pythia.event[j].isFinal()){
        if(pythia.event[j].pT() < trackLowPtCut) continue;                 //pt cut
        if(TMath::Abs(pythia.event[j].eta()) > trackEtaCut) continue;      //eta cut
	      fourvec[0]=pythia.event[j].px();
        fourvec[1]=pythia.event[j].py();
        fourvec[2]=pythia.event[j].pz();
        fourvec[3]=pythia.event[j].e();
        fastjet::PseudoJet PythiaParticle(fourvec);
        PythiaParticle.set_user_index(0);
	      if(flavour>0){
	        if (TMath::Abs(pythia.event[j].id())==mypart) {
	          PythiaParticle.set_user_index(2);
		        int MotherIndex = pythia.event[j].mother1();
            int MotherType = TMath::Abs(pythia.event[MotherIndex].id());

            while (MotherIndex > 0) {
              MotherType = TMath::Abs(pythia.event[MotherIndex].id());
              if (MotherType == 221) {  //decay
                PythiaParticle.set_user_index(11);
                break;
              } else if (MotherType == 111) {  //decay
                PythiaParticle.set_user_index(12);
                break;
              } else if (MotherType >= 1 && MotherType <= 5) {  //quark
                PythiaParticle.set_user_index(10);
                break;
              } else if (MotherType == 22) {
                MotherIndex = pythia.event[MotherIndex].mother1();
              } else {
                break;
              }
            }
          }
	      }
        fjInputs1.push_back(PythiaParticle);
	    }
    }//event size


    try{
      std::vector<fastjet::PseudoJet> PythiaJets;
      fastjet::ClusterSequence clustSeq_Sig(fjInputs1, *jetDefAKT_Sig);
      PythiaJets = sorted_by_pt(clustSeq_Sig.inclusive_jets(1.));
      
      for(Int_t n=0;n<PythiaJets.size();n++){
	      int flagtag = 0;
	      vector < fastjet::PseudoJet > constit = sorted_by_pt(PythiaJets[n].constituents());
        if(!EtaCut(PythiaJets[n], etamin_Sig, etamax_Sig)) continue;
	      for (int t=0; t < constit.size(); t++){
	        if (constit[t].user_index() != 0){
	          flagtag = 1;
	          break;
	        }
	      }
	      if (flagtag == 0) continue;
        fShapesVar[0]=PythiaJets[n].perp();
        
	      IterativeDeclustering(PythiaJets[n], jetR, weight, flavour, 0);

	      fTreeResponse->Fill();
      }
    }
    catch (fastjet::Error) { /*return -1;*/ }
  }
  //End event loop
  TFile* outFile = new TFile(Form("Output_simple_num%d_flavour%d_jetR%f_nEvents%d_lundplane_allsplits_allproc_sync%d.root",cislo,flavour,jetR,nEvents,outputnumber),"RECREATE");
  
  outFile->cd();
  fTreeResponse->Write();
  outFile->Close();
  return 0;
}
