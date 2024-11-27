#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TObject.h"
using namespace std;
//#include "tdrstyle.C"
////#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"


void macro(){

    TFile *f = TFile::Open("Output_simple_num1_flavour0_jetR0.400000_lundplane_allsplits_allproc_sync.root");

    TTree *tree = (TTree*)f->Get("variables");
    Long64_t nentries = tree->GetEntries();
    std::cout<< " Entries : " << nentries <<std::endl;

    float rg=0;
    float kt=0;
    tree->SetBranchAddress("rg", &rg);
    tree->SetBranchAddress("kg", &kt);
    
    TH2F* hLund = new TH2F ("hLund", "", 100, 1.5,5, 100, 0,3);

    for (int i=0; i< nentries; i++){
        tree->GetEntry(i);
        if (rg==0) continue;
        double rglog=TMath::Log(1./rg);
        double ktlog=TMath::Log(kt);

        hLund->Fill(rglog,ktlog);

    }
    hLund->Draw("colz");
}
