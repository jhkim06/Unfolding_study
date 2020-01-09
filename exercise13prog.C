// system include files
#include <iostream>

// from root
#include <TError.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

// from TUnfold V17.3
#include "TUnfoldDensity.h"
#include "TUnfoldBinningXML.h"

using namespace std;

int main() {
  // report warnings and stop on errors
  gErrorIgnoreLevel=kInfo+1;
  gErrorAbortLevel=kError;

  // open root file for output
  TFile *output=new TFile("exercise13.root","RECREATE");

  // read binning schemes in XML format
  TUnfoldBinning *detectorBinning,*generatorBinning;
  TDOMParser parser;
  Int_t error=parser.ParseFile("exercise13binning.xml");
  if(error) cout<<"error="<<error<<" from TDOMParser\n";
  TXMLDocument const *XMLdocument=parser.GetXMLDocument();
  detectorBinning=TUnfoldBinningXML::ImportXML(XMLdocument,"detector");
  generatorBinning=TUnfoldBinningXML::ImportXML(XMLdocument,"generator");

  detectorBinning->Write();
  generatorBinning->Write();

  // within generatorBinnig, locate signal and background
  TUnfoldBinning const *signalBinning=generatorBinning->FindNode("signal");
  TUnfoldBinning const *backgroundBinning=generatorBinning->FindNode("background");

  // book histograms
  // Data: 1d histogram designed for use with TUnfold
  TH1 *hist_data_recTUnfold=
    detectorBinning->CreateHistogram("hist_data_recTUnfold");
  // Data signal truth: 2d histogram
  TH1 *hist_data_truth2d=
    signalBinning->CreateHistogram("hist_data_truth2d",kTRUE);

  // MC: 2d histogram of migrations desined for use with TUnfold
  TH2 *hist_mc1_recgenTUnfold=
    TUnfoldBinning::CreateHistogramOfMigrations
    (generatorBinning,detectorBinning,"hist_mc1_recgen");

  // loop over data events
  {
    // open input file and  get access to the events
    TFile inputFile("data.root");
    TTree *tree;
    inputFile.GetObject("rec",tree);

    float rRec,pRec,w;
    int isTrig;
    tree->SetBranchAddress("rRec",&rRec);
    tree->SetBranchAddress("pRec",&pRec);
    tree->SetBranchAddress("w",&w);
    tree->SetBranchAddress("isTrig",&isTrig);

    output->cd();

    for(int iEvent=0;iEvent<tree->GetEntriesFast();iEvent++) {
      tree->GetEvent(iEvent);
      // fill histograms
      if(isTrig==1) {
	double aRec=rRec*TMath::Sin(pRec);
	double bRec=rRec*TMath::Cos(pRec);
	int iBinRec=detectorBinning->GetGlobalBinNumber(aRec,bRec);
	hist_data_recTUnfold->Fill(iBinRec,w);
      }
    }
    // clean up
    delete tree;
  }

  // loop over data truth events and fill histogram
  {
    TFile inputFile("datatruth.root");
    TTree *tree;
    inputFile.GetObject("gen",tree);

    // define cut
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    output->cd();
    tree->Project("hist_data_truth2d","rGen*cos(pGen):rGen*sin(pGen)",
		  genWeight);

    // clean up
    delete tree;
  }

    // loop over mc1 events
  {
    // open input file and  get access to the events
    TFile inputFile("mc1.root");
    TTree *tree;
    inputFile.GetObject("recgen",tree);

    float rRec,pRec,w;
    int isTrig;

    float rGen,pGen;
    int isBgr;

    tree->SetBranchAddress("rRec",&rRec);
    tree->SetBranchAddress("pRec",&pRec);
    tree->SetBranchAddress("w",&w);
    tree->SetBranchAddress("isTrig",&isTrig);

    tree->SetBranchAddress("rGen",&rGen);
    tree->SetBranchAddress("pGen",&pGen);
    tree->SetBranchAddress("isBgr",&isBgr);

    output->cd();

    for(int iEvent=0;iEvent<tree->GetEntriesFast();iEvent++) {
      tree->GetEvent(iEvent);
      // fill histograms
      // calculate variables
      double aRec=-2;
      double bRec=-2;
      double aGen=-1;
      double bGen=-2;
      if(isTrig==1) {
	aRec=rRec*TMath::Sin(pRec);
	bRec=rRec*TMath::Cos(pRec);
      }
      if(isBgr==0) {
	aGen=rGen*TMath::Sin(pGen);
	bGen=rGen*TMath::Cos(pGen);
      }

      // (1) locate rec bin
      int iBinRec=0;
      if(isTrig==1) {
	iBinRec=detectorBinning->GetGlobalBinNumber(aRec,bRec);
      }
      // (2) locate gen bin
      int iBinGen=0;
      if(isBgr==1) {
	// background event
	iBinGen=backgroundBinning->GetGlobalBinNumber(aRec,bRec);
      } else {
	// signal event
	iBinGen=signalBinning->GetGlobalBinNumber(aGen,bGen);
      }
      // fill histogram
      hist_mc1_recgenTUnfold->Fill(iBinGen,iBinRec,w);
    }
    // clean up
    delete tree;
  }

  // unfold
  TUnfoldDensity unfold(hist_mc1_recgenTUnfold,
			TUnfold::kHistMapOutputHoriz,
			TUnfold::kRegModeSize,
			TUnfold::kEConstraintArea,
			TUnfoldDensity::kDensityModeNone,
			generatorBinning,
			detectorBinning);
  unfold.SetInput(hist_data_recTUnfold,1.,1.);
  //TSpline *rhoVSlogtau=0;
  //unfold.ScanTau(30,0.,0.,&rhoVSlogtau);
  TGraph *lCurve=0;
  unfold.ScanLcurve(30,0.,0.,&lCurve);
  unfold.GetOutput("hist_data_unfoldSignal",
		   "data (TUnfolded);a_{TUNFOLD};b_{TUNFOLD}","signal");
  unfold.GetOutput("hist_data_unfoldBackground",0,"background");

  // close output file
  output->cd();
  output->Write();
  delete output;

  return 0;
}
