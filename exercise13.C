#include <iostream>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TUnfold.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TMarker.h>

using namespace std;

void exercise13(void) {
  //==============================================================
  // Unfolding exercise 13
  //  2D plots of the variables a and b
  //  compare bin-by-bin unfolding and regularized unfolding results
  //==============================================================

  // calculate proper uncertainties when filling weighted events
  //   into histograms
  TH1::SetDefaultSumw2();

  // change to memory
  gDirectory->cd("Rint:/");
  TDirectory *histogramDir=gDirectory;

  // book histograms here
  int nABGen=15;
  double ab0Gen=-0.9;
  double ab1Gen=0.9;
  int nABRec=15;
  double ab0Rec=-0.9;
  double ab1Rec=0.9;
  // data
  TH2D *hist_data_rec2d=new TH2D
    ("hist_data_rec2d","Data (reconstructed);a_{rec};b_{rec}",
     nABRec,ab0Rec,ab1Rec,nABRec,ab0Rec,ab1Rec);
  TH2D *hist_data_truth2d=new TH2D
    ("hist_data_truth2d","Data (truth);a_{gen};b_{gen}",
     nABGen,ab0Gen,ab1Gen,nABGen,ab0Gen,ab1Gen);
  // MC1
  TH2D *hist_mc1_rec2d=new TH2D
    ("hist_mc1_rec2d","MC1(reconstructed);a_{rec};b_{rec}",
     nABRec,ab0Rec,ab1Rec,nABRec,ab0Rec,ab1Rec);
  TH2D *hist_mc1_gen2d=new TH2D
    ("hist_mc1_gen2d","MC1(truth);a_{gen};b_{gen}",
     nABGen,ab0Gen,ab1Gen,nABGen,ab0Gen,ab1Gen);
  // derived histogram
  TH2D *hist_data_BBB2d=new TH2D
    ("hist_data_BBB2d","Data (bin-by-bin);a_{BBB};b_{BBB}",
     nABGen,ab0Gen,ab1Gen,nABGen,ab0Gen,ab1Gen);

  // loop over data events
  {
    // open input file and  get access to the events
    TFile inputFile("data.root");
    TTree *tree;
    inputFile.GetObject("rec",tree);

    // define cut
    TString recWeight="w*isTrig";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_data_rec2d","rRec*cos(pRec):rRec*sin(pRec)",recWeight);

    // clean up
    delete tree;
  }

  // loop over data truth events
  {
    // open input file and  get access to the events
    TFile inputFile("datatruth.root");
    TTree *tree;
    inputFile.GetObject("gen",tree);

    // define cut
    TString genWeight="w*(1-isBgr)";
    TString bgrWeight="w*isBgr";

    // fill histograms
    histogramDir->cd();
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

    // define cut
    TString genWeight="w*(1-isBgr)";
    TString bgrWeight="w*isBgr";
    TString recWeight="w*isTrig";
    TString recbgrWeight="w*isTrig*isBgr";
    TString recsigWeight="w*isTrig*(1-isBgr)";
    TString norecbgrWeight="w*(1-isTrig)*isBgr";
    TString norecsigWeight="w*(1-isTrig)*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_gen2d","rGen*cos(pGen):rGen*sin(pGen)",genWeight);
    tree->Project("hist_mc1_rec2d","rRec*cos(pRec):rRec*sin(pRec)",recWeight);

    // clean up
    delete tree;
  }

  // bin-by-bin unfolding
  hist_data_BBB2d->Divide(hist_mc1_gen2d,hist_mc1_rec2d);
  hist_data_BBB2d->Multiply(hist_data_rec2d);

  // read results of regularized unfolding
  TFile regularizedResults("exercise13.root");
  TH2 *hist_data_unfoldSignal;
  regularizedResults.GetObject("hist_data_unfoldSignal",
				hist_data_unfoldSignal);

  // make plot
  TCanvas *canvas=
    new TCanvas("2d Unfolding rec, truth, bin-by-bin","",900,600);
  canvas->Divide(3,2);

  canvas->cd(1);
  TVirtualPad *p1=canvas->cd(1);
  p1->SetRightMargin(0.15);
  hist_data_rec2d->Draw("COLZ");
  hist_data_rec2d->Draw("BOX SAME");

  canvas->cd(2);
  TVirtualPad *p2=canvas->cd(2);
  p2->SetRightMargin(0.15);
  hist_data_truth2d->SetMinimum(-0.1);
  hist_data_truth2d->Draw("COLZ");
  hist_data_truth2d->Draw("BOX SAME");

  canvas->cd(3);
  TVirtualPad *p3=canvas->cd(3);
  p3->SetRightMargin(0.15);
  hist_data_BBB2d->SetMinimum(-0.1);
  hist_data_BBB2d->Draw("COLZ");
  hist_data_BBB2d->Draw("BOX SAME");

  canvas->cd(4);
  TVirtualPad *p4=canvas->cd(4);
  p4->SetRightMargin(0.15);
  hist_mc1_rec2d->Draw("COLZ");
  hist_mc1_rec2d->Draw("BOX SAME");

  canvas->cd(5);
  TVirtualPad *p5=canvas->cd(5);
  p5->SetRightMargin(0.15);
  hist_mc1_gen2d->SetMinimum(-0.1);
  hist_mc1_gen2d->Draw("COLZ");
  hist_mc1_gen2d->Draw("BOX SAME");

  canvas->cd(6);
  TVirtualPad *p6=canvas->cd(6);
  p6->SetRightMargin(0.15);
  TH1 *hUnfold=hist_data_unfoldSignal->DrawCopy("COLZ");
  hUnfold->Draw("BOX SAME");

}
