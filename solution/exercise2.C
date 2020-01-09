#include <iostream>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

void exercise2(void) {
  //==============================================================
  // Unfolding exercise 2
  //  study correlation of rRec and rGen
  //  determine the RMS resolution using a profile plot
  //==============================================================

  // the event samples are divided into subsets
  //  for most exercises we divide the event samples into 30 subsets
  //  and use only the first subset
  int nSubset=30;  int iSubset=0;

  // calculate proper uncertainties when filling weighted events
  //   into histograms
  TH1::SetDefaultSumw2();

  // change to memory
  gDirectory->cd("Rint:/");
  TDirectory *histogramDir=gDirectory;

  // book histograms here
  TH2D *hist_mc1_recgen=new TH2D("hist_mc1_recgen","MC1 truth to rec migrations;r_{gen};r_{rec}",
				 25,0.,1.,25,0.,1.);
  TProfile *hist_mc1_resolution=new TProfile
    ("hist_mc1_resolution","MC1 resolution;r_{gen};r_{rec}-r_{gen}",50,0.,1.,"s");
  TH2D *hist_mc2_recgen=new TH2D("hist_mc2_recgen","MC2 truth to rec migrations;r_{gen};r_{rec}",
				 25,0.,1.,25,0.,1.);
  TProfile *hist_mc2_resolution=new TProfile
    ("hist_mc2_resolution","MC2 resolution;r_{gen};r_{rec}-r_{gen}",50,0.,1.,"s");

  // loop over mc1 events
  {
    // open input file and  get access to the events
    TFile inputFile("mc1.root");
    TTree *tree;
    inputFile.GetObject("recgen",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString recWeight="w*isTrig";
    TString bgrWeight="w*isTrig*isBgr";
    TString signalWeight="w*isTrig*(1-isBgr)";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc1_recgen","rRec:rGen",signalWeight,
		  "",nEvent,firstEvent);
    tree->Project("hist_mc1_resolution","rRec-rGen:rGen",signalWeight,
		  "profs",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // loop over mc2 events
  {
    // open input file and  get access to the events
    TFile inputFile("mc2.root");
    TTree *tree;
    inputFile.GetObject("recgen",tree);
    // determine number of events to use  
    int firstEvent=iSubset*tree->GetEntriesFast() / nSubset;
    int lastEvent=(iSubset+1)*tree->GetEntriesFast() / nSubset;
    int nEvent=lastEvent-firstEvent;

    // define cut
    TString recWeight="w*isTrig";
    TString bgrWeight="w*isTrig*isBgr";
    TString signalWeight="w*isTrig*(1-isBgr)";
    TString genWeight="w*(1-isBgr)";

    // fill histograms
    histogramDir->cd();
    tree->Project("hist_mc2_recgen","rRec:rGen",signalWeight,
		  "",nEvent,firstEvent);
    tree->Project("hist_mc2_resolution","rRec-rGen:rGen",signalWeight,
		  "profs",nEvent,firstEvent);

    // clean up
    delete tree;
  }

  // make plot
  TCanvas *canvas=
    new TCanvas("resolution studies","rRec and rGen",600,600);
  canvas->Divide(2,2);

  TVirtualPad *p1=canvas->cd(1);
  p1->SetRightMargin(0.15);
  hist_mc1_recgen->Draw("COLZ");
  hist_mc1_recgen->Draw("BOX SAME");

  canvas->cd(2);
  hist_mc1_resolution->GetYaxis()->SetRangeUser(-0.2,0.3);
  hist_mc1_resolution->SetLineColor(kBlue);
  hist_mc1_resolution->SetFillColor(kBlue-10);
  hist_mc1_resolution->SetFillStyle(1001);
  TH1 *h1=hist_mc1_resolution->DrawCopy("E2");
  hist_mc1_resolution->SetFillStyle(0);
  hist_mc1_resolution->Draw("HIST SAME");
  TLegend *legend1=new TLegend(0.35,0.7,0.9,0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->AddEntry(h1,"MC1 mean and RMS","lf");
  legend1->Draw();

  TVirtualPad *p3=canvas->cd(3);
  p3->SetRightMargin(0.15);
  hist_mc2_recgen->Draw("COLZ");
  hist_mc2_recgen->Draw("BOX SAME");

  canvas->cd(4);
  hist_mc2_resolution->GetYaxis()->SetRangeUser(-0.2,0.3);
  hist_mc2_resolution->SetLineColor(kBlue);
  hist_mc2_resolution->SetFillColor(kBlue-10);
  hist_mc2_resolution->SetFillStyle(1001);
  TH1 *h2=hist_mc2_resolution->DrawCopy("E2");
  hist_mc2_resolution->SetFillStyle(0);
  hist_mc2_resolution->Draw("HIST SAME");
  TLegend *legend2=new TLegend(0.35,0.7,0.9,0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(h2,"MC2 mean and RMS","lf");
  legend2->Draw();

}
